function autorunList(instrName)
% Automated script for testing model performance

%% Add needed paths

% Get fileparts and path to '../'

splPath = regexp(which('autorunList'),filesep,'split');
fDir    = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];

sDir    = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep,'SimData',filesep];

if exist(sDir,'dir') == 0 
    mkdir(sDir);
end

% Add dependencies
addpath([fDir,'0-SharedTools']);
addpath([fDir,'1-GroundTruth']);
addpath([fDir,'2-GaussianEst']);
addpath([fDir,'3-PiecewiseEst']);
addpath([fDir,'4-MixtureGaussian']);
addpath([fDir,'6-autorunInstructions']);

display(['Start time is ',datestr(now)]);

rng('shuffle');

%% Generate a bunch of datasets with varying lapse rates and prior widths

instrFxn    = str2func(instrName);
instrStr    = instrFxn();

% Extract settings from instruction script
domain      = instrStr.domain;
trPoolSize  = instrStr.trPoolSize;
priorF      = instrStr.priorF;
priorStruc  = instrStr.priorStruc;
lapseRate   = instrStr.lapseRate;
s1Vels      = instrStr.s1Vels;
s1Conts     = instrStr.s1Conts;
s2VelsDelta = instrStr.s2VelsDelta;
s2Conts     = instrStr.s2Conts;
model       = instrStr.model;
numTrials   = instrStr.numTrials;
numRuns     = instrStr.numRuns;
numReps     = instrStr.numReps;
tbtOn       = instrStr.tbtOn;
fitmethod   = instrStr.fitmethod;
interpOn    = instrStr.interpOn;
binScale    = instrStr.binScale;
contScF     = instrStr.contScF;
numComps    = instrStr.numComps;
estDistEval = instrStr.estDistEval;

% Update save dir with fit method
sDir = [sDir,fitmethod,filesep];

% Package up parameters
parStr.s1Vels       = s1Vels;
parStr.s1Conts      = s1Conts;
parStr.s2VelsDelta  = s2VelsDelta;
parStr.s2Conts      = s2Conts;
parStr.trPoolSize   = trPoolSize;
parStr.numTrials    = numTrials;
parStr.priorF       = priorF;
parStr.contScF      = contScF;
parStr.domain       = domain;
parStr.priorStruc   = priorStruc;
parStr.lapseRate    = lapseRate;

% Set up loop
indCnt = 1;

numSets = numel(lapseRate)*priorStruc.numPriors;
gtDatasets = cell(numSets,1);

% Loop over lapse rates
for i = 1:numel(lapseRate)
    
    % Loop over prior widths
    for j = 1:priorStruc.numPriors
        
        % Just give it a string id that's kinda unique
        priorVals = ['pInd_',num2str(j),'_',num2str(priorStruc.numPriors)];
        priorID   = [priorVals];
        fName     = ['gt_',domain,'_',priorF,'_',priorID,...
                     '_LR_',num2str(lapseRate(i))];

        gtDatasets{indCnt} = [sDir,fName];
        
        if exist([gtDatasets{indCnt},'/simData.mat'],'file') ~= 0
            
            % Check to see if existing dataset was run using the same set
            % of instructions (pass indices to only check the instructions
            % dealing with the data, not the estimation)
            temp = load([gtDatasets{indCnt},'/instrStr.mat']);
%             switch priorF
%                 case 'splitExp'
%                     test = structureComp(instrStr,temp.instrStr,[1:10]);
%                 otherwise           
%                     test = structureComp(instrStr,temp.instrStr,[1:10]);
%             end
            test = 1;
            if ~test
%                 error('Existing dataset run using different parameters, move or delete and try again.');
            else
                disp('Simulated dataset already exists, skipping.');
            end
        else

            disp(['Generating data set ',num2str(indCnt),'/',num2str(numSets)]);
            
            %%%%% should fix this bit so it's less awkward - need to go
            %%%%% into autogengroundtruth and gengroundtruth
            parStr.priorInd = j;
            parStr.LRind = i;
            
            autoGenGroundTruth(sDir,fName,parStr);
            
            save([sDir,filesep,fName,filesep,'instrStr'],'instrStr');
            
        end
        
        indCnt = indCnt + 1;
        
    end
    
end


%% Make estimates of parameters using different number of trials

%%%%%%%%%%%%%%%%%%%%
% Start parallel pool, define number of cores to use (24 on savio2)

numCores = feature('numcores');

if numCores > 4
    nCores = numCores;
    maxNumCompThreads(nCores);
    parpool(nCores);
    
    useParal = false;
else
    % Use parallel in fmincon
    useParal = true;
end

% Make set of individual worker streams instead of nested loop
streams     = makeCombos([numel(gtDatasets) numel(numTrials) numRuns]);
numStreams  = size(streams,1);
dsStream    = streams(:,1);
trStream    = streams(:,2);
rsxStream   = streams(:,3);
%%%%%%%%%%%%%%%%%%%%

% Define model estimation dirs
estDir = cell(lapseRate*priorStruc.numPriors);
indCnt = 1;

% Loop over lapse rates
for i = 1:numel(lapseRate)
    
    % Loop over prior widths
    for j = 1:priorStruc.numPriors
        
        % Just give it a string id that's kinda unique
        priorVals = ['pInd_',num2str(j),'_',num2str(priorStruc.numPriors)];
        priorID   = [priorVals];
        fName     = ['gt_',domain,'_',priorF,'_',priorID,'_'...
                     'LR_',num2str(lapseRate(i)),'_',model,'_',datestr(now,29)];
        
        estDir{indCnt} = [sDir,fName];
        
        if ~exist(estDir{indCnt})
            mkdir(estDir{indCnt});
        end 
        
        indCnt = indCnt + 1;
        
    end
    
end

% RESAMPLE AND TRIM DOWN STRUCTURES PASSED TO WORKERS
%%% Setting things up like this helps solve parallel worker
%%% failures, but be careful if re-running after cluster timeout - then you
%%% do need to reshuffle. keeping simData2 also allows for proper
%%% comparison between methods. 

% Loop over worker streams
for ws = 1:numStreams
    
    ds = dsStream(ws);
    nc = trStream(ws);
    
    simDataPath = gtDatasets{ds};
    numTrialsRS = numTrials(nc);
    
    if exist([simDataPath,filesep,'simDataN',num2str(numTrialsRS),'.mat'],'file') == 0
        
        tempStruc   = load([simDataPath,filesep,'simData']);
        
        % Resample out of parfor loop
        vStim1          = tempStruc.vStim1;
        cStim1          = tempStruc.cStim1;
        vStim2Delta     = tempStruc.vStim2Delta;
        cStim2          = tempStruc.cStim2;
        trPoolSize      = tempStruc.trPoolSize;
        simData         = tempStruc.simData;
        
%         % ID where likelihoods jump out of expected bin
%         %%% need to fix this script to work for all priors
%         if strcmp(priorF,'splitExp')
%             [simData] = IDbinJump(simData,'minDist',splitLoc);
%         end
        
        % Loop over independent samples
        for rsx = 1:numRuns
            [rsxTemp,rsxTemp2,lm1,lm2]  = choiceResampler(cStim1,cStim2,vStim1,vStim2Delta,...
                trPoolSize,simData,numTrialsRS);
            ptvs_data{rsx}      = rsxTemp;
            ptvs_datatbt{rsx}   = rsxTemp2;
            likeMeans1{rsx}     = lm1;
            likeMeans2{rsx}     = lm2;
        end
        
        simData2.vStim1         = vStim1;
        simData2.cStim1         = cStim1;
        simData2.vStim2Delta    = vStim2Delta;
        simData2.cStim2         = cStim2;
        simData2.numtrials      = numTrialsRS;
        simData2.ptvs_data      = ptvs_data;
        simData2.ptvs_datatbt   = ptvs_datatbt;
        simData2.parVec         = tempStruc.parVec;
        simData2.priorF         = simData(1).priorF;
        simData2.likeMeans1     = likeMeans1;
        simData2.likeMeans2     = likeMeans2;
        
        save([simDataPath,filesep,'simDataN',num2str(numTrialsRS)],'simData2');
        
    end
    
    if exist([estDir{ds},filesep,'initPars','.mat'],'file') == 0  
        
        load([simDataPath,filesep,'simDataN',num2str(numTrialsRS)],'simData2');
        
        % Pre-generate initial parameter vectors to focus variability on
        % samples, not initialization positions in parameter space
        numS1Vels  = numel(simData2.vStim1);
        numS2Conts = numel(simData2.cStim2);
        
        initPars.avloginit = rand(numReps,numS1Vels)*(5-(-5)) +  -5;
        % Initialize so first passes linear inequality
        offsetSlope        = repmat(0.1-0.01*(1:numS2Conts),[numReps 1]);
        initPars.hcinit    = rand(numReps,numS2Conts)*(2-0.1) + 0.1 + offsetSlope;
        initPars.sigPinit  = rand(numReps,numComps)*(1-0.1)   + 0.1;
        initPars.muPinit   = rand(numReps,numComps)*5         + 0.0;
        initPars.wPinit    = rand(numReps,numComps);
        
        % Pass number of components TO FIT
        initPars2.numComps  = numComps;
        
        save([estDir{ds},filesep,'initPars'],'initPars');
        
        clear sdStruct
        
    end
    
end

% Loop over parallel streams

estOpts = {tbtOn,fitmethod,interpOn,binScale,useParal,numReps,numComps,estDistEval};


if numCores > 4
    % Parallelize in loop for cluster
    parfor workStr = 1:numStreams
        
        dsi = dsStream(workStr);
        trj = trStream(workStr);
        rsx = rsxStream(workStr);
        
        simDataPath = gtDatasets{dsi};
        
        disp(['Fitting data set ',num2str(workStr),'/',num2str(numStreams)]);
        
        if tbtOn
            if interpOn == 1
                outFName = [model,'tbtInterp_jitter_n',num2str(numTrials(trj))];
            elseif interpOn == 2
                outFName = [model,'tbtNN_jitter_n',num2str(numTrials(trj))];
            else
                outFName = [model,'tbt_jitter_n',num2str(numTrials(trj))];
            end
        else
            if binScale
                outFName = [model,'bsc_jitter_n',num2str(numTrials(trj))];
            else
                outFName = [model,'_jitter_n',num2str(numTrials(trj))];
            end
        end
        
        tstart = tic;

        switch model
            case 'ss'
                autoEstParametersSS2006Vfix(simDataPath,outFName,rsx,numTrials(trj),estOpts,estDir{dsi});
            case 'mg'
                autoEstParametersMixGaussVfix(simDataPath,outFName,rsx,numTrials(trj),estOpts,estDir{dsi});
            case 'mgfix'
                autoEstParametersMixGaussBFit(simDataPath,outFName,rsx,numTrials(trj),estOpts,estDir{dsi});
            case 'mgfixneg'
                autoEstParametersMixGaussBFitNeg(simDataPath,outFName,rsx,numTrials(trj),estOpts,estDir{dsi});
            case 'zeroMu'
                autoEstParametersMixGaussZeroMu(simDataPath,outFName,rsx,numTrials(trj),estOpts,estDir{dsi});
            case 'zeroMuWeq'
                autoEstParametersMixGaussZeroMuWEq(simDataPath,outFName,rsx,numTrials(trj),estOpts,estDir{dsi});
        end
        
        time = toc(tstart);
        
        disp(['(Set: ',num2str(workStr),'/',num2str(numStreams),...
            ') Total time Elapsed: ',num2str(time)]);
        
    end
else
    % Parallelize in fmincon for PC
    for workStr = 1:numStreams

        dsi = dsStream(workStr);
        trj = trStream(workStr);
        rsx = rsxStream(workStr);
        
        simDataPath = gtDatasets{dsi};
        
        disp(['Fitting data set ',num2str(workStr),'/',num2str(numStreams)]);
        
        if tbtOn
            if interpOn == 1
                outFName = [model,'tbtInterp_jitter_n',num2str(numTrials(trj))];
            elseif interpOn == 2
                outFName = [model,'tbtNN_jitter_n',num2str(numTrials(trj))];
            else
                outFName = [model,'tbt_jitter_n',num2str(numTrials(trj))];
            end
        else
            if binScale
                outFName = [model,'bsc_jitter_n',num2str(numTrials(trj))];
            else
                outFName = [model,'_jitter_n',num2str(numTrials(trj))];
            end
        end
        
        tstart = tic;
        
        switch model
            case 'ss'
                autoEstParametersSS2006Vfix(simDataPath,outFName,rsx,numTrials(trj),estOpts,estDir{dsi});
            case 'mg'
                autoEstParametersMixGaussVfix(simDataPath,outFName,rsx,numTrials(trj),estOpts,estDir{dsi});
            case 'mgfix'
                autoEstParametersMixGaussBFit(simDataPath,outFName,rsx,numTrials(trj),estOpts,estDir{dsi});
            case 'mgfixneg'
                autoEstParametersMixGaussBFitNeg(simDataPath,outFName,rsx,numTrials(trj),estOpts,estDir{workStr});
            case 'zeroMu'
                autoEstParametersMixGaussZeroMu(simDataPath,outFName,rsx,numTrials(trj),estOpts,estDir{dsi});
            case 'zeroMuWeq'
                autoEstParametersMixGaussZeroMuWEq(simDataPath,outFName,rsx,numTrials(trj),estOpts,estDir{dsi});
        end
        
        time = toc(tstart);
        
        disp(['(Set: ',num2str(workStr),'/',num2str(numStreams),...
            ') Total time Elapsed: ',num2str(time)]);
        
    end
end

disp(['End time is ',datestr(now)]);

end