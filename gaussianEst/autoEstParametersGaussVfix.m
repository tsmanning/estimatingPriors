function [fitPars] = autoEstParametersGaussVfix(simDataPath,outFName,rsx,numTrials,estOpts)

% Run model fitting for Gaussian parameterization of prior
%
% Usage: [fitPars] = autoEstParametersGauss(simDataPath,outFName,numTrials)

% Setup options
tbtOn       = estOpts{1};
fitmethod   = estOpts{2};
interpOn    = estOpts{3};
binScale    = estOpts{4};
useParal    = estOpts{5};
numReps     = estOpts{6};

% Extract pars from structure
load([simDataPath,filesep,'simDataN',num2str(numTrials)]);

parVecGT    = simData2.parVec;
vStim1      = simData2.vStim1;
cStim1      = simData2.cStim1;
vStim2Delta = simData2.vStim2Delta;
cStim2      = simData2.cStim2;
ptvsFull    = simData2.ptvs_data;
ptvsFulltbt = simData2.ptvs_datatbt;

% Grab resampled set for this run
if tbtOn
    ptvs_data   = ptvsFulltbt{rsx};
    nllMethod   = 'bernoulli';
else
    ptvs_data   = ptvsFull{rsx};
    nllMethod   = 'binomial';
end

numS2Conts = numel(cStim2);
numS1Vels = numel(vStim1);

numPars = numel(parVecGT);

xfit = false;

if ~strcmp(priorF,'gauss')
    
    xfit = true;
    
    switch priorF
        case 'mixgauss'
            numPars = 19;
        case 'piece'
            numPars = 19;
    
    end
end

numUniqCont = numel(unique([cStim1 cStim2]));

% Save ID for display to terminal during parallel fitting
tempID = regexp(simDataPath,filesep,'split');
dispID = tempID{end};


%% Make anonymous function for negative log-likelihood & fit model to data

% Fix velocity-dependent component of likelihood widths
% if ~xfit
%     % Gaussian: parVec = [gvlog hc sigP]
    vfix = parVecGT(1:numS1Vels);
% else
%     % Piecewise: parVec = [avlog gvlog hc];
%     vfix = parVecGT(numS1Vels+1:2*numS1Vels);
% end

nllOpts = {'gauss',fitmethod,interpOn,nllMethod,binScale};

lossFun = @(parVec)(calculate_model_nll([vfix parVec],ptvs_data,...
                    vStim1,cStim1,vStim2Delta,cStim2,numTrials,nllOpts));

% Fminunc/fmincon
opts    = optimset('display','off','useparallel',useParal,'tolx',1e-13,...
                   'maxfunevals',1e4,'largescale', 'off');
% lb = [zeros(1,numUniqCont) 0];
% ub = inf*ones(1,numUniqCont + 1);
lb = [0.001*ones(1,numUniqCont) 1e-5];
ub = [100*ones(1,numUniqCont) 1e5];

% fminsearch
% opts    = optimset('display','iter','tolx',1e-4,'maxfunevals',15e4);

%% Fit data

indCnt  = 1;

gvlog       = nan(numReps,numS1Vels);
gvlogF      = nan(numReps,numS1Vels);
hc          = nan(numReps,numUniqCont);
hcF         = nan(numReps,numUniqCont);
sigP        = nan(numReps,1);
sigPF       = nan(numReps,1);
ptvsF       = cell(numReps,1);
nllF        = nan(numReps,1);
initConds   = nan(numReps,numPars);

while indCnt <= numReps
    
    disp(['Run ',num2str(indCnt),'/',num2str(numReps),...
          ' - ',dispID,', n=',num2str(numTrials)]);
      
    tsl = tic;
    
    %%%%%%%%%%%%%% Choose Initialization %%%%%%%%%%%%%%%
    
    % Randomize initial parameters within reasonable window
    %%% grab preallocated vals to test binomial/bernoulli
    if isfield(simData2,'hcinit')
        hcinit          = simData2.hcinit(indCnt,:);
        sigPinit        = simData2.sigPinit(indCnt,:);
    else
        error('No pre allocated parameter vector');
%         hcinit          = rand(1,numS2Conts)*(2-0.1) + 0.1;
%         sigPinit        = rand(1,1)*(1-0.1)          + 0.1;
    end
%     hcinit          = ones(1,numS2Conts);
%     sigPinit        = 10*ones(1,1);
    
    parVec0         = [hcinit sigPinit];
   
%     % Use GT for testing
%     if xfit
%         % No real ground truth for fitting Gaussian to non-Gaussian
%         parVecInit  = [parVecGT(numS1Vels + 1:end) 1];
%     else
%         parVecInit  = parVecGT;
%     end
%     parVec0                     = parVecInit.*(rand(1,numPars) + 0.5);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     [prshat,~,exFlag,outStruc]  = fminunc(lossFun,parVec0,opts);
%     [prshat,~,exFlag,outStruc]  = fminsearch(lossFun,parVec0,opts);

try
    [prshat,~,exFlag,outStruc]  = fmincon(lossFun,parVec0,[],[],[],[],lb,ub,[],opts);
catch
    %%% put below back in once done comparing bernoulli and binomial
%     try
%         hcinit          = rand(1,numS2Conts)*(2-0.1) + 0.1;
%         sigPinit        = rand(1,1)*(1-0.1)          + 0.1;
%         parVec0         = [hcinit sigPinit];
%         [prshat,~,exFlag,outStruc]  = fmincon(lossFun,parVec0,[],[],[],[],lb,ub,[],opts);
%     catch
%         % if rerolling the initial conditions fails twice in a row,
%         % just give up.
%         
        prshat   = nan(1,numel(parVec0));
        exFlag   = nan;
        outStruc.iterations = [];
        outStruc.funcCount = [];
        outStruc.constrviolation = [];
        outStruc.stepsize = [];
        outStruc.algorithm = [];
        outStruc.firstorderopt = [];
        outStruc.cgiterations = [];
        outStruc.message = 'failed to initialize pars';
%     end
end
    
    gvlogF(indCnt,:) = vfix;
    hcF(indCnt,:)    = prshat(1:numUniqCont);
    sigPF(indCnt,1)  = prshat(end);
    
    initConds(indCnt,:) = [vfix parVec0];
    os(indCnt)          = outStruc;    
    
    % Output best fitting psychometric function and nLL
    ptvsF{indCnt}   = calculate_ptvs_gauss(vStim1,cStim1,vStim2Delta,cStim2,...
                                         gvlogF(indCnt,:),hcF(indCnt,:),sigPF(indCnt,1),fitmethod);
    nllF(indCnt)    = lossFun(prshat);
                        
    % Unpack GT parameters
    if xfit
        nllOptsXfit   = {'piece',fitmethod,interpOn,nllMethod,binScale};
        lossFunXfit   = @(parVec)(calculate_model_nll([vfix parVec],ptvs_data,...
                        vStim1,cStim1,vStim2Delta,cStim2,numTrials,nllOptsXfit));
                        
        nllGT(indCnt) = lossFunXfit(parVecGT);
        
        gvlog         = parVecGT(1:numS1Vels);
        hc            = parVecGT(numS1Vels + 1:numS1Vels + numUniqCont);
        avlog         = parVecGT(numS1Vels + numUniqCont + 1:end);
    else
        gvlog         = parVecGT(1:numS1Vels);
        hc            = parVecGT(numS1Vels + 1:numS1Vels + numUniqCont);
        sigP          = parVecGT(end);
        
        freePars      = [hc sigP];
        
        nllGT(indCnt) = lossFun(freePars);
    end
        
    % Only increment if optimization didn't get stuck in a black hole
    % CAREFUL HERE FOR DEBUGGING; fminsearch exits with 0 even for 50k+
    % iterations
    %%%% turned off for binomial/bernoulli comparison
%     if exFlag > 0
        indCnt = indCnt + 1;
%     elseif exFlag == 0
%         disp('Max function evals hit, repeat loop.');
%     end
    
    t = toc(tsl);
    disp(['Time Elapsed: ',num2str(t)]);
    
end

%% Save parameter estimates for this run
parEstPath = [simDataPath,filesep,'parEsts_',outFName,'.mat'];

% Check to see if a run with this name already exists
if numel(dir(parEstPath)) ~= 0
    
    % If so, grab the appended number (if any) and increment it
    % RegEx: grab numbers preceded by 'path/outFName_' and followed by '.mat'
    lDir = ls([simDataPath,filesep,'parEsts_*']);
    
    appNumStr = regexp(lDir,['(?<=.+',outFName,'_)[\w]*(?=.mat)'],'match');
   
    if ~isempty(appNumStr)
        
        temp = cellfun(@(x) str2num(x),appNumStr);
        appNum = max(temp);
        appNum = appNum + 1;
        appNumStr = num2str(appNum);
        
    else
        
        appNumStr = '2';
        
    end
    
    parEstPath = [simDataPath,filesep,'parEsts_',outFName,'_',appNumStr,'.mat'];
    
end

if xfit
    save(parEstPath,'numTrials','gvlog','hc','avlog','initConds','gvlogF',...
        'hcF','sigPF','ptvsF','ptvs_data','nllF','nllGT','os','parVecGT');
    
    fitPars.gvlog       = gvlog;
    fitPars.hc          = hc;
    fitPars.avlog       = avlog;
else
    save(parEstPath,'numTrials','gvlog','hc','sigP','initConds','gvlogF',...
        'hcF','sigPF','ptvsF','ptvs_data','nllF','nllGT','os','parVecGT');

    fitPars.gvlog       = gvlog;
    fitPars.hc          = hc;
    fitPars.sigP        = sigP;
end

fitPars.parVecGT    = parVecGT;
fitPars.numTrials   = numTrials;
fitPars.initConds   = initConds;
fitPars.gvlogF      = gvlogF;
fitPars.hcF         = hcF;
fitPars.sigPF       = sigPF;
fitPars.ptvsF       = ptvsF;
fitPars.ptvs_data   = ptvs_data;
fitPars.nllF        = nllF;
fitPars.nllGT       = nllGT;
fitPars.fminuncOut  = os;

end
