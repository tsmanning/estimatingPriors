function [fitPars] = autoEstParametersMixGaussBFitNeg(simDataPath,outFName,rsx,numTrials,estOpts,estDir)

% Run parameter estimation using mixture of Gaussian approch to prior
% estimation fitting ONLY the weights of an N component mixture
%
% Usage: [fitPars] = autoEstParametersMixGauss(simDataPath,outFName,rsx,numTrials,estOpts)

% Setup options
tbtOn       = estOpts{1};
fitmethod   = estOpts{2};
interpOn    = estOpts{3};
binScale    = estOpts{4};
useParal    = estOpts{5};
numReps     = estOpts{6};
numComp     = estOpts{7};

% Extract pars from structure
load([simDataPath,filesep,'simDataN',num2str(numTrials)]);
load([estDir,filesep,'initPars']);

parVecGT    = simData2.parVec;
priorF      = simData2.priorF;
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

numPars = numS1Vels + numS2Conts + 3*numComp;

numCompGT   = (numel(parVecGT) - numS2Conts)/3;

xfit = false;

% Flag that this is fitting a different gt function than the method and
% define number of parameters in gt function
if ~strcmp(priorF,'mixgauss') || (numPars ~= numel(parVecGT))
    
    xfit = true;
    
%     switch priorF
%         case 'gauss'
%             numPars = 14;
%         case 'piece'
%             numPars = 19;
%         case 'mixgauss'
% %             numPars = numel(parVecGT);
% 
%     end
end

numUniqCont = numel(unique([cStim1 cStim2]));

% Save ID for display to terminal during parallel fitting
tempID = regexp(simDataPath,filesep,'split');
dispID = tempID{end};


%% Make anonymous function for negative log-likelihood & fit model to data

nllOpts = {'mixgauss',fitmethod,interpOn,nllMethod,binScale};

vfix    = parVecGT(1:numS1Vels);
% Find close form solution for why 0.5 is the minimum?
sigfix  = 0.5*ones(1,numComp);
% Make the first component negative to avoid bias prior peak, set last
% component at the highest test velocity
mubnds  = getLogXform([0 vStim1*vStim2Delta(end)],0.3);
mufix   = linspace(mubnds(1),mubnds(2),numComp-1);
mufix   = [mufix(1)-diff(mufix(1:2)) mufix];

% Allow only 1 sigma and the weights to be fit freely
%%%% parVec = [hc sigma(1) weights];
lossFun = @(parVec)(calculate_model_nll([vfix parVec(1:numS2Conts) sigfix mufix parVec(numS2Conts+1:end)],ptvs_data,...
                    vStim1,cStim1,vStim2Delta,cStim2,numTrials,nllOpts));
     
% parVec = [1x6 gvlog, 1x7 hc, 1xNumComp];
                
% Fminunc/fmincon
opts    = optimset('display','off','useparallel',useParal,'tolx',1e-13,...
                   'maxfunevals',1e4,'largescale', 'off');
               
%%%% adding constraint that best fit sigma must be at least as big as the smallest interval between test speeds
% smallInt = min(diff(getLogXform(vStim1(1)*vStim2Delta,0.3)))/vfix(1);
lb = [0.001*ones(1,numUniqCont) zeros(1,numComp)];

% lb = [0.001*ones(1,numUniqCont) 1e-5*ones(1,numComp) -eps*ones(1,numComp) zeros(1,numComp)];
ub = [100*ones(1,numUniqCont) ones(1,numComp)];               
               
Aeq = [zeros(1,numUniqCont) ones(1,numComp)];
beq = 1;

% Make sure likelihood width is a monotonic decreasing function of contrast
A = [-eye(numUniqCont-1) zeros(numUniqCont-1,1)] + [zeros(numUniqCont-1,1) eye(numUniqCont-1)];
A = [A zeros(numUniqCont-1,numComp)];
b = zeros(numUniqCont-1,1);


%% Fit data

indCnt  = 1;

gvlog       = nan(numReps,numS1Vels);
gvlogF      = nan(numReps,numS1Vels);
hc          = nan(numReps,numUniqCont);
hcF         = nan(numReps,numUniqCont);
sigP        = nan(numReps,numComp);
sigPF       = nan(numReps,numComp);
muP         = nan(numReps,numComp);
muPF        = nan(numReps,numComp);
wP          = nan(numReps,numComp);
wPF         = nan(numReps,numComp);
ptvsF       = cell(numReps,1);
nllF        = nan(numReps,1);
initConds   = nan(numReps,numPars);

while indCnt <= numReps
    
    disp(['Run ',num2str(indCnt),'/',num2str(numReps),...
          ' - ',dispID,', n=',num2str(numTrials)]);
      
    tsl = tic;
    
    % Grab preallocated vals to test binomial/bernoulli
    if isfield(initPars,'hcinit')
        hcinit          = initPars.hcinit(indCnt,:);
        sigPinit        = initPars.sigPinit(indCnt,:);
        muPinit         = initPars.muPinit(indCnt,:);
        wPinit          = initPars.wPinit(indCnt,:);
        
        parVec0         = [hcinit wPinit];
    else
        error('No pre allocated parameter vector');
    end

    % Run optimization
    try
        [prshat,~,exFlag,outStruc]  = fmincon(lossFun,parVec0,A,b,Aeq,beq,lb,ub,[],opts);
    catch
        prshat                   = nan(1,numel(parVec0));
        exFlag                   = nan;
        outStruc.iterations      = [];
        outStruc.funcCount       = [];
        outStruc.constrviolation = [];
        outStruc.stepsize        = [];
        outStruc.algorithm       = [];
        outStruc.firstorderopt   = [];
        outStruc.cgiterations    = [];
        outStruc.message         = 'failed to initialize pars';
    end
    
    % Extract best fit pars
    gvlogF(indCnt,:) = vfix;
    hcF(indCnt,:)    = prshat(1:numUniqCont);
%     sigPF(indCnt,:)  = prshat(numUniqCont + 1:numUniqCont + numComp);
%     muPF(indCnt,:)   = prshat(numUniqCont + numComp + 1:numUniqCont + 2*numComp);
%     wPF(indCnt,:)    = prshat(numUniqCont + 2*numComp + 1:numUniqCont + 3*numComp);
    sigPF(indCnt,:)  = sigfix;
    muPF(indCnt,:)   = mufix;
    wPF(indCnt,:)    = prshat(numUniqCont + 1:end);
    
    initConds(indCnt,:) = [vfix parVec0(1:numUniqCont) sigfix mufix parVec0(numUniqCont + 1:end)];
    os(indCnt)          = outStruc;
    
    % Output best fitting psychometric function and nLL
    ptvsF{indCnt}   = calculate_ptvs_mixgauss(vStim1,cStim1,vStim2Delta,cStim2,...
                                                gvlogF(indCnt,:),hcF(indCnt,:),...
                                                sigPF(indCnt,:),muPF(indCnt,:),...
                                                wPF(indCnt,:),fitmethod);
    nllF(indCnt)    = lossFun(prshat);
    
    % Unpack GT parameters
    if xfit
        switch priorF
%             case 'gauss'
%                 nllOptsXfit   = {'gauss',fitmethod,interpOn,nllMethod,binScale};
%                 lossFunXfit   = @(parVec)(calculate_model_nll([vfix parVec],ptvs_data,...
%                     vStim1,cStim1,vStim2Delta,cStim2,numTrials,nllOptsXfit));
%                 
%                 nllGT(indCnt) = lossFunXfit(parVecGT);
%                 
%                 gvlog         = parVecGT(1:numS1Vels);
%                 hc            = parVecGT(numS1Vels + 1:numS1Vels + numUniqCont);
%                 sigP          = parVecGT(end);
                
            case 'expon'
                nllOptsXfit   = {'expon',fitmethod,interpOn,nllMethod,binScale};
                lossFunXfit   = @(parVec)(calculate_model_nll([vfix parVec],ptvs_data,...
                    vStim1,cStim1,vStim2Delta,cStim2,numTrials,nllOptsXfit));
                
                nllGT(indCnt) = lossFunXfit(parVecGT);
                
                gvlog         = parVecGT(1:numS1Vels);
                hc            = parVecGT(numS1Vels + 1:2*numS1Vels);
                avlog         = parVecGT(2*numS1Vels + 1:end);
        
            case 'mixgauss'
                %%%% this assumes that the GT function is the mixgauss
                %%%% version with free pars = 3*components
                nllOptsXfit   = {'mixgauss',fitmethod,interpOn,nllMethod,binScale};
                lossFunXfit   = @(parVec)(calculate_model_nll(parVec,ptvs_data,...
                    vStim1,cStim1,vStim2Delta,cStim2,numTrials,nllOptsXfit));
                
                nllGT(indCnt) = lossFunXfit(parVecGT);
                
                numCompGT   = (numel(parVecGT) - numS2Conts - numS1Vels)/3;
                
                gvlog         = parVecGT(1:numS1Vels);
                hc            = parVecGT(numS1Vels + 1:numS1Vels + numUniqCont);
                sigP          = parVecGT(numS1Vels + numUniqCont + 1:numS1Vels + numUniqCont + numCompGT);
                muP           = parVecGT(numS1Vels + numUniqCont + numCompGT + 1:numS1Vels + numUniqCont + 2*numCompGT);
                wP            = parVecGT(numS1Vels + numUniqCont + 2*numCompGT + 1:numS1Vels + numUniqCont + 3*numCompGT);
                
        end
    else
        gvlog         = parVecGT(1:numS1Vels);
        hc            = parVecGT(numS1Vels + 1:numS1Vels + numUniqCont);
        sigP          = parVecGT(numS1Vels + numUniqCont + 1:numS1Vels + numUniqCont + numComp);
        muP           = parVecGT(numS1Vels + numUniqCont + numComp + 1:numS1Vels + numUniqCont + 2*numComp);
        wP            = parVecGT(numS1Vels + numUniqCont + 2*numComp + 1:numS1Vels + numUniqCont + 3*numComp);
        
        freePars      = [parVecGT(numS1Vels + 1) parVecGT(end-(numCompGT-1):end)];
        
        nllGT(indCnt) = lossFun(freePars);
    end 
    
    %%%% turn this off for now, prob not needed.
%     if exFlag > 0
        indCnt = indCnt + 1;
%     elseif exFlag == 0
%         disp('Max function evals hit, repeat loop.');
%     end

    t = toc(tsl);
    disp(['Time Elapsed: ',num2str(t)]);
    
end


%% Save parameter estimates for this run
parEstPath = [estDir,filesep,'parEsts_',outFName,'_',num2str(rsx),'.mat'];

% % Check to see if a run with this name already exists
% if numel(dir(parEstPath)) ~= 0
%     
%     % If so, grab the appended number (if any) and increment it
%     % RegEx: grab numbers preceded by 'path/outFName_' and followed by '.mat'
%     lDir = ls([estDir,filesep,'parEsts_*']);
%     
%     appNumStr = regexp(lDir,['(?<=.+',outFName,'_)[\w]*(?=.mat)'],'match');
%    
%     if ~isempty(appNumStr)
%         
%         temp = cellfun(@(x) str2num(x),appNumStr);
%         appNum = max(temp);
%         appNum = appNum + 1;
%         appNumStr = num2str(appNum);
%         
%     else
%         
%         appNumStr = '2';
%         
%     end
%     
%     parEstPath = [estDir,filesep,'parEsts_',outFName,'_',appNumStr,'.mat'];
%     
% end

if xfit
    switch priorF
%         case 'gauss'
%             save(parEstPath,'numTrials','gvlog','hc','sigP','avlog','initConds','gvlogF',...
%                             'hcF','avlogF','ptvsF','ptvs_data','nllF','nllGT','os',...
%                             'parVecGT','estOpts','rsx');
%             
%             fitPars.gvlog       = gvlog;
%             fitPars.hc          = hc;
%             fitPars.sigP        = sigP;
            
        case 'expon'
            save(parEstPath,'numTrials','gvlog','hc','avlog','initConds','gvlogF',...
                            'hcF','avlogF','ptvsF','ptvs_data','nllF','nllGT','os',...
                            'parVecGT','estOpts','rsx');
            
            fitPars.gvlog       = gvlog;
            fitPars.hc          = hc;
            fitPars.avlog       = avlog;
            
        case 'mixgauss'
            save(parEstPath,'numTrials','gvlog','hc','sigP','muP','wP','initConds','gvlogF',...
                            'hcF','sigPF','muPF','wPF','ptvsF','ptvs_data','nllF','nllGT',...
                            'os','parVecGT','estOpts','rsx');

            fitPars.gvlog       = gvlog;
            fitPars.hc          = hc;
            fitPars.sigP        = sigP;
            fitPars.muP         = muP;
            fitPars.wP          = wP;
                
    end
else
    save(parEstPath,'numTrials','gvlog','hc','sigP','muP','wP','initConds','gvlogF',...
                    'hcF','sigPF','muPF','wPF','ptvsF','ptvs_data','nllF','nllGT',...
                    'os','parVecGT','estOpts','rsx');

    fitPars.gvlog       = gvlog;
    fitPars.hc          = hc;
    fitPars.sigP        = sigP;
    fitPars.muP         = muP;
    fitPars.wP          = wP;
end

fitPars.parVecGT    = parVecGT;
fitPars.numTrials   = numTrials;
fitPars.initConds   = initConds;
fitPars.gvlogF      = gvlogF;
fitPars.hcF         = hcF;
fitPars.sigPF       = sigPF;
fitPars.muPF        = muPF;
fitPars.wPF         = wPF;
fitPars.ptvsF       = ptvsF;
fitPars.ptvs_data   = ptvs_data;
fitPars.nllF        = nllF;
fitPars.nllGT       = nllGT;
fitPars.fminuncOut  = os;

end