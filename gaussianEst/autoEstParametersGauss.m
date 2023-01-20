function [fitPars] = autoEstParametersGauss(simDataPath,outFName,rsx,numTrials)

% Run model fitting for Gaussian parameterization of prior
%
% Usage: [fitPars] = autoEstParametersGauss(simDataPath,outFName,numTrials)

% Extract pars from structure
load([simDataPath,filesep,'simData2']);

parVecGT    = simData2.parVec;
vStim1      = simData2.vStim1;
cStim1      = simData2.cStim1;
vStim2Delta = simData2.vStim2Delta;
cStim2      = simData2.cStim2;
ptvsFull    = simData2.ptvs_data;

% Grab resampled set for this run
ptvs_data   = ptvsFull{rsx};

numS2Conts = numel(cStim2);
numS1Vels = numel(vStim1);

numPars = numel(parVecGT);

xfit = false;

if numPars == 19
    numPars = 14;
    xfit = true;
end

numUniqCont = numel(unique([cStim1 cStim2]));

% Save ID for display to terminal during parallel fitting
tempID = regexp(simDataPath,filesep,'split');
dispID = tempID{end};


%% Make anonymous function for negative log-likelihood & fit model to data
lossFun = @(parVec)(calculate_model_nll_gauss(parVec,ptvs_data,vStim1,cStim1,vStim2Delta,cStim2,numTrials));

% Fminunc/fmincon
opts    = optimset('display','off','tolx',1e-13,...
                   'maxfunevals',1e4,'largescale', 'off');
lb = [zeros(1,numS1Vels) zeros(1,numUniqCont) 0];
ub = inf*ones(1,numS1Vels + numUniqCont + 1);

% fminsearch
% opts    = optimset('display','iter','tolx',1e-4,'maxfunevals',15e4);

%% Fit data

numReps = 12;
indCnt  = 1;

gvlog       = nan(numReps,numS1Vels);
hc          = nan(numReps,numUniqCont);
sigP        = nan(numReps,1);
ptvsF       = cell(numReps,1);
nllF        = nan(numReps,1);
initConds   = nan(numReps,numPars);

while indCnt <= numReps
    
    disp(['Run ',num2str(indCnt),'/',num2str(numReps),...
          ' - ',dispID,', n=',num2str(numTrials)]);
      
    tsl = tic;
    
    %%%%%%%%%%%%%% Choose Initialization %%%%%%%%%%%%%%%
    
    % Randomize initial parameters within reasonable window
    gvloginit       = rand(1,numS1Vels)*(2-0.1)  + 0.1;
    hcinit          = rand(1,numS2Conts)*(2-0.1) + 0.1;
    sigPinit        = rand(1,1)*(1-0.1)         + 0.1;
    
    parVec0         = [gvloginit hcinit sigPinit];
   
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
    [prshat,~,exFlag,outStruc]  = fmincon(lossFun,parVec0,[],[],[],[],lb,ub,[],opts);
    
    gvlogF(indCnt,:) = prshat(1:numS1Vels);
    hcF(indCnt,:)    = prshat(numS1Vels + 1:numS1Vels + numUniqCont);
    sigPF(indCnt,1)  = prshat(end);
    
    initConds(indCnt,:) = parVec0;
    os(indCnt)          = outStruc;    
    
    % Output best fitting psychometric function and nLL
    ptvsF{indCnt}   = calculate_ptvs_gauss(vStim1,cStim1,vStim2Delta,cStim2,...
                                         gvlogF(indCnt,:),hcF(indCnt,:),sigPF(indCnt,1));
    nllF(indCnt)    = lossFun(prshat);
                        
    % Unpack GT parameters
    if xfit
        lossFunXfit   = @(parVec)(calculate_model_nll(parVec,ptvs_data,...
                                  vStim1,cStim1,vStim2Delta,cStim2,numTrials));
        nllGT(indCnt) = lossFunXfit(parVecGT);
        
        avlog         = parVecGT(1:numS1Vels);
        gvlog         = parVecGT(numS1Vels + 1:2*numS1Vels);
        hc            = parVecGT(2*numS1Vels + 1:numPars);
    else
        gvlog         = parVecGT(1:numS1Vels);
        hc            = parVecGT(numS1Vels + 1:numS1Vels + numUniqCont);
        sigP          = parVecGT(end);
        nllGT(indCnt) = lossFun(parVecGT);
    end
        
    % Only increment if optimization didn't get stuck in a black hole
    % CAREFUL HERE FOR DEBUGGING; fminsearch exits with 0 even for 50k+
    % iterations
    if exFlag > 0
        indCnt = indCnt + 1;
    elseif exFlag == 0
        disp('Max function evals hit, repeat loop.');
    end
    
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
        
        appNum = str2num(appNumStr{end});
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