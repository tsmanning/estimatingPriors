function [fitPars] = autoEstParametersSS2006Vfix(simDataPath,outFName,rsx,numTrials)

% Run model fitting for piecewise estimation of prior
%
% Usage: [fitPars] = autoEstParametersSS2006(simDataPath,outFName,numTrials)

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

if numPars == 8
    numPars = 13;
    xfit = true;
end

numUniqCont = numel(unique([cStim1 cStim2]));

% Save ID for display to terminal during parallel fitting
tempID = regexp(simDataPath,filesep,'split');
dispID = tempID{end};


%% Make anonymous function for negative log-likelihood & fit model to data
vfix = parVecGT(1:numS1Vels);
lossFun = @(parVec)(calculate_model_nll(...
                    [parVec(1:numS1Vels) vfix parVec(numS1Vels+1:numS1Vels+numUniqCont)],...
                    ptvs_data,vStim1,cStim1,vStim2Delta,cStim2,numTrials));

% Fminunc/fmincon
opts    = optimset('display','off','tolx',1e-13,...
                   'maxfunevals',1e4,'largescale', 'off');
lb = [-inf*ones(1,numS1Vels) zeros(1,numUniqCont)];
ub = inf*ones(1,numS1Vels + numUniqCont);

% fminsearch
% opts    = optimset('display','iter','tolx',1e-4,'maxfunevals',15e4);

%% Fit data

numReps = 5;
% numReps = 1;
indCnt  = 1;

avlogF      = nan(numReps,numS1Vels);
gvlogF      = nan(numReps,numS1Vels);
hcF         = nan(numReps,numUniqCont);
ptvsF       = cell(numReps,1);
nllF        = nan(numReps,1);
initConds   = nan(numReps,numPars);

while indCnt <= numReps
    
    disp(['Run ',num2str(indCnt),'/',num2str(numReps),...
          ' - ',dispID,', n=',num2str(numTrials)]);
      
    tsl = tic;
    
    %%%%%%%%%%%%%% Choose Initialization %%%%%%%%%%%%%%%
    
   % Randomize initial parameters within reasonable window
   avloginit       = rand(1,numS1Vels)*(5-(-5)) +  -5;
   gvloginit       = rand(1,numS1Vels)*(2-0.1)  + 0.1;
   hcinit          = rand(1,numS2Conts)*(2-0.1) + 0.1;
   
   parVec0         = [avloginit gvloginit hcinit];
    
%      % Use GT for testing
%      if xfit
%          parVecInit  = [-(1/parVecGT(end)^2)*log(1 + vStim1/0.3) parVecGT(1:numS1Vels+numS2Conts)];
%      else
%          parVecInit  = parVecGT;
%      end
%      parVec0                     = parVecInit.*(rand(1,numPars) + 0.5);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [prshat,~,exFlag,outStruc]  = fminunc(lossFun,parVec0,opts);
%     [prshat,~,exFlag,outStruc]  = fminsearch(lossFun,parVec0,opts);
    
    avlogF(indCnt,:) = prshat(1:numS1Vels);
    gvlogF(indCnt,:) = prshat(numS1Vels + 1:2*numS1Vels);
    hcF(indCnt,:)    = prshat(2*numS1Vels + 1:numPars);

    initConds(indCnt,:) = parVec0;
    os(indCnt)          = outStruc;
    
    % Output best fitting psychometric function and nLL
    ptvsF{indCnt}     = calculate_ptvs(vStim1,cStim1,vStim2Delta,cStim2,...
                                       avlogF(indCnt,:),gvlogF(indCnt,:),hcF(indCnt,:));
    nllF(indCnt)      = lossFun(prshat);
    
    % Unpack GT parameters
    if xfit
        lossFunXfit   = @(parVec)(calculate_model_nll_gauss(parVec,ptvs_data,...
                                  vStim1,cStim1,vStim2Delta,cStim2,numTrials));
        nllGT(indCnt) = lossFunXfit(parVecGT);
        
        gvlog         = parVecGT(1:numS1Vels);
        hc            = parVecGT(numS1Vels + 1:numS1Vels + numUniqCont);
        sigP          = parVecGT(end);
        avlog         = -(1/sigP^2)*log(1 + vStim1/0.3);
    else
        avlog         = parVecGT(1:numS1Vels);
        gvlog         = parVecGT(numS1Vels + 1:2*numS1Vels);
        hc            = parVecGT(2*numS1Vels + 1:numPars);
        nllGT(indCnt) = lossFun(parVecGT);
    end
    
    % Only increment if optimization didn't get stuck in a black hole
    % Fminsearch almost always hits ceiling so don't even bother.
    if (exFlag > 0) || ~(strcmp(os(indCnt),'quasi-newton'))
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
    save(parEstPath,'numTrials','gvlog','hc','sigP','avlog','initConds','gvlogF',...
        'hcF','avlogF','ptvsF','ptvs_data','nllF','nllGT','os','parVecGT');
    
    fitPars.gvlog       = gvlog;
    fitPars.hc          = hc;
    fitPars.sigP        = sigP;
else
    save(parEstPath,'numTrials','gvlog','hc','avlog','initConds','gvlogF',...
        'hcF','avlogF','ptvsF','ptvs_data','nllF','nllGT','os','parVecGT');

    fitPars.gvlog       = gvlog;
    fitPars.hc          = hc;
    fitPars.avlog       = avlog;
end

fitPars.parVecGT    = parVecGT;
fitPars.numTrials   = numTrials;
fitPars.initConds   = initConds;
fitPars.gvlogF      = gvlogF;
fitPars.hcF         = hcF;
fitPars.sigPF       = avlogF;
fitPars.ptvsF       = ptvsF;
fitPars.ptvs_data   = ptvs_data;
fitPars.nllF        = nllF;
fitPars.nllGT       = nllGT;
fitPars.fminuncOut  = os;

end
