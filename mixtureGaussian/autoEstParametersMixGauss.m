function [fitPars] = autoEstParametersMixGauss(simDataPath,outFName,rsx,numTrials,estOpts)

% Run parameter estimation using mixture of Gaussian approch to prior
% estimation
%
% Usage: [fitPars] = autoEstParametersMixGauss(simDataPath,outFName,rsx,numTrials,estOpts)

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

numPars = numel(parVecGT);

xfit = false;

if ~strcmp(priorF,'mog')
    
    xfit = true;
    
    switch priorF
        case 'gauss'
            numPars = 14;
        case 'piece'
            numPars = 19;
    
    end
end

numUniqCont = numel(unique([cStim1 cStim2]));

% Save ID for display to terminal during parallel fitting
tempID = regexp(simDataPath,filesep,'split');
dispID = tempID{end};


%% Make anonymous function for negative log-likelihood & fit model to data

vfix = parVecGT();




% Fminunc/fmincon
opts    = optimset('display','off','useparallel',useParal,'tolx',1e-13,...
                   'maxfunevals',1e4,'largescale', 'off');


%% Fit data




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
    save(parEstPath,'numTrials','gvlog','hc','sigP','avlog','initConds','gvlogF',...
        'hcF','avlogF','ptvsF','ptvs_data','nllF','nllGT','os','parVecGT','estOpts');
    
    fitPars.gvlog       = gvlog;
    fitPars.hc          = hc;
    fitPars.sigP        = sigP;
else
    save(parEstPath,'numTrials','gvlog','hc','avlog','initConds','gvlogF',...
        'hcF','avlogF','ptvsF','ptvs_data','nllF','nllGT','os','parVecGT','estOpts');

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