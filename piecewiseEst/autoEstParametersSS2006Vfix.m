function [fitPars] = autoEstParametersSS2006Vfix(simDataPath,outFName,rsx,numTrials,estOpts,estDir)

% Run model fitting for piecewise estimation of prior
%
% Usage: [fitPars] = autoEstParametersSS2006(simDataPath,outFName,numTrials)

% Setup options
tbtOn       = estOpts{1};
fitmethod   = estOpts{2};
interpOn    = estOpts{3};
binScale    = estOpts{4};
useParal    = estOpts{5};
numReps     = estOpts{6};

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

numParsGT = numel(parVecGT);
numPars = numS2Conts + 2*numS1Vels;

xfit = false;

if ~strcmp(priorF,'expon') && ~strcmp(priorF,'splitExp')
    
    xfit = true;
    
%     switch priorF
%         case 'gauss'
%             numPars = 14;
%         case 'mixgauss'
%             numPars = 19;
%     end
end

numUniqCont = numel(unique([cStim1 cStim2]));

% Save ID for display to terminal during parallel fitting
tempID = regexp(simDataPath,filesep,'split');
dispID = tempID{end};


%% Make anonymous function for negative log-likelihood & fit model to data

% Fix velocity-dependent component of likelihood widths
vfix = parVecGT(1:numS1Vels);

nllOpts = {'expon',fitmethod,interpOn,nllMethod,binScale};

lossFun = @(parVec)(calculate_model_nll([vfix parVec],ptvs_data,...
                    vStim1,cStim1,vStim2Delta,cStim2,numTrials,nllOpts));

% Fminunc/fmincon
opts    = optimset('display','off','useparallel',useParal,'tolx',1e-13,...
                   'maxfunevals',1e4,'largescale', 'off');
% lb = [-inf*ones(1,numS1Vels) zeros(1,numUniqCont)];
% ub = inf*ones(1,numS1Vels + numUniqCont);
lb = [0.001*ones(1,numUniqCont) -100*ones(1,numS1Vels)];
ub = 100*ones(1,numS1Vels + numUniqCont);

% Make sure likelihood width is a monotonic decreasing function of contrast
A = [-eye(numUniqCont-1) zeros(numUniqCont-1,1)] + [zeros(numUniqCont-1,1) eye(numUniqCont-1)];
% Pad with zeros to cover non-contrast pars
A = [A zeros(numUniqCont-1,numS1Vels)];
b = zeros(numUniqCont-1,1);


%% Fit data

indCnt  = 1;

avlogF      = nan(numReps,numS1Vels);
gvlogF      = nan(numReps,numS1Vels);
hcF         = nan(numReps,numUniqCont);
ptvsF       = cell(numReps,1);
nllF        = nan(numReps,1);
initConds   = nan(numReps,numPars);


% Loop over number of optimization runs requested
while indCnt <= numReps
    
    disp(['Run ',num2str(indCnt),'/',num2str(numReps),...
          ' - ',dispID,', n=',num2str(numTrials)]);
      
    tsl = tic;
    
    %%%%%%%%%%%%%% Choose Initialization %%%%%%%%%%%%%%%
    
   % Randomize initial parameters within reasonable window
   %%% grab preallocated vals to test binomial/bernoulli
   if isfield(initPars,'hcinit')
       avloginit       = initPars.avloginit(indCnt,:);
       hcinit          = initPars.hcinit(indCnt,:);
   else
       error('No pre allocated parameter vector');
   end
   
   parVec0         = [hcinit avloginit];
    
%      % Use GT for testing
%      if xfit
%
%          parVecInit  = [-(1/parVecGT(end)^2)*log(1 + vStim1/0.3) parVecGT(1:numS1Vels+numS2Conts)];
%      else
%          parVecInit  = parVecGT;
%      end
%      parVec0                     = parVecInit.*(rand(1,numPars) + 0.5);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    [prshat,~,exFlag,outStruc]  = fmincon(lossFun,parVec0,A,b,[],[],lb,ub,[],opts);
catch
    %%% put below back in once done comparing bernoulli and binomial
%     try
%         avloginit       = rand(1,numS1Vels)*(5-(-5)) +  -5;
%         hcinit          = rand(1,numS2Conts)*(2-0.1) + 0.1;
%         parVec0         = [avloginit hcinit];
%         [prshat,~,exFlag,outStruc]  = fmincon(lossFun,parVec0,[],[],[],[],lb,ub,[],opts);
%     catch
        % if rerolling the initial conditions fails twice in a row,
        % just give up.
        
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

    gvlogF(indCnt,:)  = vfix;
    hcF(indCnt,:)     = prshat(1:numUniqCont);
    avlogF(indCnt,:)  = prshat(numUniqCont + 1:end);

    initConds(indCnt,:) = [vfix hcinit avloginit];
    os(indCnt)          = outStruc;
    
    % Output best fitting psychometric function and nLL
    ptvsF{indCnt}     = calculate_ptvs(vStim1,cStim1,vStim2Delta,cStim2,...
                                       avlogF(indCnt,:),gvlogF(indCnt,:),hcF(indCnt,:),fitmethod,interpOn);
    nllF(indCnt)      = lossFun(prshat);
    
    % Unpack GT parameters
    if xfit
        nllOptsXfit   = {'mixgauss',fitmethod,interpOn,nllMethod,binScale};
        lossFunXfit   = @(parVec)(calculate_model_nll([vfix parVec],ptvs_data,...
                        vStim1,cStim1,vStim2Delta,cStim2,numTrials,nllOptsXfit));
                    
        nllGT(indCnt) = lossFunXfit(parVecGT(numS1Vels+1:end));
        
        gvlog         = parVecGT(1:numS1Vels);
        hc            = parVecGT(numS1Vels + 1:numS1Vels + numUniqCont);
        
        numComp       = (numel(parVecGT) - numS1Vels - numUniqCont)/3;
        
        sigP          = parVecGT(numS1Vels + numUniqCont + 1:numS1Vels + numUniqCont + numComp);
        muP           = parVecGT(numS1Vels + numUniqCont + numComp + 1:numS1Vels + numUniqCont + 2*numComp);
        wP            = parVecGT(numS1Vels + numUniqCont + 2*numComp + 1:numS1Vels + numUniqCont + 3*numComp);
        
        % Local slopes are just the weighted sum of the component slopes
        dNorm         = @(x,mu,sig) -(1/sig^2)*(x-mu);
        x             = getLogXform(vStim1,0.3);

        avlog         = zeros(1,numS1Vels);
        for ii = 1:numComp
            avlog     = avlog + wP(ii)*dNorm(x,muP(ii),sigP(ii));
        end
    else
        gvlog         = parVecGT(1:numS1Vels);
        hc            = parVecGT(numS1Vels + 1:numS1Vels + numUniqCont);
        avlog         = parVecGT(numS1Vels + numUniqCont + 1:end);
        freePars      = [hc avlog];
        
        nllGT(indCnt) = lossFun(freePars);
    end
    
%     if (exFlag > 0) || ~(strcmp(os(indCnt),'quasi-newton'))
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
    save(parEstPath,'numTrials','gvlog','hc','sigP','muP','wP','avlog','initConds','gvlogF',...
        'hcF','avlogF','ptvsF','ptvs_data','nllF','nllGT','os','parVecGT','estOpts','rsx');
    
    fitPars.gvlog       = gvlog;
    fitPars.hc          = hc;
    fitPars.sigP        = sigP;
    fitPars.muP         = muP;
    fitPars.wP          = wP;
else
    save(parEstPath,'numTrials','gvlog','hc','avlog','initConds','gvlogF',...
        'hcF','avlogF','ptvsF','ptvs_data','nllF','nllGT','os','parVecGT','estOpts','rsx');

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
