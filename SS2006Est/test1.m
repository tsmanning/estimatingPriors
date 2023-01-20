% Set number of initial starting positions for parameter estimation
global numReps
numReps  = 1;

% Set number of times to bootstrap
numBoot  = 1;

% Use parallel comp toolbox?
global useParal
useParal = false;

% Define MoG method/number of components
global numComp
numComp  = 1;

% Set whether or not to run fitting
runFit = true;
% runFit = false;

% Fitting verbosity
global verbo
verbo = 'iter';
% verbo = 'off';

splPath = regexp(which('MogVsPiece'),filesep,'split');
fDir    = [filesep,fullfile(splPath{1:numel(splPath)-2}),filesep];
sDir    = [fDir,'SS2006Data',filesep];

load([sDir,'s1']);

[A1,A2,A3,A4] = piece_est_noGV(s1);

% [s1_priorPW,vels] = buildPiecewisePrior(slopes_s1);
[s1_priorPW,vels] = recoverPriorF_SSData2(A1,A3,A2,A4);


%% helper functions

function [bestPrSlopes,bestPrGv,bestPrHc,bestNll,outMat] = piece_est_noGV(dataMat)

% Setup
global numReps
global useParal
global verbo

cStim2      = unique( dataMat(4,:) );
numUniqCont = numel(cStim2);
vStim1      = unique( dataMat(1,:) );
numS1Vels   = numel(vStim1);
numPars     = numUniqCont + numS1Vels;

gvFix   = 0.25;
GTgvlog = gvFix*ones(1,numel(vStim1));

% Make anonymous function for negative log-likelihood & fit model to data
method  = 'piece';
% lossFun = @(parVec)(calcPFxn_piece([gvFix*ones(1,numS1Vels) parVec],dataMat));
lossFun = @(parVec)(cNLL_pw(parVec,dataMat,vStim1,cStim2,GTgvlog));

% Setup fmincon constraints               
opts    = optimset('useparallel',useParal,'Algorithm', 'interior-point', 'Display', verbo, ...
        'MaxFunEvals', 5000, 'MaxIter', 500, 'GradObj', 'off');

% lb = [0.001*ones(1,numUniqCont) ...
%       -100*ones(1,numS1Vels)];
% ub = [100*ones(1,numUniqCont) ...
%       100*ones(1,numS1Vels)];
lb = [-100*ones(1,numS1Vels) ...
      0.001*ones(1,numUniqCont)];
ub = [100*ones(1,numS1Vels) ...
      100*ones(1,numUniqCont)];
  
% Make sure likelihood width is a monotonic decreasing function of contrast
% A = [-eye(numUniqCont-1) zeros(numUniqCont-1,1)] + [zeros(numUniqCont-1,1) eye(numUniqCont-1)];
% Pad with zeros to cover non-contrast pars
% A = [A zeros(numUniqCont-1,numS1Vels)];
% b = zeros(numUniqCont-1,1);

% Generate random set of initial conditions
% offsetSlope        = repmat(0.1-0.01*(1:numUniqCont),[numReps 1]);
% initPars.hcinit    = rand(numReps,numUniqCont)*(2-0.1) + 0.1 + offsetSlope;
% initPars.avloginit = rand(numReps,numS1Vels)*(5-(-5))  + -5;

% Lock to Ben's starting points each time
initPars.hcinit    = 1.2*ones(numReps,numUniqCont);
initPars.avloginit = -2.2*ones(numReps,numS1Vels);

% Estimation loop

indCnt  = 1;

avlogF      = nan(numReps,numS1Vels);
gvlogF      = nan(numReps,numS1Vels);
hcF         = nan(numReps,numUniqCont);
nllF        = nan(numReps,1);
initConds   = nan(numReps,numPars);

while indCnt <= numReps
    
    avloginit       = initPars.avloginit(indCnt,:);
    hcinit          = initPars.hcinit(indCnt,:);
    
%     parVec0         = [hcinit avloginit];
    parVec0         = [avloginit hcinit];

    try
%         [prshat,~,~,outStruc]  = fmincon(lossFun,parVec0,A,b,[],[],lb,ub,[],opts);
        [prshat,~,~,outStruc]  = fmincon(lossFun,parVec0,[],[],[],[],lb,ub,[],opts);
    catch
        prshat                   = nan(1,numel(parVec0));
        outStruc.iterations      = [];
        outStruc.funcCount       = [];
        outStruc.constrviolation = [];
        outStruc.stepsize        = [];
        outStruc.algorithm       = [];
        outStruc.firstorderopt   = [];
        outStruc.cgiterations    = [];
        outStruc.message         = 'failed to initialize pars';
    end
    
    gvlogF(indCnt,:)  = gvFix*ones(1,numS1Vels);
%     hcF(indCnt,:)     = prshat(1:numUniqCont);
%     avlogF(indCnt,:)  = prshat(numUniqCont + 1:end);
    avlogF(indCnt,:)  = prshat(1:numS1Vels);
    hcF(indCnt,:)     = prshat(numS1Vels + 1:end);
    
    initConds(indCnt,:) = parVec0;
    os(indCnt)          = outStruc;
    
    % Output best fitting psychometric function and nLL
    nllF(indCnt)      = lossFun(prshat);
    
    indCnt = indCnt + 1;
    
end

% Select best fit & nll to output
[~,bestInd]      = min(nllF);
bestNll          = nllF(bestInd);
bestPrSlopes     = avlogF(bestInd,:);
bestPrHc         = hcF(bestInd,:);
bestPrGv         = gvlogF(bestInd,:);

outMat.avlogF    = avlogF;
outMat.gvlogF    = gvlogF;
outMat.hcF       = hcF;
outMat.nllF      = nllF;
outMat.initConds = initConds;

end

function [nll] = calcModelNLL(pars,method,dataMat)

% Get 2AFC responses
choices = dataMat(9,:);

numV1 = numel(unique( dataMat(1,:) ));
numC2 = numel(unique( dataMat(4,:) ));

% Calculate p("v2">"v1") for current model
switch method
    case 'piece'
        sPars{1} = pars(1                :numV1);            % gvlog
%         sPars{2} = pars(numV1 + 1        :numV1 + numC2);    % hc
%         sPars{3} = pars(numV1 + numC2 + 1:2*numV1 + numC2);  % avlog
        sPars{3} = pars(numV1 + 1        :2*numV1);    % hc
        sPars{2} = pars(2*numV1 + 1:2*numV1 + numC2);  % avlog
        
        [pFxn] = calcPFxn_piece(sPars,dataMat);
        
    case 'mog'
        % Assuming comps have all their parameters free
        numComp = (numel(pars) - numV1 - numC2)/3;
        
        sPars{1} = pars(1                            :numV1);                     % gvlog
        sPars{2} = pars(numV1 + 1                    :numV1 + numC2);             % hc
        sPars{3} = pars(numV1 + numC2 + 1            :numV1 + numC2 + numComp);   % w
        sPars{4} = pars(numV1 + numC2 + numComp + 1  :numV1 + numC2 + 2*numComp); % sig
        sPars{5} = pars(numV1 + numC2 + 2*numComp + 1:end);                       % mu
        
        [pFxn] = calcPFxn_mog(sPars,dataMat);
        
end

% Don't let likelihood go to zero
minlikli = eps;

% Compute log-likelihood of data given current model for each trial
% nllVec = -( choices.*log(pFxn) + (1-choices).*log(1-pFxn) );
% nll    = sum(max(nllVec,minlikli));
nll = -sum( choices.*log(pFxn) + (1-choices).*log(1-pFxn) );

end

function [nll] = calcPFxn_piece(pars,dataMat)

% Calculate psych function values given a set of observer/exp parameters

% Extract pars
% gvlog = pars{1};
% hc    = pars{2};
% avlog = pars{3};

gvlog = pars(1:6);
hc    = pars(13:19);
avlog = pars(7:12);

vStim1 = dataMat(1,:);
cStim1 = dataMat(2,:);
vStim2 = dataMat(3,:);
cStim2 = dataMat(4,:);

% Get logxforms
rf_vlog   = log(1 + vStim1/0.3);
test_vlog = log(1 + vStim2/0.3);

% Index into contrast pars
[~,~,tmpCntVec] = unique([cStim1;cStim2]);
tmpCntMat       = reshape(tmpCntVec,[2 size(dataMat,2)]);
tsCont          = tmpCntMat(1,:);
rfCont          = tmpCntMat(2,:);

% Index into velocity pars; use nearest neighbor Ref vel to index Test vel
[uniqVels,~,rfVel] = unique(vStim1);
% distVals  = abs(repmat(uniqVels',[1 size(dataMat,2)]) - ...
%                 repmat(vStim2,[numel(uniqVels) 1]));
% [~,tsVel] = min(distVals);

% Just make test index equal reference like scratch script
tsVel = rfVel;
    
% Grab likelihood sigmas from exp parameters
thisRSigL    = gvlog(rfVel).*hc(rfCont);
thisTSigL    = gvlog(tsVel).*hc(tsCont);

% Use nearest neighbor Ref vel to index into fitted slopes
refSlope  = avlog(rfVel);
testSlope = avlog(tsVel);

% Calculate p(v2_hat>v1_hat) with piecewise method (using linear
% interpolation)
rfMu    = rf_vlog + refSlope.*thisRSigL.^2;
rfVar   = thisRSigL.^2;

testMu  = test_vlog + testSlope.*thisTSigL.^2;
testVar = thisTSigL.^2;

Z    = (testMu - rfMu) ./ sqrt(rfVar + testVar);

% pFxn = min(max(normcdf(Z), eps), 1 - eps);
pFxn = normcdf(Z);

nll = -sum(dataMat(9,:) .* log(pFxn) + (1 - dataMat(9,:)) .* log(1 - pFxn));

end

function nll = cNLL_pw(parVec,subjData,vStim1,cStim2,GTgvlog)

%% Likelihood of data given current model pars

avlog = parVec(1:numel(vStim1));
hc    = parVec(numel(vStim1) + 1:end);
gvlog = GTgvlog;

% moving test slope outside of for loop for speed when using
% interpolation. Actually, nothing else
% here really needs a for loop, just used it for clarity when debugging

% Define test velocities
testVelLin    = subjData(3,:);

% Convert to normalized log space
testVelLog  = log(1 + testVelLin/0.3);

% for each trial, calculate mus and vars of map sampling for ref and
% test
for x = 1:size(subjData,2)
    
    %% reference
    rfCont = subjData(2,x);
    
    % Convert to normalized vlog space
    rfVelLin(x) = subjData(1,x);
    rfVelLog(x)  = log(1 + rfVelLin(x)/0.3);
    
    % mean and var of standard posterior in tranformed velocity space
    % hc = all contrasts; assume numel(cStim2) == numel(hc) and
    % cStim1 is a subset of cStim2
    refSlope(x) = avlog(rfVelLin(x) == vStim1);
    
    refMu(x)   = rfVelLog(x) + refSlope(x)*(gvlog(rfVelLin(x) == vStim1)*hc(rfCont == cStim2))^2;
    refVar(x)  = (gvlog(rfVelLin(x) == vStim1)*hc(rfCont == cStim2))^2;
    
    %% test
    testCont = subjData(4,x);
    testSlope(x) = avlog(rfVelLin(x) == vStim1);
    
    testMu(x)      = testVelLog(x) + testSlope(x)*(gvlog(rfVelLin(x) == vStim1)*hc(testCont == cStim2))^2;
    testVar(x)     = (gvlog(rfVelLin(x) == vStim1)*hc(testCont == cStim2))^2;
    
end

% naeker way
Z = (testMu - refMu) ./ sqrt(refVar + testVar);

rho = normcdf(Z);

% Compute log-likelihood - log of bernoulli distribution
nll = -sum(subjData(9,:) .* log(rho) + ...
    (1 - subjData(9,:)) .* log(1 - rho));

end
