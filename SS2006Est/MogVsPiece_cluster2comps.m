%% General settings

% Set number of initial starting positions for parameter estimation
global numReps
numReps  = 15;

% Set number of times to bootstrap
numBoot  = 15;

% Use parallel comp toolbox?
global useParal
% useParal = false;
useParal = true;

% Define MoG method/number of components
global numComp
numComp  = 2;

% Set whether or not to run fitting
% runpwFit = true;
% runmogFit = true;

% Fitting verbosity
global verbo
% verbo = 'iter';
verbo = 'off';

% Set resampling method
% rsMethod = 'simple';
rsMethod = 'triplet';

% Set plot on or off
plotOn = true;

% Assign multiple slopes per ref or use nearest neighbor
global slAs
% slAs = 'mspf';
slAs = 'nn';

% Set save string
ds = datestr(now,30);
% ds = '2021-08-16_moreps';
% saveStrPW = ['_',ds,'_',slAs,'_',rsMethod];
% saveStrPW = ['_','2021-07-27','_',slAs,'_',rsMethod];
saveStr   = ['_',ds,'_',num2str(numComp),'comps_',slAs,'_',rsMethod];

%% Load in Data & setup

% Their data arrangement:
% 11 x n mat where n is the number of trials
% Row 1: Ref Speed (deg/s)
% Row 2: Ref Cont (%)
% Row 3: Test Speed (deg/s)
% Row 4: Test Cont (%)
% Row 9: Response (Test>Ref? - 0/1)

splPath = regexp(which('MogVsPiece'),filesep,'split');
fDir    = [filesep,fullfile(splPath{1:numel(splPath)-2}),filesep];
sDir    = [fDir,'SS2006Data',filesep];
figDir    = [fDir,'SS2006Data/figures',filesep];

load([sDir,'s1']);
load([sDir,'s2']);

%%%% should prob fix this later to not assume same for both subs?
V1s = unique( s1(1,:) );
numV1 = numel(V1s);
C1s = unique( s1(2,:) );
numC1 = numel(C1s);
C2s = unique( s1(4,:) );
numC2 = numel(C2s);

% Add dependencies
addpath([fDir,'0-SharedTools']);
addpath([fDir,'4-MixtureGaussian']);

rng('shuffle');

%% Resample for bootstrap
switch rsMethod
    
    case 'simple'
        % Simple resampling with replacement
        for ii = 1:numBoot
            s1boot(ii).s1 = datasample(s1,size(s1,2),2,'Replace',true);
            s2boot(ii).s2 = datasample(s2,size(s2,2),2,'Replace',true);
        end
        
    case 'triplet'
        % Resampling within a [V1,C1,C2] triplet
        for ii = 1:numBoot
            
            s1boot(ii).s1 = [];
            s2boot(ii).s2 = [];
            
            for jj = 1:numV1
                for kk = 1:numC1
                    for ll = 1:numC2
                        
                        theseIndsS1 = (s1(1,:) == V1s(jj)) & ...
                                      (s1(2,:) == C1s(kk)) & ...
                                      (s1(4,:) == C2s(ll));
                        
                        theseIndsS2 = (s2(1,:) == V1s(jj)) & ...
                                      (s2(2,:) == C1s(kk)) & ...
                                      (s2(4,:) == C2s(ll));
                        
                        rsTrialss1    = datasample(s1(:,theseIndsS1),size(s1(:,theseIndsS1),2),2,'Replace',true);
                        rsTrialss2    = datasample(s2(:,theseIndsS2),size(s2(:,theseIndsS2),2),2,'Replace',true);
                        
                        s1boot(ii).s1 = [s1boot(ii).s1 rsTrialss1];
                        s2boot(ii).s2 = [s2boot(ii).s2 rsTrialss2];
                        
                    end
                end
            end
            
        end
        
        save([sDir,'fitData',saveStr,'bootData'],'s1boot','s2boot');
        disp(['Saving bootstrapped data to: ','./fitData',saveStr,'bootData.mat']);
end

%% MoG

nCores = feature('numcores');
maxNumCompThreads(nCores);
parpool(nCores);

% True
% Run first subject
disp('Starting first sub (MoG)');
tic
[bestSigPFMoG_s1,bestMuPFMoG_s1,bestWPFMoG_s1,...
    bestGvMoG_s1,bestHcMoG_s1,bestNllMoG_s1,outMat_mog1] = MoG_est_OnlySig(s1);
disp('Finished first sub (MoG)');
toc

% Run second subject
disp('Starting second sub (MoG)');
tic
[bestSigPFMoG_s2,bestMuPFMoG_s2,bestWPFMoG_s2,...
    bestGvMoG_s2,bestHcMoG_s2,bestNllMoG_s2,outMat_mog2] = MoG_est_OnlySig(s2);
disp('Finished second sub (MoG)');
toc

% Resampled
s1mogsig = nan(numBoot,numComp);
s1mogmu  = nan(numBoot,numComp);
s1mogw   = nan(numBoot,numComp);
s1moggv  = nan(numBoot,numV1);
s1moghc  = nan(numBoot,numC2);
s1mognll = nan(numBoot,1);
s2mogsig = nan(numBoot,numComp);
s2mogmu  = nan(numBoot,numComp);
s2mogw   = nan(numBoot,numComp);
s2moggv  = nan(numBoot,numV1);
s2moghc  = nan(numBoot,numC2);
s2mognll = nan(numBoot,1);

for ii = 1:numBoot
    disp(['First Sub Bootstrap (MoG) ',num2str(ii),'/',num2str(numBoot)]);
    [s1mogsig(ii,:),s1mogmu(ii,:),s1mogw(ii,:),s1moggv(ii,:),s1moghc(ii,:),s1mognll(ii,:)] = MoG_est_OnlySig(s1boot(ii).s1);
    disp(['Second Sub Bootstrap (MoG) ',num2str(ii),'/',num2str(numBoot)]);
    [s2mogsig(ii,:),s2mogmu(ii,:),s2mogw(ii,:),s2moggv(ii,:),s2moghc(ii,:),s2mognll(ii,:)] = MoG_est_OnlySig(s2boot(ii).s2);
end


mogDat.bestSigPFMoG_s1  = bestSigPFMoG_s1;
mogDat.bestMuPFMoG_s1   = bestMuPFMoG_s1;
mogDat.bestWPFMoG_s1    = bestWPFMoG_s1;
mogDat.bestGvMoG_s1     = bestGvMoG_s1;
mogDat.bestHcMoG_s1     = bestHcMoG_s1;
mogDat.bestNllMoG_s1    = bestNllMoG_s1;
mogDat.outMat_mog1      = outMat_mog1;

mogDat.bestSigPFMoG_s2  = bestSigPFMoG_s2;
mogDat.bestMuPFMoG_s2   = bestMuPFMoG_s2;
mogDat.bestWPFMoG_s2    = bestWPFMoG_s2;
mogDat.bestGvMoG_s2     = bestGvMoG_s2;
mogDat.bestHcMoG_s2     = bestHcMoG_s2;
mogDat.bestNllMoG_s2    = bestNllMoG_s2;
mogDat.outMat_mog2      = outMat_mog2;

mogDat.s1mogsig = s1mogsig;
mogDat.s1mogmu  = s1mogmu;
mogDat.s1mogw   = s1mogw;
mogDat.s1moggv  = s1moggv;
mogDat.s1moghc  = s1moghc;
mogDat.s1mognll = s1mognll;
mogDat.s2mogsig = s2mogsig;
mogDat.s2mogmu  = s2mogmu;
mogDat.s2mogw   = s2mogw;
mogDat.s2moggv  = s2moggv;
mogDat.s2moghc  = s2moghc;
mogDat.s2mognll = s2mognll;

save([sDir,'fitData',saveStr,'mogDat'],'mogDat');
disp(['Saving fitted parameters to: ','./fitData',saveStr,'mogDat.mat']);

%% Helper functions

function [bestSigPF,bestMuPF,bestWPF,bestGv,bestHc,bestNll,outMat] = MoG_est(dataMat)

% Setup
global numReps
global useParal
global numComp
global verbo

numUniqCont = numel(unique( dataMat(4,:) ));
numS1Vels   = numel(unique( dataMat(1,:) ));
numPars     = numUniqCont + numS1Vels + 3*numComp;

% Make anonymous function for negative log-likelihood & fit model to data
method  = 'mog';
lossFun = @(parVec)(calcModelNLL(parVec,method,dataMat));
                    
% Fminunc/fmincon
opts    = optimset('display',verbo,'useparallel',useParal,'tolx',1e-13,...
                   'maxfunevals',1e4,'largescale', 'off');
               
% parVec = [gv hc w sig mu];
               
% Best fit sigma must be at least as big as the smallest interval between test speeds
smallInt = min(diff( getLogXform(unique(dataMat(3,:)),0.3) ));
% lb = [0.001*ones(1,numS1Vels) ...
%       0.001*ones(1,numUniqCont) ...
%       zeros(1,numComp) ...
%       smallInt*ones(1,numComp) ...
%       -eps*ones(1,numComp)];
% ub = [100*ones(1,numS1Vels) ...
%       100*ones(1,numUniqCont) ...
%       ones(1,numComp) ...
%       1e5*ones(1,numComp) ...
%       100*ones(1,numComp)];    

% Clamp mu to only fit sigma
lb = [0.001*ones(1,numS1Vels) ...
      0.001*ones(1,numUniqCont) ...
      zeros(1,numComp) ...
      smallInt*ones(1,numComp) ...
      -2*eps*ones(1,numComp)];
ub = [100*ones(1,numS1Vels) ...
      100*ones(1,numUniqCont) ...
      ones(1,numComp) ...
      1e5*ones(1,numComp) ...
      2*eps*ones(1,numComp)];   
               
% All weights must sum to 1
Aeq = [zeros(1,numS1Vels + numUniqCont) ones(1,numComp) zeros(1,2*numComp)];
beq = 1;

% Make sure likelihood width is a monotonic decreasing function of contrast
A = [-eye(numUniqCont-1) zeros(numUniqCont-1,1)] + [zeros(numUniqCont-1,1) eye(numUniqCont-1)];
A = [zeros(numUniqCont-1,numS1Vels) A zeros(numUniqCont-1,3*numComp)];
b = zeros(numUniqCont-1,1);

% Generate random set of initial conditions
offsetSlope        = repmat(0.1-0.01*(1:numUniqCont),[numReps 1]);
initPars.hcinit    = rand(numReps,numUniqCont)*(2-0.1) + 0.1 + offsetSlope;
initPars.gvloginit = rand(numReps,numS1Vels)*(0.25) +  0.3;
initPars.sigPinit  = rand(numReps,numComp)*(1-0.1)   + 0.1;
% initPars.muPinit   = rand(numReps,numComp)*5         + 0.0;
% initPars.wPinit    = rand(numReps,numComp);

% Only fit sigma
initPars.muPinit   = zeros(numReps,numComp);
initPars.wPinit    = ones(numReps,numComp);

% Estimation loop

indCnt  = 1;

gvlogF      = nan(numReps,numS1Vels);
hcF         = nan(numReps,numUniqCont);
sigPF       = nan(numReps,numComp);
muPF        = nan(numReps,numComp);
wPF         = nan(numReps,numComp);
% ptvsF       = cell(1,numReps);
nllF        = nan(numReps,1);
initConds   = nan(numReps,numPars);

while indCnt <= numReps

    % Randomize parameter initialization and package up
    gvinit          = initPars.gvloginit(indCnt,:);
    hcinit          = initPars.hcinit(indCnt,:);
    sigPinit        = initPars.sigPinit(indCnt,:);
    muPinit         = initPars.muPinit(indCnt,:);
    wPinit          = initPars.wPinit(indCnt,:);
    
    parVec0         = [gvinit hcinit wPinit sigPinit muPinit];
    
%     try
        [prshat,~,~,outStruc]  = fmincon(lossFun,parVec0,A,b,Aeq,beq,lb,ub,[],opts);
%     catch
%         prshat                   = nan(1,numel(parVec0));
%         outStruc.iterations      = [];
%         outStruc.funcCount       = [];
%         outStruc.constrviolation = [];
%         outStruc.stepsize        = [];
%         outStruc.algorithm       = [];
%         outStruc.firstorderopt   = [];
%         outStruc.cgiterations    = [];
%         outStruc.message         = 'failed to initialize pars';
%     end

    gvlogF(indCnt,:) = prshat(1                                      :numS1Vels);
    hcF(indCnt,:)    = prshat(numS1Vels + 1                          :numS1Vels + numUniqCont);
    
    wPF(indCnt,:)    = prshat(numS1Vels + numUniqCont + 1            :numS1Vels + numUniqCont + numComp); 
    sigPF(indCnt,:)  = prshat(numS1Vels + numUniqCont + numComp + 1  :numS1Vels + numUniqCont + 2*numComp);
    muPF(indCnt,:)   = prshat(numS1Vels + numUniqCont + 2*numComp + 1:end);
    
    initConds(indCnt,:) = parVec0;
    os(indCnt)          = outStruc;
    
%     ptvsF{indCnt}   = calcPFxn_mog(prshat,dataMat);
    nllF(indCnt)    = lossFun(prshat);
    
    indCnt = indCnt + 1;

end

% Select best fit & nll to output
[~,bestInd]      = min(nllF);
bestNll          = nllF(bestInd);
bestSigPF        = sigPF(bestInd,:);
bestMuPF         = muPF(bestInd,:);
bestWPF          = wPF(bestInd,:);
bestGv           = gvlogF(bestInd,:);
bestHc           = hcF(bestInd,:);
% bestPars         = pFxn(bestInd);

outMat.gvlogF    = gvlogF;
outMat.hcF       = hcF;
outMat.sigPF     = sigPF;
outMat.muPF      = muPF;
outMat.wPF       = wPF;
% outMat.ptvsF     = pFxn;
outMat.nllF      = nllF;
outMat.initConds = initConds;

end

function [bestSigPF,bestMuPF,bestWPF,bestGv,bestHc,bestNll,outMat] = MoG_est_OnlySig(dataMat)

% Setup
global numReps
global useParal
global numComp
global verbo

numUniqCont = numel(unique( dataMat(4,:) ));
numS1Vels   = numel(unique( dataMat(1,:) ));
numPars     = numUniqCont + numComp;

gvFix = 0.25;
nonlinOn = true;

% Make anonymous function for negative log-likelihood & fit model to data
method  = 'mog';
lossFun = @(parVec)(calcModelNLL([gvFix*ones(1,numS1Vels) ...
                                  parVec(1:numUniqCont) ...
                                  ones(1,numComp)/numComp ...
                                  parVec(numUniqCont + 1:numUniqCont + numComp) ...
                                  zeros(1,numComp)],method,dataMat));
                    
% Fminunc/fmincon
opts = optimset('useparallel',useParal,'Algorithm', 'interior-point', 'Display', verbo, ...
                'MaxFunEvals', 5000, 'MaxIter', 500, 'GradObj', 'off');
                           
% Best fit sigma must be at least as big as the smallest interval between test speeds
smallInt = min(diff( getLogXform(unique(dataMat(3,:)),0.3) ));
  
% Only fit sigma
lb = [0.001*ones(1,numUniqCont) smallInt*ones(1,numComp)];
ub = [100*ones(1,numUniqCont) 1e5*ones(1,numComp)];

% Make sure likelihood width is a monotonic decreasing function of contrast
if nonlinOn
A = [-eye(numUniqCont-1) zeros(numUniqCont-1,1)] + [zeros(numUniqCont-1,1) eye(numUniqCont-1)];
A = [A zeros(numUniqCont-1,numComp)];
b = zeros(numUniqCont-1,1);
end

% Generate random set of initial conditions
% offsetSlope        = repmat(0.1-0.01*(1:numUniqCont),[numReps 1]);
% initPars.hcinit    = rand(numReps,numUniqCont)*(2-0.1) + 0.1 + offsetSlope;
% initPars.sigPinit  = rand(numReps,numComp)*(1-0.1)   + 0.1;

% Estimation loop

indCnt  = 1;

gvlogF      = nan(numReps,numS1Vels);
hcF         = nan(numReps,numUniqCont);
sigPF       = nan(numReps,numComp);
muPF        = nan(numReps,numComp);
wPF         = nan(numReps,numComp);
nllF        = nan(numReps,1);
initConds   = nan(numReps,numPars);

while indCnt <= numReps

    % Randomize parameter initialization and package up
%     hcinit          = initPars.hcinit(indCnt,:);
%     sigPinit        = initPars.sigPinit(indCnt,:);
    offsetSlope        = 0.1-0.01*(1:numUniqCont);
    hcinit    = rand(1,numUniqCont)*(2-0.1) + 0.1 + offsetSlope;
    sigPinit  = rand(1,numComp)*(1-0.1)   + 0.1;
    
    parVec0   = [hcinit sigPinit];
    
    try
        if nonlinOn
            [prshat,~,ef,outStruc]  = fmincon(lossFun,parVec0,A,b,[],[],lb,ub,[],opts);
        else
            [prshat,~,ef,outStruc]  = fmincon(lossFun,parVec0,[],[],[],[],lb,ub,[],opts);
        end
    catch
        ef = 0;
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

    hcF(indCnt,:)    = prshat(1:numUniqCont);
    sigPF(indCnt,:)  = prshat(numUniqCont + 1:numUniqCont + numComp);
    
    gvlogF(indCnt,:) = gvFix*ones(1,numS1Vels);
    wPF(indCnt,:)    = ones(1,numComp)/numComp; 
    muPF(indCnt,:)   = zeros(1,numComp);
    
    initConds(indCnt,:) = parVec0;
    os(indCnt)          = outStruc;
    
    nllF(indCnt)    = lossFun(prshat);
    
    if ef ~= 0
        indCnt = indCnt + 1;
    end
    
end

% Select best fit & nll to output
[~,bestInd]      = min(nllF);
bestNll          = nllF(bestInd);
bestSigPF        = sigPF(bestInd,:);
bestMuPF         = muPF(bestInd,:);
bestWPF          = wPF(bestInd,:);
bestGv           = gvlogF(bestInd,:);
bestHc           = hcF(bestInd,:);

outMat.gvlogF    = gvlogF;
outMat.hcF       = hcF;
outMat.sigPF     = sigPF;
outMat.muPF      = muPF;
outMat.wPF       = wPF;
outMat.nllF      = nllF;
outMat.initConds = initConds;

end

function [bestPrSlopes,bestPrGv,bestPrHc,bestNll,outMat] = piece_est(dataMat)

% Setup
global numReps
global useParal
global verbo

numUniqCont = numel(unique( dataMat(4,:) ));
numS1Vels   = numel(unique( dataMat(1,:) ));
numPars     = numUniqCont + 2*numS1Vels;

% Make anonymous function for negative log-likelihood & fit model to data
method  = 'piece';
lossFun = @(parVec)(calcModelNLL(parVec,method,dataMat));

% Setup fmincon constraints
opts    = optimset('display',verbo,'useparallel',useParal,'tolx',1e-13,...
                   'maxfunevals',1e4,'largescale', 'off');

lb = [0.001*ones(1,numS1Vels) ...
      0.001*ones(1,numUniqCont) ...
      -100*ones(1,numS1Vels)];
ub = [100*ones(1,numS1Vels) ...
      100*ones(1,numUniqCont) ...
      100*ones(1,numS1Vels)];
% ub = [100*ones(1,numS1Vels) ...
%       100*ones(1,numUniqCont) ...
%       0*ones(1,numS1Vels)];

% Make sure likelihood width is a monotonic decreasing function of contrast
A = [-eye(numUniqCont-1) zeros(numUniqCont-1,1)] + [zeros(numUniqCont-1,1) eye(numUniqCont-1)];
% Pad with zeros to cover non-contrast pars
A = [zeros(numUniqCont-1,numS1Vels) A zeros(numUniqCont-1,numS1Vels)];
b = zeros(numUniqCont-1,1);

% Generate random set of initial conditions
offsetSlope        = repmat(0.1-0.01*(1:numUniqCont),[numReps 1]);
initPars.hcinit    = rand(numReps,numUniqCont)*(2-0.1) + 0.1 + offsetSlope;
initPars.gvloginit = rand(numReps,numS1Vels)*(0.25)    + 0.3;
initPars.avloginit = rand(numReps,numS1Vels)*(5-(-5))  + -5;

% Estimation loop

indCnt  = 1;

avlogF      = nan(numReps,numS1Vels);
gvlogF      = nan(numReps,numS1Vels);
hcF         = nan(numReps,numUniqCont);
% pFxn_PW        = cell(1,numReps);
nllF        = nan(numReps,1);
initConds   = nan(numReps,numPars);

while indCnt <= numReps
    
    gvloginit       = initPars.gvloginit(indCnt,:);
    avloginit       = initPars.avloginit(indCnt,:);
    hcinit          = initPars.hcinit(indCnt,:);
    
    parVec0         = [gvloginit hcinit avloginit];
    
    try
        [prshat,~,~,outStruc]  = fmincon(lossFun,parVec0,A,b,[],[],lb,ub,[],opts);
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
    
    gvlogF(indCnt,:)  = prshat(1:numS1Vels);
    hcF(indCnt,:)     = prshat(numS1Vels + 1:numS1Vels + numUniqCont);
    avlogF(indCnt,:)  = prshat(numS1Vels + numUniqCont + 1:end);

    initConds(indCnt,:) = parVec0;
    os(indCnt)          = outStruc;
    
    % Output best fitting psychometric function and nLL
%     pFxn{indCnt}      = calcPFxn_piece(prshat,dataMat);
    nllF(indCnt)      = lossFun(prshat);
    
    indCnt = indCnt + 1;
    
end

% Select best fit & nll to output
[~,bestInd]      = min(nllF);
bestNll          = nllF(bestInd);
bestPrSlopes     = avlogF(bestInd,:);
bestPrGv         = gvlogF(bestInd,:);
bestPrHc         = hcF(bestInd,:);

outMat.avlogF    = avlogF;
outMat.gvlogF    = gvlogF;
outMat.hcF       = hcF;
% outMat.ptvsF     = pFxn;
outMat.nllF      = nllF;
outMat.initConds = initConds;

end

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
% lossFun = @(parVec)(calcModelNLL([gvFix*ones(1,numS1Vels) parVec],method,dataMat));
lossFun = @(parVec)(cNLL_pw(parVec,dataMat,vStim1,cStim2,GTgvlog));

% Setup fmincon constraints               
opts    = optimset('useparallel',useParal,'Algorithm', 'interior-point', 'Display', verbo, ...
        'MaxFunEvals', 5000, 'MaxIter', 500, 'GradObj', 'off');

lb = [0.001*ones(1,numUniqCont) ...
      -100*ones(1,numS1Vels)];
ub = [100*ones(1,numUniqCont) ...
      100*ones(1,numS1Vels)];
% lb = [-100*ones(1,numS1Vels) ...
%       0.001*ones(1,numUniqCont)];
% ub = [100*ones(1,numS1Vels) ...
%       100*ones(1,numUniqCont)];

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
initPars.hcinit    = 1*ones(numReps,numUniqCont);
initPars.avloginit = -2*ones(numReps,numS1Vels);

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
    
    parVec0         = [hcinit avloginit];
%     parVec0         = [avloginit hcinit];

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
    hcF(indCnt,:)     = prshat(1:numUniqCont);
    avlogF(indCnt,:)  = prshat(numUniqCont + 1:end);
%     avlogF(indCnt,:)  = prshat(1:numS1Vels);
%     hcF(indCnt,:)     = prshat(numS1Vels + 1:end);

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
        sPars{2} = pars(numV1 + 1        :numV1 + numC2);    % hc
        sPars{3} = pars(numV1 + numC2 + 1:2*numV1 + numC2);  % avlog
       
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

function [pFxn] = calcPFxn_piece(pars,dataMat)

% Calculate psych function values given a set of observer/exp parameters

% Extract pars
gvlog = pars{1};
hc    = pars{2};
avlog = pars{3};

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

end

function [pFxn] = calcPFxn_mog(pars,dataMat)

% Calculate psych function values given a set of observer/exp parameters

% Assumes dataMat is a mxn matrix where first 4 rows hold stimulus values
% and n is number of trials

numTrials = size(dataMat,2);

% Extract pars (assuming all components free version)
gvlog = pars{1};
hc    = pars{2};
wP    = pars{3};
sigP  = pars{4};
muP   = pars{5};

vStim1 = dataMat(1,:);
cStim1 = dataMat(2,:);
vStim2 = dataMat(3,:);
cStim2 = dataMat(4,:);

% Get logxforms
rf_vlog   = log(1 + vStim1/0.3);
test_vlog = log(1 + vStim2/0.3);

% Assign index to each unique contrast
[~,~,tmpCntVec] = unique([cStim1;cStim2]);
tmpCntMat       = reshape(tmpCntVec,[2 size(dataMat,2)]);
rfCont          = tmpCntMat(1,:);
tsCont          = tmpCntMat(2,:);

% Assign index to each unique velocity; use nearest neighbor Ref vel to 
% index Test vel
[uniqVels,~,rfVel] = unique(vStim1);
distVals  = abs(repmat(uniqVels',[1 size(dataMat,2)]) - ...
                repmat(vStim2,[numel(uniqVels) 1]));
[~,tsVel] = min(distVals);

% Grab likelihood sigmas from exp parameters
thisRSigL    = gvlog(rfVel).*hc(rfCont);
thisTSigL    = gvlog(tsVel).*hc(tsCont);

% Calculate p(v2_hat>v1_hat) with MoG estimation
% MoG post pars are (nxm) where n = # components, m = # likelihoods
[alphas1,muTildes1,wTildes1] = getMogPostPars(wP,sigP,muP,thisRSigL,rf_vlog);
[alphas2,muTildes2,wTildes2] = getMogPostPars(wP,sigP,muP,thisTSigL,test_vlog);

% Loop over trials (or figure out how to reorg pars, since you'd need a 3D
% mat)
pFxn = nan(1,numTrials);

% estDistEval = 'original';
estDistEval = 'rederived';

for ii = 1:numTrials
    switch estDistEval
        case 'original'
            denom   = sqrt((thisTSigL(ii)^2)*(alphas2(:,ii).^2) + (thisRSigL(ii)^2)*(alphas1(:,ii).^2));
            
            Z       = (muTildes2(:,ii)  + alphas2(:,ii)*test_vlog(ii) ...
                     - muTildes1(:,ii)' - alphas1(:,ii)'*rf_vlog(ii)   )./denom;
            
            pFxn(ii) = wTildes2(:,ii)'*normcdf(Z)*wTildes1(:,ii);
            
        case 'rederived'
            denom   = sqrt((thisTSigL(ii).^2).*sum(wTildes2(:,ii).*alphas2(:,ii)).^2 + ...
                           (thisRSigL(ii).^2).*sum(wTildes1(:,ii).*alphas1(:,ii)).^2);
            
            Z       = ( sum(wTildes2(:,ii).*(muTildes2(:,ii) + alphas2(:,ii).*test_vlog(ii))) - ...
                        sum(wTildes1(:,ii).*(muTildes1(:,ii) + alphas1(:,ii).*rf_vlog(ii)  )) )./denom;
            pFxn(ii) = normcdf(Z);
            
    end
end

end

function nll = cNLL_pw(parVec,subjData,vStim1,cStim2,GTgvlog)

%% Likelihood of data given current model pars

% global slAs
slAs = 'nn';

% avlog = parVec(1:numel(vStim1));
% hc    = parVec(numel(vStim1) + 1:end);
avlog = parVec(numel(cStim2) + 1:end);
hc    = parVec(1:numel(cStim2));
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
    switch slAs
        case 'mspf'
            testSlope(x) = avlog(rfVelLin(x) == vStim1);
            testMu(x)    = testVelLog(x) + testSlope(x)*(gvlog(rfVelLin(x) == vStim1)*hc(testCont == cStim2))^2;
            testVar(x)   = (gvlog(rfVelLin(x) == vStim1)*hc(testCont == cStim2))^2;
        case 'nn'
            [~,thisInd]  = min(abs(testVelLin(x) - vStim1)); 
            
            testSlope(x) = avlog(thisInd);
            testMu(x)    = testVelLog(x) + testSlope(x)*(gvlog(thisInd)*hc(testCont == cStim2))^2;
            testVar(x)   = (gvlog(thisInd)*hc(testCont == cStim2))^2;
    end
    
end

% naeker way
Z = (testMu - refMu) ./ sqrt(refVar + testVar);

rho = normcdf(Z);

% Compute log-likelihood - log of bernoulli distribution
nll = -sum(subjData(9,:) .* log(rho) + ...
    (1 - subjData(9,:)) .* log(1 - rho));

end

    