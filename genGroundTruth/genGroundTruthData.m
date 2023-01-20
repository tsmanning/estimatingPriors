function [simData] = genGroundTruthData(priorF,priorPar,like1Par,like2Par,...
                                domain,numTr,lapseRate)

% Generates a 2AFC data set based on Bayesian estimation without making any
% underlying assumptions
%
% For ONE unique triplet of [C_ref,V_ref,C_test] and several different
% V_test; this function looped over in autoGenGroundTruth to generate full
% set
%
% Usage: [pSDT,simProps,numTrials] = genGroundTruthData(priorF,priorPar,like1Par,like2Par,numTr)
%   - Likelihood is assumed to be Gaussian with muL and sigL
%   - priorF is either 'expon' or 'gauss' according to PDF desired
%   - priorPar should be 2x1 or 1x2 vector of parameters (mu and sigma or
%     multiplier and offset of exponential)
%   - like1Par should be 2x1 or 1x2 vector of mu and sigma
%   - like2Par should be a 2x1 or 1x2 cell array of vectors containing
%     desired means and sigmas
%   - numTr determines number of unique trials run for each Stim 2
%
% Example call: genGroundTruthData('gauss',[0 1],[5 1],{0:10,1},'lin',135,1);
%               genGroundTruthData('expon',[-7 0],[2 0.25],{2:0.5:8,0.38},'log',80,1);

%% Make sure inputs are in correct dimensions

muL1 = like1Par(1);
sigL1 = like1Par(2);

muL2 = like2Par{1};
sigL2 = like2Par{2};

% mu and sigma should both be row vectors
if (numel(sigL2) ~= 1) && (numel(muL2) ~= numel(sigL2))
   error('Means and sigmas of 2nd stimulus likelihoods must be same size.') 
end

if ~isrow(muL2)
    muL2 = muL2';
end

if ~isrow(sigL2)
    sigL2 = sigL2';
end

% If only one sigma for second likelihood, assume sigma is the same across
% all likelihoods
if numel(sigL2) == 1
    [s1,s2] = size(muL2);
    sigL2 = sigL2*ones(s1,s2);
end

% If we're in the log domain (Stocker & Simoncelli), transform the means,
% but not the sigmas (these are already assumed to be in log space)

if strcmp(domain,'log')
    
    % Log Normalization value if using (maybe don't hardcode in future
    % iterations)
    x0 = 0.3;
    
    muL1 = getLogXformSigned(muL1,x0);
    muL2 = getLogXformSigned(muL2,x0);
    
end


%% Define functions and axes, shuffle RNG seed

numS2Means = numel(muL2);

% Make axes along domain in which observer encodes stims:
% smaller mu - 3xlargest SD : larger mu + 3x largest SD 
x = linspace(min([0 muL1 muL2]) - 3*max([sigL1 sigL2]), ...
             max([0 muL1 muL2]) + 3*max([sigL1 sigL2]),800);
         
% x = linspace(-5,5,800);
         
numSamps = numel(x);

% Assert background number of reps before subsampling
%%%%% this should just be numTr, since subsampling is done outside this
%%%%% script
numReps = numTr;

% Can be used for multiple means if x is a row vector and mu is a column vector
gauF = @(x,mu,sig) (1./(sig*sqrt(2*pi))).*exp(-0.5*((x-mu)./sig).^2);

expF = @(x,a,b) exp(a*x + b);

% Shuffle random number generator
rng('shuffle');


%% Define prior distribution

% Do it numerically for now
%%%%% NOTE: THESE PRIORS AREN'T NORMALIZED TO 1
switch priorF
    
    case 'gauss'
        
        muP  = priorPar(1);
        sigP = priorPar(2);
                
        % Gaussian Prior
        priorVals = gauF(x,muP,sigP);
        
    case 'expon'
        
        aP = priorPar(1);
        bP = priorPar(2);        
        
        % Exponential
        priorVals = expF(x,aP,bP);
                
    case 'splitExp'
        
        % For a prior composed of two exponentials with different powers
        aP        = priorPar{1};
        bP        = priorPar{2}; 
        splitLoc  = getLogXform(priorPar{3},0.3);
        
        % right now setup to take a single split, and single b val.
        priorVals = makePiecewisePrior(x,aP(1),aP(2),bP,splitLoc);
        
    case 'mixgauss'
        
        % For a prior composed of multiple Gaussian components
        parMat  = priorPar{1};
        
        gam   = parMat(1,:);
        mu      = getLogXformSigned(parMat(2,:),0.3);
        w       = parMat(3,:);
        
        priorVals = zeros(1,numel(x));
        
        for ii = 1:numel(w)
            priorVals = priorVals + w(ii)*gauF(x,mu(ii),gam(ii));
        end
        
    case 'genGauss'
        
        % For a generalized Gaussian prior
        genNorm = @(x,u,a,b) b/(2*a)*(1/gamma(1/b))*exp(-(abs(x-u)/a).^b);
        
        uP = priorPar(1);
        aP = priorPar(2);
        bP = priorPar(3);
        
        priorVals = genNorm(x,uP,aP,bP);
        
end


%% Get closed form solutions for psychometric curves (or at least our high n numeric approx)
switch priorF
    
    case 'gauss'
        
        % mean and SIGMA of posterior in tranformed velocity space
        [muP1,~] = findGaussProd(muL1,sigL1,priorPar(1),priorPar(2));
        alpha1   = (priorPar(2)^2)/(sigL1^2 + priorPar(2)^2);
        sigP1    = (alpha1)*(sigL1);
        
        [muP2,~] = findGaussProd(muL2,sigL2,priorPar(1),priorPar(2));
        alpha2   = (priorPar(2)^2)./(sigL2.^2 + priorPar(2)^2);
        sigP2    = (alpha2).*(sigL2);
        
        % loop through P2 values
        for m = 1:numel(muP2)
            
            Z = (muP2(m) - muP1) ./ ...
                 sqrt(sigP2(m)^2 + sigP1^2);
            closedForm_pSDT(m) = normcdf(Z);
            
        end
        
    case 'expon'
        
        % mean and var of posterior in tranformed velocity space        
        muP1   = muL1 + priorPar(1)*(sigL1).^2;
        sigP1  = sigL1;
        
        muP2   = muL2 + priorPar(1)*(sigL2).^2;
        sigP2  = sigL2;
        
        % loop through P2 values
        for m = 1:numel(muP2)
            
            Z = (muP2(m) - muP1) ./ ...
                 sqrt(sigP2(m)^2 + sigP1^2);
             
            closedForm_pSDT(m) = normcdf(Z);
        end
        
    case 'splitExp'
        
        %%%% just assume a 2 slope prior for now    
        
        % Set mean of reference MAP estimate
        if muL1 < splitLoc
            muP1   = muL1 + aP(1)*(sigL1).^2;  
        else
            muP1   = muL1 + aP(2)*(sigL1).^2;  
        end
        
        % Set means of test MAPs estimate
        sl1 = sum(muL2 < splitLoc);
        sl2 = sum(muL2 >= splitLoc);
        
        priorVec = [aP(1)*ones(1,sl1) aP(2)*ones(1,sl2)];
        
        muP2   = muL2 + priorVec.*(sigL2).^2;
        
        sigP1  = sigL1;
        sigP2  = sigL2;
        
        % loop through P2 values
%         for m = 1:numel(muP2)
%             
%             Z = (muP2(m) - muP1) ./ ...
%                  sqrt(sigP2(m)^2 + sigP1^2);
%              
%             closedForm_pSDT(m) = normcdf(Z);
%         end

        %%% WHY NOT LOOP THROUGH THIS TOO??
        Z = (muP2 - muP1) ./ sqrt(sigP2.^2 + sigP1^2);
        
        closedForm_pSDT = normcdf(Z);
        
    case 'mixgauss'
        
        % need to work out
        
        closedForm_pSDT = nan;
        
end


%% Define likelihood distribution

% Do it numerically for now
%
% n = num Trials, m = num Samples, p = num s2 means
%
% Generate set of noisy estimates for each stimulus (nx1 or nxp for s2).
% Using repmat for both mu and sigma if sigma is dependent on mu (e.g.
% weber's law in linear domain)
likeMeans1          = muL1 + sigL1*randn([numReps,1]);
likeMeans2          = repmat(muL2,[numReps,1]) + ...
                      repmat(sigL2,[numReps,1]).*randn([numReps,numS2Means]);

% Define function in lin or log
likeF               = @(x,mu,sigma) gauF(x,mu,sigma);

% Generate set of likelihood functions 
% ( n x m x p: m = numel(x), p = numel(muL2) )
lFxns1              = likeF(x,likeMeans1,sigL1);

% Could probably get this out of a loop with some adjustment to the
% anonymous function
lFxns2              = nan(numReps,numSamps,numS2Means);
for i = 1:numS2Means
    lFxns2(:,:,i)   = likeF(x,likeMeans2(:,i),sigL2(:,i));
end


%% Calculate set of posteriors, MAP estimates, and PDF

% Do it numerically for now, set it run analytically if using gaussians
% and exponentials

% Calculate posteriors
postFxns1       = repmat(priorVals,[numReps,1]).*lFxns1;
% Get index of MAP
[~,MAPinds1]    = max(postFxns1,[],2);
% Get MAP
MAPs1           = x(MAPinds1)';

% ... And again for the set of stim 2
postFxns2       = repmat(priorVals,[numReps,1,numS2Means]).*lFxns2;
[~,MAPinds2]    = max(postFxns2,[],2);
MAPs2           = squeeze(x(MAPinds2));

%%%%%%%%%%%%%%%%% calculate if Bayes estimator is mean and not MAP
ESTs1           = sum(postFxns1.*repmat(x(1:numSamps),[numReps 1]),2)./sum(postFxns1,2);

ESTs2           = sum(postFxns2.*repmat(x(1:numSamps),[numReps 1 numS2Means]),2)./sum(postFxns2,2);
ESTs2           = squeeze(ESTs2);


%% Simulate 2AFC experiment using binomial distribution

% Generate long binary list of which stimulus was chosen
% s2Choice        = MAPs2 > MAPs1;
s2Choice        = ESTs2 > ESTs1;
numEntries      = numel(s2Choice);

% Make estimate of "ground truth" psychometric curve
pSDTFull        = sum(s2Choice)./numReps;

% Random sampling of linear indices for "lapse trials" without replacement
lapseInds       = datasample(1:numEntries,round(lapseRate*numEntries),'Replace',false);

% Assign s2 > s1 based on coin flip for these trials
choices                 = rand(numel(lapseInds),1) > 0.5;
s2ChoiceL               = s2Choice;
s2ChoiceL(lapseInds)    = choices;

% Subsample according to desired number of trials to simulate (something
% roughly human - numTr*numS2Means: 80*11 = 880 trials; 135*11 = 1485 trials)
% simChoices = nan(numTr,numS2Means);
%%%%%%%%%%% should just eliminate this, since subsampling is done outside
%%%%%%%%%%% this script

for i = 1:numS2Means
    randTrInds      = datasample(1:numReps,numTr,'Replace',false);
    simChoices(:,i) = s2ChoiceL(randTrInds,i);
end

% Calculate proportion of S2 choices in that sub sample
simProps = sum(simChoices)./numTr;
numTrials = numTr*ones(numS2Means,1);


%% Collect into output structure
simData.cf_pSDT         = closedForm_pSDT;  % Calculated directly with gs/ss
simData.pSDT            = pSDTFull;         % No lapse trials
simData.simChoices      = simChoices;       % Matrix of all binary choices
simData.simProps        = simProps;         % pSDT with lapse rate
simData.numTrials       = numTrials;        % Number of trials simulated
simData.s1Vels          = muL1;
simData.s2Vels          = muL2;

% Ground truth information
simData.domain      = domain;
simData.priorF      = priorF;
simData.priorPar    = priorPar;
simData.s1Vel       = like1Par(1);
simData.s1SD        = like1Par(2);
simData.s2Vel       = like2Par{1};
simData.s2SD        = like2Par{2};
simData.ESTs1       = ESTs1;
simData.ESTs2       = ESTs2;
simData.MAPs1       = MAPs1;
simData.MAPs2       = MAPs2;

% Save likelihood means (with additive noise) for debugging
simData.likeMeans1  = likeMeans1;
simData.likeMeans2  = likeMeans2;

end