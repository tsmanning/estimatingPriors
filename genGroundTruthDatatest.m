function [simData] = genGroundTruthDatatest(priorF,priorPar,like1Par,like2Par,...
                                domain,numTr,lapseRate)


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

x0 = 0.3;

muL1 = getLogXformSigned(muL1,x0);
muL2 = getLogXformSigned(muL2,x0);


%% Define functions and axes, shuffle RNG seed

numS2Means = numel(muL2);

% Make axes along domain in which observer encodes stims:
% smaller mu - 3xlargest SD : larger mu + 3x largest SD 
x = linspace(min([muL1 muL2]) - 3*max([sigL1 sigL2]), ...
             max([muL1 muL2]) + 3*max([sigL1 sigL2]),600);
         
numSamps = numel(x);

% Assert background number of reps before subsampling
numReps = numTr;

% Can be used for multiple means if x is a row vector and mu is a column vector
gauF = @(x,mu,sig) (1./(sig*sqrt(2*pi))).*exp(-0.5*((x-mu)./sig).^2);

% Shuffle random number generator
rng('shuffle');
        
% For a prior composed of multiple Gaussian components
parMat  = priorPar{1};

gamma   = parMat(1,:);
mu      = getLogXformSigned(parMat(2,:),0.3);
w       = parMat(3,:);

priorVals = zeros(1,numel(x));

for ii = 1:numel(w)
    priorVals = priorVals + w(ii)*normpdf(x,mu(ii),gamma(ii));
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

% Calculate posteriors
postFxns1       = repmat(priorVals,[numReps,1]).*lFxns1;

% ... And again for the set of stim 2
postFxns2       = repmat(priorVals,[numReps,1,numS2Means]).*lFxns2;

MAPs1           = sum(postFxns1.*repmat(x(1:numSamps),[numReps 1]),2)./sum(postFxns1,2);

MAPs2           = sum(postFxns2.*repmat(x(1:numSamps),[numReps 1 numS2Means]),2)./sum(postFxns2,2);
MAPs2           = squeeze(MAPs2);


%% Simulate 2AFC experiment using binomial distribution

% Generate long binary list of which stimulus was chosen
s2Choice        = MAPs2 > MAPs1;

% Make estimate of "ground truth" psychometric curve
pSDTFull        = sum(s2Choice)./numReps;


%% Collect into output structure
simData.pSDT            = pSDTFull; 


end