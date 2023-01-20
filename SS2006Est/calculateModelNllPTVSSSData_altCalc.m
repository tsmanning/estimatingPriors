function nll = calculateModelNllPTVSSSData_altCalc(parVec,rawDat,vStim1)

%% Likelihood of data given current model pars

a       = parVec(1:numel(vStim1));
gvlog   = parVec(numel(vStim1) + 1:2*numel(vStim1));
sigVals = parVec(2*numel(vStim1) + 1:end);

% Pair/interpolate slope values to each stimulus velocity
% These two are 1x5760 (num trials in original exp)
refV = rawDat(1,:);
testV = rawDat(3,:);
refSlope  = interp1(vStim1,a,refV, 'linear');
testSlope = interp1(vStim1,a,testV, 'linear', 'extrap');

% Index contrast-dependent likelihood fits
[~,~,CindVec] = unique([rawDat(2,:) rawDat(4,:)]);
[~,~,VindVec] = unique(rawDat(1,:));
refSig  = sigVals(CindVec(1:5760)).*gvlog(VindVec);
testSig = sigVals(CindVec(5761:end)).*gvlog(VindVec);

% Get standard scores for each difference in v_hat between stims on each trial
Z = (testV - refV + refSlope .* refSig.^2 - testSlope .* testSig.^2) ./ ...
    sqrt(refSig .^ 2 + testSig .^ 2);

% Get psychometric function by feeding standard scores into normCDF
EPS = 1e-12;
rho = min(max(normcdf(Z), EPS), 1 - EPS);

% Compute log-likelihood of each trial response using Bernoulli
% distribution and sum over all trials
choices = rawDat(9,:);
nll = -sum(choices .* log(rho) + (1 - choices) .* log(1 - rho));


end