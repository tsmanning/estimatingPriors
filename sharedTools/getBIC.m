function [bic] = getBIC(numPars,numPts,nll)

% Calculate bayesian information criterion for a fit

bic = numPars.*log(numPts) + 2*nll;


end