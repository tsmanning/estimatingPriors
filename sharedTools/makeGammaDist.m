function [gd_pdf] = makeGammaDist(x,k,th)

% Make a discretely sampled gamma distribution given a support and two
% parameters

gd_pdf = (1/gamma(k))*(1/(th^k))*x.^(k-1).*exp(-x/th);

end