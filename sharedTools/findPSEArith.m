function [PSE,bias,biasPer] = findPSEArith(vStim1,cStim1,cStim2,method,pars)

% Arithmetically estimate PSE when close form solution exists
%
% Usage: [PSE,bias,biasPer] = findPSEArith(vStim1,cStim1,cStim2,method,pars)

% Make a matrix of biases/PSEs for each combination of [vStim1,cStim1,cStim2]
% e.g. 6x7x2 in S&S 2006 

numV1 = numel(vStim1);
numC1 = numel(cStim1);
numC2 = numel(cStim2);

if ~iscolumn(vStim1)
    vStim1 = vStim1';
end


% Get inds for cStim1 (assuming cStim1 is a subset of cStim2)
for i = 1:numC1
    C1Ind(i) = find(cStim1(i) == cStim2);
end


% Grab likelihood parameters
gvlog   = pars.gvlog;
hc      = pars.hc;

% Reshape pars into desired output matrix size
vRef    = repmat(vStim1,[1 numC2 numC1]);

varRef  = repmat(reshape(gvlog'*hc(C1Ind),[numV1 1 numC1]),[1 7 1]);
varTest = repmat(gvlog'*hc,[1,1,2]);


% Calculate biases
switch method
    
    case 'piecewise'
        
        avlog = pars.avlog;
        
        %%% setting aTest to avlog for now, but may not capture best fitting PSE
        %%% like the numeric method does for the nearest neighbor approach. that
        %%% method can also produce multiple PSEs
        aRef    = repmat(avlog',[1 numC2 numC1]);
        aTest   = repmat(avlog',[1 numC2 numC1]);
        
        bias    = aRef.*varRef - aTest.*varTest;
        PSE     = vRef + bias;
        
        biasPer = PSE./vRef - 1;

    case 'gaussian'
        
        sigP    = pars.sigP;
        
        % PSE is a fixed ratio based on variances of likelihoods and prior
        
        biasPer = (sigP^2 + varTest)./(sigP^2 + varRef);
        
        PSE     = biasPer*vRef;
        
        bias    = PSE - vRef;
        
    case 'mixGauss'
        
        error('No current analytic solution for PSE with mixture of Gaussians method.');

end