function [alphas,muTildes,wTildes,xHat] = getMogPostPars(prior_mus,prior_sigs,prior_ws,lik_mu,lik_sig)

% For a given mixture of Gaussians prior (prior_prior_mus,prior_sigs,prior_ws) and likelihood (lik_mu,lik_sig), 
% returns parameters for posterior and BLS estimate

% Output mat dim 1: component n of MoG 
% Output mat dim 2: stimulus m

% Make sure prior pars are column vecs, likelihood pars are rows
if ~iscolumn(prior_ws)
    prior_ws = prior_ws';
end
if ~iscolumn(prior_sigs)
    prior_sigs = prior_sigs';
end
if ~iscolumn(prior_mus)
    prior_mus = prior_mus';
end

if ~isrow(lik_sig)
    lik_sig = lik_sig';
end
if ~isrow(lik_mu)
    lik_mu = lik_mu';
end

% Standard Normal
phi        = @(x) (1/sqrt(2*pi))*exp(-0.5*(x.^2));

% Mixture component shrinkage factors
alphas     = (prior_sigs.^2) ./ (prior_sigs.^2 + lik_sig.^2);

% Mixture component means
muTildes   = ( (lik_sig.^2) ./ (lik_sig.^2 + prior_sigs.^2) ).*prior_mus;

% Adjusted mixture component weights (note this is x-dependent!)
wTildes    = (prior_ws./sqrt(prior_sigs.^2 + lik_sig.^2)).* ...
             phi((lik_mu-prior_mus) ./ sqrt(prior_sigs.^2 + lik_sig.^2));
wTildes    = wTildes./sum(wTildes,1);

% Get BLS estimate
xHat       = sum(wTildes.*(alphas.*lik_mu + muTildes));

end