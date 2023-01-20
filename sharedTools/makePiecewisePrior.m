function [priorF] = makePiecewisePrior(x,a1,a2,b,split)

% Make a prior compose of multiple components (assume everything's
% exponential for now, can come back and add Gaussian, etc later on)

%% Split dependent variable if needed

a1Inds = (x <= split);


%% Get prior distribution
% For simple two-piece prior, can just get away with negating logical
% indexing

expF = @(x,a,b) exp(a*x + b);

b1 = b;
b2 = (a1-a2)*split + b;

priorF = [expF(x(a1Inds),a1,b1) expF(x(~a1Inds),a2,b2)];



end