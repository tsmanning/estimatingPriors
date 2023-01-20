function [nll,ptvs,nllPart] = calculate_model_nll(parVec,ptvs_meas,vStim1,cStim1,vStim2Delta,...
                                   cStim2,numtrials,nllOpts)

%% evaluate psychometric functions predicted by model

% modelType: 'piece', 'gauss', 'mixGauss'
% fitmethod: 'numerical' or 'closed'
% interpOn: 0/1 - interpolate piecewise log-slopes
% nllMethod: 'bernoulli' or 'binomial'
% binScale: 0/1 - Scale binomial NLL to match bernoulli

modelType = nllOpts{1};
fitmethod = nllOpts{2};
interpOn  = nllOpts{3};
nllMethod = nllOpts{4};
binScale  = nllOpts{5};
estDistEval = nllOpts{6};

switch modelType
    
    case 'expon'
        
        gvlog = parVec(1:numel(vStim1));
        hc    = parVec(numel(vStim1) + 1:numel(vStim1) + numel(cStim2));
        avlog = parVec(numel(vStim1) + numel(cStim2) + 1:end);
        
        ptvs = calculate_ptvs(vStim1,cStim1,vStim2Delta,cStim2,avlog,gvlog,hc,fitmethod,interpOn);
        
    case 'gauss'
        
        gvlog = parVec(1:numel(vStim1));
        hc    = parVec(numel(vStim1) + 1:numel(vStim1) + numel(cStim2));
        sigP  = parVec(end);
        
        ptvs = calculate_ptvs_gauss(vStim1,cStim1,vStim2Delta,cStim2,gvlog,hc,sigP,fitmethod);
        
    case 'mixgauss'
        
        numComps = (numel(parVec) - numel(vStim1) - numel(cStim2))/3;
        
        gvlog = parVec(1:numel(vStim1));
        hc    = parVec(numel(vStim1) + 1:numel(vStim1) + numel(cStim2));
        sigP  = parVec(numel(vStim1) + numel(cStim2) + 1:numel(vStim1) + numel(cStim2) + numComps);
        muP   = parVec(numel(vStim1) + numel(cStim2) + numComps + 1:numel(vStim1) + numel(cStim2) + 2*numComps);
        wP    = parVec(numel(vStim1) + numel(cStim2) + 2*numComps + 1:numel(vStim1) + numel(cStim2) + 3*numComps);
        
%         ptvs = calculate_ptvs_mixgauss(vStim1,cStim1,vStim2Delta,cStim2,gvlog,hc,sigP,muP,wP,fitmethod,estDistEval);
        ptvs = calcPFxnMoG(vStim1,cStim1,vStim2Delta,cStim2,gvlog,hc,sigP,muP,wP,fitmethod,estDistEval);
end

% Don't let likelihood go to zero
minlikli = eps;

switch nllMethod
    case 'binomial'

        for n = 1:numel(ptvs)
            
            % evaluate the probability of the data, given the model, and the
            % binominal distribution
            if ~binScale
                nllPart(n) = binopdf(ptvs_meas(n)*numtrials,numtrials,ptvs(n));
            else
                nllPart(n) = binopdf(ptvs_meas(n)*numtrials,numtrials,ptvs(n))./...
                    nchoosek(numtrials,ptvs_meas(n)*numtrials);
            end
        end
        
        nll = -sum(log(max(minlikli,nllPart)));
        
    case 'bernoulli'
        
        % Compute log-likelihood using Bernoulli trial-by-trial method
        for n = 1:numel(ptvs)
            
            choices    = ptvs_meas{n};
            nllPart(n) = -sum( choices.*log(ptvs(n)) + (1-choices).*log(1-ptvs(n)) );
            
        end
        
        nll = sum(max(nllPart,minlikli));

end
