function [coinFlip,cumGauss] = getModelFitBounds(datamat)

% Get NLL bounds to contextualize Bayesian ideal observer fits


% Get NLL for coinflip
coinFlip = -sum( log(0.5) * size(datamat,2) );

% Get NLL for cumulative Gaussian
v1s = unique(datamat(1,:));
c1s = unique(datamat(2,:));
c2s = unique(datamat(4,:));

numV1 = numel(v1s);
numC1 = numel(c1s);
numC2 = numel(c2s);

NLL = 0;

for ii = 1:numV1
    
    for kk = 1:numC2
        
        for jj = 1:numC1
            
            % Find all columns for this pFxn
            theseInds = datamat(1,:) == v1s(ii) & ...
                        datamat(2,:) == c1s(jj) & ...
                        datamat(4,:) == c2s(kk);
            
            % Grab the set of responses
            theseResps = datamat(9,theseInds);
            
            if ~isempty(theseResps)
                
                % Grab test stimulus velocities and indices
                [theseVels,~,tv] = unique(datamat(3,theseInds));
                
                pResp = nan(1,numel(theseVels));
                pRespTr = nan(1,numel(theseVels));
                
                % Calculate proportion of responses
                for hh = 1:numel(theseVels)
                    vInds = tv == hh;
                    pResp(hh)   = sum(theseResps(vInds))/sum(vInds);
                    pRespTr(hh) = sum(vInds);
                end
                
                % Fit data with cumulative Gaussian
                thisDataFit = fitCumGauss(theseVels,pResp,pRespTr);
                
                % Get nll for this set
                NLL = -sum( theseResps.*log(thisDataFit(tv)) + ...
                           (1-theseResps).*log(max(1-thisDataFit(tv),eps)) ) + NLL;

            end
            
        end
        
        
    end
    
end

cumGauss = NLL;

end

