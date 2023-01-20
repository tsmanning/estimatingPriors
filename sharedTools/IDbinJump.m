function [simData] = IDbinJump(simData,binDef,split)

% ID MAP estimates that jump outside the expected bin under the S&S
% piecewise model
%
% Usage: [] = IDbinJump(simData)

% Find midpoint on log scale between reference velocities

s1Vels      = unique(arrayfun(@(x) x.s1Vel,simData));
s1VelsLog   = getLogXform(s1Vels,0.3);
numRefVels  = numel(s1Vels);
numTestVels = numel(simData(1).s2Vels);

MPs         = s1VelsLog(1:end-1) + diff(s1VelsLog)/2;

binEdges    = MPs;

% Get unique SD values

SDs         = unique(arrayfun(@(x) x.s2SD,simData));
numSDs      = numel(SDs);

numTr       = numel(simData(1).likeMeans1);

%% ID trials where the likelihood was jittered out of the bin

for ii = 1:numel(simData)
    
    % Get indices of likelihood EVs and SDs
    velInd  = find(simData(ii).s1Vel == s1Vels);
%     sdInd   = find(simData(ii).s2SD == SDs);
    
    % Get slopes for this split exponent
    %%% this would need to be fixed for non-split expon case
    
    slopes  = simData(ii).priorPar{1};
    
    if simData(ii).s1Vel < split
        slope   = slopes(1);
    else
        slope   = slopes(2);
    end
    
    % Get amount of bias for each posterior
    biasRef     = slope*(simData(ii).s1SD)^2;
    biasTest    = slope*(simData(ii).s2SD)^2;
    
    %% Check for means crossing over into other bins unexpectedly
    
    % Reference Stims
    if velInd == 1
        
        % First bin
        oobInds1 = simData(ii).likeMeans1 > binEdges(1);
        
    elseif velInd == numRefVels
        
        % Last bin
        oobInds1 = simData(ii).likeMeans1 < binEdges(end);
        
    else
        
        % Middle bins
        oobInds1 = (simData(ii).likeMeans1 < binEdges(velInd-1)) | ...
                   (simData(ii).likeMeans1 > binEdges(velInd));
        
    end
    
    % Test Stims (which are expected to fall in certain bins, though not
    % always the reference bin)
    s2NoiseBins = nan(numTr,numTestVels);
    s2Bins      = nan(1,numTestVels);
    tvels       = getLogXform(simData(ii).s2Vel,0.3);
    
    for jj = 1:numRefVels
       
        % Recover bin indices for noisy and noiseless test stimuli
        
        if jj == 1
            
            % Handle first sample of prior
            theseLogicInds = simData(ii).likeMeans2 <= binEdges(1);
            tli2           = tvels <= binEdges(1);
            
        elseif jj == numRefVels
            
            % Handle last sample of prior
            theseLogicInds = simData(ii).likeMeans2 > binEdges(end);
            tli2           = tvels > binEdges(end);
            
        else
            
            % Handle other cases
            theseLogicInds = (simData(ii).likeMeans2 > binEdges(jj-1)) & ...
                             (simData(ii).likeMeans2 < binEdges(jj));
            tli2           = (tvels > binEdges(jj-1)) & ...
                             (tvels < binEdges(jj));
            
        end
        
        s2NoiseBins(theseLogicInds) = jj;
        s2Bins(tli2) = jj;
        
    end
    
    % ID trials for which bin indices aren't equal
    s2bExt = repmat(s2Bins,[numTr,1]);
    oobInds2 = s2bExt ~= s2NoiseBins;
    
    %% Check if biases would push posterior beyond lower bound of nearest
    % neighbor bin (only occurs when slope is decreasing as a function of
    % speed)
    
    % Reference
    if velInd > 1
        
        edgeDistLB1  = simData(ii).likeMeans1-binEdges(velInd-1);    
%         overLimit1 = ( (abs(slopes(2)) > abs(slopes(1)) ) & ...
%                        (edgeDistLB1 + biasRef <= 0) ) | ...
%                      oobInds1;
        overLimit1 = (edgeDistLB1 + biasRef <= 0) | ...
                     oobInds1;
                 
    else
        
        overLimit1  = oobInds1;     
        
    end
    
    % Test
    %%%%%% need to index into nearest neighbor bins here
    tBin         = [-5 binEdges];
    edgeDistLB2  = simData(ii).likeMeans2 - tBin(s2bExt);
    
%     overLimit2   = ( (abs(slopes(2)) > abs(slopes(1)) ) & ...
%                      (edgeDistLB2 + biasTest <= 0) ) | ...
%                      oobInds2;
    overLimit2   =  (edgeDistLB2 + biasTest <= 0) | ...
                     oobInds2;
           
    % Mark trials to include (i.e. eliminate anything testing positive
    % above)
    switch binDef
        
        case 'midpoint'
            % Use midpoint between reference velocities as the bin edges
           
            simData(ii).incInds1 = ~oobInds1;
            simData(ii).incInds2 = ~oobInds2;
            
        case 'minDist'
            % Find alternative cutoff where large likelihood widths result in
            % alternative MAP - just use the ground truth slope as the cutoff?
            
            simData(ii).incInds1 = ~overLimit1;
            simData(ii).incInds2 = ~overLimit2;
            
    end
    
%     if ii == 22
%         keyboard
%     end
    
end

% keyboard;

end