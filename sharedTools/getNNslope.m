function [avlogInds] = getNNslope(vStim1,vStim2Delta)

% Return the indices of the nearest neighbor reference velocity for each
% test velocity used in a 2AFC speed judgment experiment
%
% Usage: [avlogInds] = getNNslope(vStim1,vStim2Delta)
%
% output: nxm matrix of indices, where:
% n = num reference vels, m = num test vels

% ID which slope estimate to use for each vStim2;
numVS1  = numel(vStim1);
numVS2  = numel(vStim2Delta);

vs1Log  = getLogXform(vStim1,0.3);
vs2Log  = getLogXform(vStim1' * vStim2Delta,0.3);

notches = vs1Log(1:end-1) + diff(vs1Log)/2;

avlogInds = nan(numVS1,numVS2);

for i = 1:numVS1
    
    % Find which bin each test velocity falls in (where each bin is defined
    % by the midpoint between reference velocities)
    
    if i == 1
        
        % Handle first sample of prior
        theseLogicInds = vs2Log <= notches(1);
        
    elseif i == numVS1
        
        % Handle last sample of prior
        theseLogicInds = vs2Log > notches(end);
        
    else
        
        % Handle other cases
        theseLogicInds = (vs2Log > notches(i-1)) & ...
                         (vs2Log < notches(i));
        
    end
    
    avlogInds(theseLogicInds) = i;
    
end

end