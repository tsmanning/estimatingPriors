function [outMat] = makeCombos(levels)

% Basically copy fullfact for installs without stats toolbox

if isrow(levels)
    levels = levels';
end

if size(levels,2) ~= 1
    error('Input argument must be a vector.');
end

numLevels = numel(levels);
ind = 1;

while ind <= numLevels
    
    if ind == 1
        numCombs = levels(ind);
    else
        numCombs = numCombs * levels(ind);
    end
    
    ind = ind + 1;
    
end

outMat = nan(numCombs,numLevels);

for ind = 1:numLevels
    
    levelVec = [1:levels(ind)]';
    
    % Number of element-wise repeats
    elemMulti = prod(levels(1:ind-1));
    
    % Number of vector-wise repeats
    matMulti  = prod(levels(ind+1:end));
    
    outMat(:,ind) = repmat(repelem(levelVec,elemMulti,1),[matMulti 1]);

end

end