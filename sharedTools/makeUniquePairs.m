function [thesePairs] = makeUniquePairs(m)

% For a set of elements of length m, return a 2xn matrix of each unique
% pairing of elements

thesePairs = [];

for i = 1:m-1

    thesePairs = [thesePairs;[i*ones(m-i,1) [i+1:m]']];

end

end