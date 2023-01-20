function [vTilde] = getLogXformSigned(v,v0)

% Get log transformed speed
%
% Usage: [vTilde] = getLogXform(v,v0)

vTilde = sign(v).*log(1 + abs(v)/v0);

end