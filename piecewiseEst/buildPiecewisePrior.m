function [prior,vels] = buildPiecewisePrior(a)

% Reconstruct estimated prior from a set of variables
%%% also returns support vels if you give it an empty set of slopes

%% Define velocity support
v        = [0.5 1 2 4 8 12];
vLog     = getLogXform(v,0.3);

% Set midpoints between reference velocities as bin edges
vLogMids = 0.5*(vLog(1:end-1) + vLog(2:end));

% Difference between reference velocities and bin edges
vDiff    = vLogMids - vLog(1:end-1);

% Lower bin edges
vStarts  = [(vLog(1) - vDiff(1)) vLogMids];
% vStarts  = [0 vLogMids];

% Upper bin edges
vEnds    = [vLogMids (vLog(end) + vDiff(end))];

% Bin sizes
dx       = vEnds - vStarts;

%% Estimate prior as patchwork of connected exponents

numPts = numel(v) + 1;

% Numerically integrate to recover non-normalized prior
if ~isempty(a)
    for ii = 1:numPts
        
        if ii == 1
            p(ii) = 1;
        else
            p(ii) = p(ii-1) + a(ii-1)*dx(ii-1);
        end
        
    end
    
    % Exponentiate and normalize prior
    pExp  = exp(p);
    
    AUC   = sum( dx.*(0.5*(pExp(2:end) + pExp(1:end-1))) );
    prior = pExp/AUC;
else
    prior = [];
end

vels  = [vStarts vEnds(end)];

%% debug
if 0
gauF = @(x,mu,sig) (1/(sig*sqrt(2*pi)))*exp(-0.5*((x-mu)/sig).^2);    
% vels2 = [0 vels];
% gt = gauF(vels2,0,1);
% gt = gt/sum(0.5*[vels2(2)-vels2(1) dx].*(gt(2:end) + gt(1:end-1)));

gt = gauF(vels,0,1);
gt = gt/sum(0.5*dx.*(gt(2:end) + gt(1:end-1)));

for ii = 1:numel(v)
    vLab{ii} = num2str(v(ii));
end

figure;
hold on;
plot(vels,prior,'r');
% plot(vels2,gt,'--k');
plot(vels,gt,'--k');
set(gca,'xtick',getLogXform(v,0.3),'xticklabel',vLab,'yscale','lin');
end

end