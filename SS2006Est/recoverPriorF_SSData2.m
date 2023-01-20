function [pFVec,vVec] = recoverPriorF_SSData2(avlogF,hcF,gvlogF,nllF)
% function [pFVec,vLog] = recoverPriorF_SSData2(avlogF,hcF,gvlogF,nllF)
% For a set of input slopes, recover estimate of prior shape

%% Load sim data and fits

% % Grab Gaussian fits
% load([saveDir,filesep,fName,'_gs'],'nllF','sigPF','gvlogF','hcF');
% [~,nllFgs] = min(nllF);
% sigPFgs    = sigPF(nllFgs);
% g_gs       = gvlogF(nllFgs,:);
% hc_gs      = hcF(nllFgs,:);

% Plot sample with best likelihood
[~,bInd] = min(nllF);

a       = avlogF(bInd,:);
hc      = hcF(bInd,:);
g       = gvlogF(bInd,:);

%% Setup pars

% Likelihoods of sigmas (make it not-speed dependent for now)
likeSigs = abs(g'*hc);

% Speeds & such
v       = [0.5 1 2 4 8 12];
vLog    = getLogXform(v,0.3);

% c       = [0.05 0.075 0.1 0.2 0.4 0.5 0.8];
% for i = 1:numel(c)
%     cLab{i} = num2str(c(i));
% end

vLogMids    = 0.5*(vLog(1:end-1) + vLog(2:end));

vDiffsUp    = -vLog(1:end-1) + vLogMids;
vDiffsDown  =  vLog(2:end)   - vLogMids;

vStarts  = [(vLog(1) - vDiffsDown(1)) vLogMids];
vEnds    = [vLogMids (vLog(end) + vDiffsUp(end))];

% vLogFull = linspace(vStarts(1),vEnds(end),50);
% vLogInit = linspace(0,vStarts(1),10);
% vSpacing = diff(vLogFull(1:2));

vVals = [vStarts' vLog' vEnds'];

vVec = unique(vVals);

% Anonymous functions
expF = @(x,a,b) exp(a*x + b);
gauF = @(x,mu,sig) (1/(sig*sqrt(2*pi)))*exp(-0.5*((x-mu)/sig).^2);

% Solve for next exponential intercept to patch piecewise bits together
b2 = @(a1,a2,b1,x) (a1-a2)*x + b1; 

%% Estimate prior as patchwork of connected exps

patchFxn = nan(numel(v),3);
vLab = cell(numel(v),1);

for i = 1:numel(v)
    
    if i == 1
        b(i)      = 0;
    else
        b(i)      = b2(a(i-1),a(i),b(i-1),vStarts(i));
    end
    
    thisb         = b(i);
        
    patchFxn(i,:) = expF(vVals(i,:),a(i),thisb);
    
    vLab{i}       = num2str(v(i));
    
end

rarVec = [1 7 2 8 3 9 4 10 5 11 6 12 18];

pFVec = patchFxn(rarVec)';

% Numerically integrate to normalize AUC to 1 (trapezoid method)
AUC = sum( (vVec(2:end) - vVec(1:end-1)).*(0.5*(pFVec(2:end) + pFVec(1:end-1))) );

pFVec = pFVec/AUC;


%% Plot prior

if 0
% titleC = regexp(fName,'_','split');
% temp = titleC{2};

% contInd = contL;

% Semilog
f1 = figure;
f1.Position = [100 100 625 525];
hold on;
plot(vVec,pFVec,'--r','linewidth',3);
scatter(vLog,pFVec(2:2:12),100,'r','filled');
set(gca,'xtick',vLog,'xticklabel',vLab,'yscale','lin','ylim',[0 5],...
    'xlim',[vStarts(1) vEnds(end)],'fontsize',20);
ylabel('p(vel)');
xlabel('Stimulus Velocity');

for i = 1:7
    
%     thisSig = likeSigs(1,i);
    
    xAx = linspace(vLog(1),vLog(end),50);
    pF  = gauF(xAx,vLog(3),likeSigs(1,i));
    plot(xAx,pF,'b','linewidth',2);
    
end

legend({'Piecewise Fit','Reference Velocities','Likelihoods'});

% Log
f2 = figure;
f2.Position = [745 100 625 525];
hold on;
plot(vVec,pFVec,'--r','linewidth',3);
scatter(vLog,pFVec(2:2:12),100,'r','filled');
set(gca,'xtick',vLog,'xticklabel',vLab,'yscale','log','ylim',[0 5],...
    'xlim',[vStarts(1) vEnds(end)],'fontsize',20);
ylabel('p(vel)');
ylabel('p(vel)');
xlabel('Stimulus Velocity');
% title([titleC{1},' (',titleC{2},')']);

for i = 1:7
    
%     thisSig = likeSigs(1,i);
    
    xAx = linspace(vLog(1),vLog(end),50);
    pF  = gauF(xAx,vLog(3),likeSigs(1,i));
    plot(xAx,pF,'b','linewidth',2);
    
end

legend({'Piecewise Fit','Reference Velocities','Likelihoods'},'location','southwest');
end

end