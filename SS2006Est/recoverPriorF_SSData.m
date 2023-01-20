function [pFVec,vLog] = recoverPriorF_SSData(saveDir,fName,contL)

% For a set of input slopes, recover estimate of prior shape

%% Load sim data and fits

% Grab Gaussian fits
load([saveDir,filesep,fName,'_gs'],'nllF','sigPF','gvlogF','hcF');
[~,nllFgs] = min(nllF);
sigPFgs    = sigPF(nllFgs);
g_gs       = gvlogF(nllFgs,:);
hc_gs      = hcF(nllFgs,:);

% Grab piecewise fits
load([saveDir,filesep,fName]);

% Plot sample with best likelihood
[~,bInd] = min(nllF);

a       = -avlogF(bInd,:);
hc      = hcF(bInd,:);
g       = gvlogF(bInd,:);

%% Setup pars

% Likelihoods of sigmas (make it not-speed dependent for now)
likeSigs = abs(g'*hc);

% Speeds & such
v       = [0.5 1 2 4 8 12];
vLog    = getLogXform(v,0.3);

c       = [0.05 0.075 0.1 0.2 0.4 0.5 0.8];
for i = 1:numel(c)
    cLab{i} = num2str(c(i));
end

vLogMids    = 0.5*(vLog(1:end-1) + vLog(2:end));

vDiffsUp    = -vLog(1:end-1) + vLogMids;
vDiffsDown  =  vLog(2:end)   - vLogMids;

vStarts  = [(vLog(1) - vDiffsDown(1)) vLogMids];
vEnds    = [vLogMids (vLog(end) + vDiffsUp(end))];

vLogFull = linspace(vStarts(1),vEnds(end),50);
vLogInit = linspace(0,vStarts(1),10);
vSpacing = diff(vLogFull(1:2));

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

% Gaussian prior
gauFit = gauF(vVec,0,sigPFgs)/sum(gauF(vVec,0,sigPFgs));


%% Plot prior

titleC = regexp(fName,'_','split');
temp = titleC{2};

contInd = contL;

% Semilog
f1 = figure;
f1.Position = [100 100 625 525];
hold on;
plot(vVec,pFVec,'--r','linewidth',3);
scatter(vLog,pFVec(2:2:12),100,'r','filled');
plot(vVec,gauFit,'--k','linewidth',3);
set(gca,'xtick',vLog,'xticklabel',vLab,'yscale','lin','ylim',[0 5],...
    'xlim',[vStarts(1) vEnds(end)],'fontsize',20);
ylabel('p(vel)');
xlabel('Stimulus Velocity');
title([titleC{1},' (',titleC{2},'; cont = ',cLab{contInd},')']);

for i = 1:6
    
    thisSig = likeSigs(i,contInd);
    
    xAx = linspace(vLog(i) - 3*thisSig,vLog(i) + 3*thisSig,50);
    pF  = gauF(xAx,vLog(i),likeSigs(i,contInd));
    plot(xAx,pF/(sum(pF)*diff(xAx(1:2))),'b','linewidth',2);
    
end

legend({'Piecewise Fit','Reference Velocities','Gaussian Fit','Likelihoods (piecewise)'});

% Log
f2 = figure;
f2.Position = [745 100 625 525];
hold on;
plot(vVec,pFVec,'--r','linewidth',3);
scatter(vLog,pFVec(2:2:12),100,'r','filled');
plot(vVec,gauFit,'--k','linewidth',3);
set(gca,'xtick',vLog,'xticklabel',vLab,'yscale','log','ylim',[0 5],...
    'xlim',[vStarts(1) vEnds(end)],'fontsize',20);
ylabel('p(vel)');
ylabel('p(vel)');
xlabel('Stimulus Velocity');
title([titleC{1},' (',titleC{2},')']);

% for i = 1:6
%     
%     thisSig = likeSigs(i,contInd);
%     
%     xAx = linspace(vLog(i) - 3*thisSig,vLog(i) + 3*thisSig,50);
%     pF  = gauF(xAx,vLog(i),likeSigs(i,contInd));
%     plot(xAx,pF/(sum(pF)*diff(xAx(1:2))),'b','linewidth',2);
%     
% end

% legend({'Piecewise Fit','Reference Velocities','Likelihoods'},'location','southwest');
legend({'Piecewise Fit','Reference Velocities','Gaussian Fit'},'location','southwest');


%% Plot likelihood pars

gVec = [g g_gs];
hVec = [hc hc_gs];
    
if regexp(fName,'con') ~= 0
    yming = 0;
    yminhc = 0;
else
    negIndsg  = gVec < 0;
    negIndshc = hVec < 0;
    
    if sum(negIndsg) ~= 0
        yming = 1.25*min(gVec(negIndsg));
    else
        yming = 0;
    end
    if sum(negIndshc) ~= 0 
        yminhc = 1.25*min(hVec(negIndshc));
    else
        yminhc = 0;
    end
end

f3 = figure;
f3.Position = [1200 100 825 900];
subplot(2,1,1);
hold on;
scatter(log(c),hc,100,'k','filled');
scatter(log(c),hc_gs,100,'k');
set(gca,'xlim',log([c(1) c(end)]),'xtick',log(c),'xticklabel',cLab,'ylim',[yminhc 1.25*max(hVec)],...
    'fontsize',20);
title([titleC{1},' (',titleC{2},')']);
xlabel('Contrast (%)');
ylabel('h(c) (AU)');
legend({'Piecewise','Gaussian'});

subplot(2,1,2);
hold on
scatter(vLog,g,100,'k','filled');
scatter(vLog,g_gs,100,'k');
set(gca,'xlim',[vLog(1) vLog(end)],'xtick',vLog,'xticklabel',vLab,'ylim',[yming 1.25*max(gVec)],'fontsize',20);
xlabel('Speed (\circ/s)');
ylabel('g(v) (AU)');

end