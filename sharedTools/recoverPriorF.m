function [pFVec,vLog] = recoverPriorF(method,priorF,priorW,lr,nTr,contL,tbtOn,interpOn)

% For a set of input slopes, recover estimate of prior shape

%% Load sim data and fits

nummeth = 'closed';
if tbtOn && ~interpOn
    suffix = 'TBT';
elseif tbtOn && interpOn == 1
    suffix = 'TBTinterp';
elseif tbtOn && interpOn == 2
    suffix = 'TBTNN';
elseif ~tbtOn && ~interpOn
    suffix = 'biScl';
end

splPath = regexp(which('recoverPriorF'),filesep,'split');
sDir = [filesep,fullfile(splPath{1:numel(splPath)-2}),filesep,'SimData',filesep,nummeth,filesep];

% method = 'ss';
% priorF = 'gauss';
% priorW = '0.5';
% lr     = '0';
% nTr    = '250';

parEstStr = ['parEsts_',method,'_n',nTr,'_collect',suffix,'.mat'];

% Load GT/exp pars
load([sDir,'gt_',method,'_log_',priorF,'_',priorW,'_LR_',lr,filesep,'simData.mat']);
% Load fir pars
load([sDir,'gt_',method,'_log_',priorF,'_',priorW,'_LR_',lr,filesep,parEstStr]);

% Plot sample with best likelihood
% [~,bInd] = min(nllF);
bInd = find(nllF==median(nllF));

a       = avlogF(bInd,:);
hc      = hcF(bInd,:);
g       = gvlogF(bInd,:);
sigGT   = str2num(priorW);

%% Setup pars

% Likelihoods of sigmas (make it not-speed dependent for now)
likeSigs = g(1)*hc;

% Speeds & such
v       = [0.5 1 2 4 8 12];
vLog    = getLogXform(v,0.3);

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


%% Plot prior

contInd = contL;

% Semilog
f1 = figure;
f1.Position = [100 100 625 525];
hold on;
plot([vLogInit vLogFull],gauF([vLogInit vLogFull],0,sigGT)/(vSpacing*sum(gauF(vLogFull,0,sigGT))),...
    'k','linewidth',3);
plot(vVec,pFVec,'--r','linewidth',3);
scatter(vLog,pFVec(2:2:12),100,'r','filled');
set(gca,'xtick',vLog,'xticklabel',vLab,'yscale','lin','ylim',[0 5],...
    'xlim',[vStarts(1) vEnds(end)],'fontsize',20);
ylabel('p(vel)');
xlabel('Stimulus Velocity');

for i = 1:6
    
    thisSig = likeSigs(contInd);
    
    xAx = linspace(vLog(i) - 3*thisSig,vLog(i) + 3*thisSig,50);
    pF  = gauF(xAx,vLog(i),likeSigs(contInd));
    plot(xAx,pF/(sum(pF)*diff(xAx(1:2))),'b','linewidth',3);
    
end

legend({'Ground Truth','Piecewise Fit','Reference Velocities','Likelihoods'});

% Log
logTicks = [1e-14 1e-12 1e-10 1e-8 1e-6 1e-4 1e-2 1];

f2 = figure;
f2.Position = [745 100 625 525];
hold on;
plot([vLogInit vLogFull],gauF([vLogInit vLogFull],0,sigGT)/(vSpacing*sum(gauF(vLogFull,0,sigGT))),...
    'k','linewidth',3);
plot(vVec,pFVec,'--r','linewidth',3);
scatter(vLog,pFVec(2:2:12),100,'r','filled');
set(gca,'xtick',vLog,'xticklabel',vLab,'yscale','log','ylim',[10e-14 5],...
    'xlim',[vStarts(1) vEnds(end)],'ytick',logTicks,'fontsize',20);
ylabel('p(vel)');
ylabel('p(vel)');
xlabel('Stimulus Velocity');

% for i = 1:6
%     
%     thisSig = likeSigs(contInd);
%     
%     xAx = linspace(vLog(i) - 3*thisSig,vLog(i) + 3*thisSig,50);
%     pF  = gauF(xAx,vLog(i),likeSigs(contInd));
%     plot(xAx,pF/(sum(pF)*diff(xAx(1:2))),'b','linewidth',2);
%     
% end

% legend({'Ground Truth','Piecewise Fit','Reference Velocities'},'location','southwest');

end