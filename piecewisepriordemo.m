%% Generate figure showing piecewise prior estimation

clear all
close all

% Define functions
gauF  = @(x,mu,sigma) (1/(sigma*sqrt(2*pi)))*exp(-0.5*((x-mu)/sigma).^2);
dGauF = @(x,mu,sigma) -((x-mu)/(sqrt(2*pi)*sigma^3)).*exp(-0.5*((x-mu)/sigma).^2);
lineF = @(x,m,a,b) m*(x-a) + b;
expF  = @(x,a,b) exp(a*x+b);

x            = getLogXform(0:0.001:18,0.3);
likeMeansLin = [0.5 1 2 4 8 12];
likeMeansLog = getLogXform(likeMeansLin,0.3);

xlabs = cell(numel(likeMeansLog),1);
for i = 1:numel(likeMeansLog)
    xlabs{i,1} = num2str(likeMeansLin(i));
end

% Make a prior
ws   = [0.6 0.2 0.2];
sigs = [1 0.6 0.6];
mus  = [0 0 1];

priorSlopes = zeros(1,numel(likeMeansLog));
priorVals   = zeros(1,numel(likeMeansLog));
gt          = zeros(1,numel(x));
mogest      = zeros(1,numel(likeMeansLog)+2);

priorStarts     = getLogXform([0.25  0.75 1.5  3  6 10],0.3);
priorEnds       = getLogXform([0.75  1.5    3  6 10 15],0.3);

for ii = 1:numel(ws)
    priorSlopes     = priorSlopes + ws(ii)*dGauF(likeMeansLog,mus(ii),sigs(ii));
    priorVals       = priorVals   + ws(ii)*gauF(likeMeansLog,mus(ii),sigs(ii));
    gt              = gt          + ws(ii)*gauF(x,mus(ii),sigs(ii));
    mogest          = mogest      + ws(ii)*gauF([0 priorStarts priorEnds(end)],mus(ii),sigs(ii));
end

% Get endpoints of piecewise estimates
priorValsStart  = priorSlopes.*(likeMeansLog - priorStarts);
priorValsEnd    = priorSlopes.*(priorEnds - likeMeansLog);
priorValsOff    = mean([priorValsStart;priorValsEnd]);
priorValsStart  = priorValsStart - priorValsOff;
priorValsEnd    = priorValsEnd   - priorValsOff;

%% Plot

% Piecewise
f1 = figure;
f1.Position = [100 100 550 500];
hold on;

% Get colors
colorMap = colororder;

% Plot GT prior
plot(x,gt,'color',colorMap(1,:),'linewidth',4);
set(gca,'xtick',likeMeansLog,'xticklabel',xlabs,'xlim',getLogXform([0 12],0.3),...
    'ylim',[0 0.5],'fontsize',20);
ylabel('p(Velocity)');
xlabel('Velocity (\circ/s)');

% Plot set of piecewise components
for i = 1:6
    plot([priorStarts(i) priorEnds(i)],...
        lineF([priorStarts(i) priorEnds(i)],priorSlopes(i),likeMeansLog(i),priorVals(i)),...
        '*--k','linewidth',6,'markersize',15);
%     plot([priorStarts(i) priorEnds(i)],...
%         expF([priorStarts(i) priorEnds(i)],log(priorSlopes(i)),priorStarts(i)),...
%         '*--k','linewidth',6);
end

legend({'Ground truth prior','Piecewise estimate'})

% MoG
f2 = figure;
f2.Position = [100 100 550 500];
hold on;

% Plot GT prior
plot(x,gt,'color',colorMap(1,:),'linewidth',4);
plot([0 priorStarts priorEnds(end)],mogest,'--k','linewidth',6);
set(gca,'xtick',likeMeansLog,'xticklabel',xlabs,'xlim',getLogXform([0 12],0.3),...
    'ylim',[0 0.5],'fontsize',20);
ylabel('p(Velocity)');
xlabel('Velocity (\circ/s)');

% Plot set of gaussian components
for ii = 1:numel(ws)
    plot(x,ws(ii)*gauF(x,mus(ii),sigs(ii)),'color',[0.3 0.3 0.3],'linestyle','--','linewidth',3);
end

legend({'Ground truth prior','Mixture of Gaussians estimate','Gaussian Components'})
