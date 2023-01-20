% Plot median JSD of all within-prior pairs across prior parameter and
% trial count (i.e. Discriminability)

close all
clear all

% model = 'ss';
% model = 'mg';

% priorDist = 'expon';
% priorDist = 'OneGauss';
priorDist = 'MixGauss';

% dim = '';
% dim = 'sig';
dim = 'mu';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

topDir = '/media/tyler/Data/MATLAB/cooperLab/2-Modeling_Simulations/BayesModelComp/';

load([topDir,'SimData/closed/','ss',priorDist,'.mat'])
outmat1 = outmat;
load([topDir,'SimData/closed/','mg',priorDist,'.mat'])
outmat2 = outmat;

switch priorDist
    case 'expon'
        xlab = 'Log slope';
        xticks = arrayfun(@(x) x.avlog(1),runStruct(:,1,1));
        defMat = outmat;
        xlims = [-8 -1.5];
        
    case 'OneGauss'
        xlab = 'Gaussian \sigma';
        xticks = arrayfun(@(x) x.sigP(1),runStruct(:,1,1));
        defMat = outmat;
        xlims = [0.5 1.5];
        
    case 'MixGauss'
        switch dim
            case 'sig'
                xticks = arrayfun(@(x) x.sigP(1),runStruct(:,1,1));
                xticks = xticks([1:5]);
                xlab = 'Gaussian \sigma';
                defMat = outmat;
                xlims = [0.5 1.5];
                
            case 'mu'
                xticks = arrayfun(@(x) x.muP(1),runStruct(:,1,1));
                xticks = xticks([3 6:9]);
                xlab = 'Gaussian \mu';
                defMat = outmat2;
                xlims = [0 1];
        end
end

ylab = 'Trial count';

c = colororder;

fig = figure;
fig.Position = [100 100 800 750];
hold on;

numTrials = {'900','4200','8400','16800'};
% numTrials = {'4200','8400','16800'};

y1    = outmat1.fitPrPar;
y2    = outmat2.fitPrPar;
nany = find(isnan(y1));

for ii = 1:numel(y1)

    if ~isnan(y1(ii))
        switch priorDist
            case 'expon'
                % ss
                fill([-8 -8 y1(ii) y1(ii)],[ii ii+0.4 ii+0.4 ii],c(ii,:));
                % mg
                fill([-8 -8 y2(ii) y2(ii)],[ii-0.4 ii ii ii-0.4],c(ii,:));
                
%             case 'OneGauss'
%                 fill([0.5 0.5 y(ii) y(ii)],[ii-0.4 ii+0.4 ii+0.4 ii-0.4],c(ii,:));
%                 fill([0.5 0.5 y(ii) y(ii)],[ii-0.4 ii+0.4 ii+0.4 ii-0.4],c(ii,:));
                
            case 'MixGauss'
                switch dim
                    case 'sig'
                        fill([0.5 0.5 y1(ii) y1(ii)],[ii ii+0.4 ii+0.4 ii],c(ii,:));
                        fill([0.5 0.5 y2(ii) y2(ii)],[ii-0.4 ii ii ii-0.4],c(ii,:));
                        
                    case 'mu'
                        fill([0 0 y1(ii) y1(ii)],[ii ii+0.4 ii+0.4 ii],c(ii,:));
                        fill([0 0 y2(ii) y2(ii)],[ii-0.4 ii ii ii-0.4],c(ii,:));
                end
        end
    else
        text((xlims(2)+xlims(1))/2,ii,'N.S. in range','fontsize',30,'horizontalalignment','center');
    end
    
end

set(gca,'xlim',xlims,'xtick',xticks,'ylim',[0 numel(y1)]+0.5,'ytick',1:numel(y1),'yticklabel',numTrials,'fontsize',30);
xlabel(xlab);
ylabel(ylab);

saveas(fig,[topDir,'Figures/VVSSabstract/summaryfigs/','comb',priorDist,dim,'_disc.svg']);
