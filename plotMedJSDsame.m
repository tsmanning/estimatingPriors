% Plot median JSD of all within-prior pairs across prior parameter and
% trial count (i.e. ACCURACY of fits)

close all
clear all

% model = 'ss';
model = 'mg';

% priorDist = 'expon';
priorDist = 'OneGauss';
% priorDist = 'MixGauss';

dim = '';
% dim = 'sig';
% dim = 'mu';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

topDir = '/media/tyler/Data/MATLAB/cooperLab/2-Modeling_Simulations/BayesModelComp/';

load([topDir,'SimData/closed/',model,priorDist,'.mat'])

numResamps = 100;

numPairs = numel(runStruct(1,1,1).JSDdist);

for ii = 1:numResamps
    
    resampMeds(:,:,ii)  = squeeze(arrayfun(@(x) median(datasample(x.JSDdist,numPairs),'omitnan'),runStruct));
    
end

q25 = quantile(resampMeds,0.25,3);
q75 = quantile(resampMeds,0.75,3);

medsTrue = squeeze(arrayfun(@(x) median(x.JSDdist,'omitnan'),runStruct));

switch priorDist
    case 'expon'
        medinds = 1:size(runStruct,1);
        
        pars = squeeze(arrayfun(@(x) x.avlog(1),runStruct));
        pars = pars(medinds,1);
        
        xlab = 'Prior log-slope';
        
    case 'OneGauss'
        medinds = 1:size(runStruct,1);
        
        pars = squeeze(arrayfun(@(x) x.sigP(1),runStruct));
        pars = pars(medinds,1);
        
        xlab = 'Prior \sigma';
        
    case 'MixGauss'
        switch dim
            case 'sig'
                medinds = 1:5;
                
                pars = squeeze(arrayfun(@(x) x.sigP(1),runStruct));
                pars = pars(medinds,1);
                
                xlab = 'Prior \sigma';
                
            case 'mu'
                medinds = [3 6:9];
                
                pars = squeeze(arrayfun(@(x) x.muP(1),runStruct));
                pars = pars(medinds,1);
                
                xlab = 'Prior \mu';
        end
end

colormat = colororder;

fig = figure;
fig.Position = [100 100 660 600];
hold on;

for ii = 1:size(medsTrue,2)
    
    p{ii} = plot(pars,medsTrue(medinds,ii),'color',colormat(ii,:),'linewidth',4);
    fill([pars' flipud(pars)'],[q25(medinds,ii)' flipud(q75(medinds,ii))'],colormat(ii,:),'FaceAlpha',0.3,'linestyle','none');
    
end

set(gca,'yscale','log','ylim',[1e-4 0.2],'ytick',[1e-4 1e-1],'fontsize',20,'plotboxaspectratio',[1 1 1],...
    'xtick',pars,'xlim',[pars(1) pars(end)]);
xlabel(xlab);
ylabel('Median JSD');
legend([p{1} p{2} p{3} p{4}],{'900 trials','4,200 trials','8,400 trials','16,800 trials'},'location','southeast');

saveas(fig,[topDir,'Figures/VVSSabstract/summaryfigs/',model,priorDist,dim,'.svg']);

