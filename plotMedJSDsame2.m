% Plot median JSD of all within-prior pairs across prior parameter and
% trial count (i.e. ACCURACY of fits)

close all
clear all

% model = 'ss';
% model = 'mg';

priorDist = 'expon';
% priorDist = 'OneGauss';
% priorDist = 'MixGauss';

dim = '';
% dim = 'sig';
% dim = 'mu';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

topDir = '/media/tyler/Data/MATLAB/cooperLab/2-Modeling_Simulations/BayesModelComp/';

load([topDir,'SimData/closed/','ss',priorDist,'.mat']);
rs1 = runStruct;
load([topDir,'SimData/closed/','mg',priorDist,'.mat']);
rs2 = runStruct;

numResamps = 100;

%%%%%%%%%%%% Calculate JSD between each fit and the ground truth

jsd1 = squeeze(arrayfun(@(x) getJSDiv(repmat(x.gtPrior,[size(x.fitPriors,1) 1]),x.fitPriors),rs1,'uniformoutput',0));
jsd2 = squeeze(arrayfun(@(x) getJSDiv(repmat(x.gtPrior,[size(x.fitPriors,1) 1]),x.fitPriors),rs2,'uniformoutput',0));

numPairs = numel(jsd1{1,1});

for ii = 1:numResamps
    
    resampMeds1(:,:,ii)  = cellfun(@(x) median(datasample(x,numPairs),'omitnan'),jsd1);
    resampMeds2(:,:,ii)  = cellfun(@(x) median(datasample(x,numPairs),'omitnan'),jsd2);
    
end

medsTrue1 = cellfun(@(x) median(x,'omitnan'),jsd1);
medsTrue2 = cellfun(@(x) median(x,'omitnan'),jsd2);

q25_1 = medsTrue1 - quantile(resampMeds1,0.25,3);
q75_1 = quantile(resampMeds1,0.75,3) - medsTrue1;

q25_2 = medsTrue2 - quantile(resampMeds2,0.25,3);
q75_2 = quantile(resampMeds2,0.75,3) - medsTrue2;

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
fig.Position = [100 100 660 650];
hold on;

offset1 = -0.125;
offset2 = 0.125;

for ii = 1:numel(pars)
    tix{ii} = num2str(pars(ii));
end

for ii = 1:size(medsTrue1,2)
    
    thisOffset1 = offset1;
    thisOffset2 = offset2;
    
%     p{ii} = plot(pars,medsTrue(medinds,ii),'color',colormat(ii,:),'linewidth',4);
%     fill([pars' flipud(pars)'],[q25(medinds,ii)' flipud(q75(medinds,ii))'],colormat(ii,:),'FaceAlpha',0.3,'linestyle','none');

%     p1 = errorbar([1:5]+thisOffset1,medsTrue1(medinds,ii),q25_1(medinds,ii),q75_1(medinds,ii),'.','color',colormat(ii,:),'linewidth',3,'MarkerSize',40);
%     p2 = errorbar([1:5]+thisOffset2,medsTrue2(medinds,ii),q25_2(medinds,ii),q75_2(medinds,ii),'o','color',colormat(ii,:),'linewidth',3,'MarkerSize',10);
    
    p1 = scatter([1:5]+thisOffset1,medsTrue1(medinds,ii),300,colormat(ii,:),'filled');
    p2 = scatter([1:5]+thisOffset2,medsTrue2(medinds,ii),300,colormat(ii,:),'linewidth',4);

end

set(gca,'yscale','log','ylim',[1e-5 1],'ytick',[1e-5 1],'fontsize',25,'plotboxaspectratio',[1 1 1],...
    'xtick',1:5,'xticklabels',tix,'xlim',[0.5 5.5]);
xtickangle(30); 
xlabel(xlab);
ylabel('Median JSD');
% legend([p1 p2],{'Piecewise','Mix. of Gauss.'},'location','northwest');

saveas(fig,[topDir,'Figures/VVSSabstract/summaryfigs/','comb',priorDist,dim,'.svg']);

