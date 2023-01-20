% Compare BIC and bias for two different approaches to fitting prior

% clear all
% close all


%% Load them in

numPriors = size(rs1,1);
numLRs    = size(rs1,2);
numTrCnts = size(rs1,3);

colorMat = colororder;

%% Bayesian Information Criterion comparison

offset = [-0.25 -0.125 0.125 0.25];

lrInd = 1;

comp = 3;

switch comp
    case 1
        xlab = 'Prior \sigma';
        xtix = {'0.5','0.75','1.0','1.25','1.5'};
        titles = 'Gaussian Prior \Delta\sigma';
        theseInds = [1:5];
        
    case 2
        xlab = 'Prior \mu';
        xtix = {'0','0.25','0.5','0.75','1.0'};
        titles = 'Gaussian Prior \Delta\mu';
        theseInds = [3 6:9];
        
    case 3
        xlab = 'Prior log-slope';
        xtix = {'-8','-6.375','-4.75','-3.125','-1.5'};
        titles = 'Exponential Prior \Delta log-slope';
        theseInds = [1:5];
end
    
bic1 = squeeze(arrayfun(@(x) x.bic,rs1,'uniformoutput',0));
bic2 = squeeze(arrayfun(@(x) x.bic,rs2,'uniformoutput',0));

bic1b = cell2mat(bic1(:));
bic2b = cell2mat(bic2(:));

deltaBIC = reshape(median(bic1b - bic2b(:,1:size(bic1b,2)),2),numPriors,4);

% Plot scatter
bicFig = figure;
if numPriors == 7
    bicFig.Position = [100 100 950 800];
else
    bicFig.Position = [100 100 630 600];
end
hold on;
      
for ii = 1:5
    
    thisOffset = ii + offset;
    
    
    
    for jj = 1:numTrCnts
          
        scatter(thisOffset(jj),deltaBIC(theseInds(ii),jj),150,colorMat(jj,:),'linewidth',7);
        
    end
    
end

fill([0.5 0.5 5.5 5.5],[-10 10 10 -10],[0 0 0],'FaceAlpha',0.3,'EdgeAlpha',0.3);
    plot([0.5 5.5],[0 0],'k','linewidth',2)
    
set(gca,'xtick',1:5,'xticklabel',xtix,'xlim',[0.5 5.5],'ylim',[-150 100],'fontsize',30,...
            'plotboxaspectratio',[1 1 1]);
        
xtickangle(30);        
xlabel(xlab);
title(titles);

legend({'900 trials','4200 trials','8400 trials','16800 trials'},'location','southeast');

