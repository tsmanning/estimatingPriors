
close all
clear all

x = getLogXform(0.5:0.01:15,0.3);
expF = @(x,a,b) exp(a*x + b);

a = linspace(-8,-1.5,5);

xtixLin = [0.5 1 2 4 8 12];
xtix = getLogXform(xtixLin,0.3);

for i = 1:numel(xtix)
    xtixLab{i} = num2str(xtixLin(i));
end

% for i = 1:5
%     
%     f{i} = figure;
%     thisfig = f{i};
%     thisfig.Position = [100+100*i 100 635 610];
%     
%     plot(x,expF(x,a(i),0),'k','linewidth',5);
%     
%     set(gca,'ytick',[0 1],'xtick',xtix,'xticklabel',xtixLab,'fontsize',30,'plotboxaspectratio',[1 1 1]);
%     
% end

f6 = figure;
f6.Position = [100+100*i 100 635 610];
hold on

for i = 1:5
    
    prior = expF(x,a(i),0);
    auc = sum( 0.01.*(0.5*(prior(2:end) + prior(1:end-1))));
    prior = prior/auc;
    
    plot(x,prior,'k','linewidth',5);
    
    set(gca,'ytick',[0 1],'xlim',getLogXform([0.5 12],0.3),'xtick',xtix,'xticklabel',xtixLab,'fontsize',30,'plotboxaspectratio',[1 1 1]);
    
end

f7 = figure;
f7.Position = [100+100*i 100 635 610];
hold on

for i = 1:5
    
    prior = expF(x,a(i),0);
    auc = sum( 0.01.*(0.5*(prior(2:end) + prior(1:end-1))));
    prior = prior/auc;
    
    plot(x,prior,'k','linewidth',5);
    
    set(gca,'ytick',[1e-6 1],'ylim',[1e-7 1e1],'xlim',getLogXform([0.5 12],0.3),'xtick',xtix,'xticklabel',xtixLab,...
        'fontsize',30,'plotboxaspectratio',[1 1 1],'yscale','log');
    
end