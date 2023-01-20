clear all
close all

%% Define simulation pars

userDef = 1;
restrictMeasNoise = 1;
usePsychFxnApprox = 1;

figDir = '/media/tyler/Data/MATLAB/cooperLab/2-Modeling_Simulations/BayesModelComp/scratchFigs/testNumeric/';

numObs = 1;

% Prior/Likelihood
plotFig1 = 1;
% Single stimulus posteriors
plotFig2 = 2;
% Joint posteriors
plotFig3 = 1;
% Tbayes
plotFig4 = 2;
% Psychometric function
plotFig5 = 1;

for nOb = 1:numObs

disp(['Running Observer: ',num2str(nOb),'/',num2str(numObs)]);
    
if restrictMeasNoise
    figSuffix = num2str(nOb);
else
    figSuffix = [num2str(nOb),'unrestr'];
end

%% Define observer pars

% Prior
% priorType = 'single Gaussian';
priorType = 'MoG';

% pmvAna = 0;
pmvAna = 1;

if userDef
    
    % Measurement noise
    sig1 = 0.4;
    sig2 = 0.2;
    mu1  = 0.2;

    switch priorType
        case 'MoG'
%                     gam = [0.5 0.75 1.5];
%                     muP  = [0 0 0];
%                     w   = [0.3 0.3 0.4];
%             
                    gam = [1.5 0.5 0.5];
                    muP  = [0 0.9 -0.9];
                    w   = [0.1 0.45 0.45];
            
%             gam = [0.5 0.36];
%             muP = [0.5 -0.5];
%             w   = [0.5 0.5];
            
%             gam = [0.5 0.5];
%             muP = [0 0];
%             w   = [0.5 0.5];
            
        case 'single Gaussian'
            gam = [0.5];
            muP = [0];
            w   = [1];
            
    end
    
    numComp = numel(gam);
    
else
    % Create random observer/stimulus scenario to test equations
    if restrictMeasNoise
        sig1 = 0.5*rand(1);
        sig2 = 0.5*rand(1);
        mu1  = rand(1);
        
        numComp = randi(5,1);
        gam     = 0.5*rand(1,numComp)+0.5;
        muP     = rand(1,numComp);
        w       = rand(1,numComp);
        
        w = w/sum(w);
    else
        sig1 = rand(1);
        sig2 = rand(1);
        mu1  = rand(1);
        
        numComp = randi(5,1);
        gam     = 2*rand(1,numComp);
        muP     = rand(1,numComp);
        w       = rand(1,numComp);
        
        w = w/sum(w);
    end
end


%% Define simulation pars

numRuns = 120000;


%% Define support

minVal = -5*max(gam);
maxVal = 5*max(gam);

numVals = 7500;
x = linspace(minVal,maxVal,numVals);
dx = diff(x(1:2));


%% Define mog prior + cumulative

% prior = zeros(1,numVals);
% 
% for ii = 1:numel(gam)
%    
%     % density function is cut in half when zero-centered, but this doesn't
%     % always work... FIX!
%     prior = prior + w(ii)*normpdf(x,muP(ii),gam(ii));
%     
% end
% 
% priorPDF = prior/sum(prior);

% cPrior = cumsum(priorPDF)';


%% Generate random samples from prior

% if ~pmvAna
% % Need to do this numerically unfortunately - no clear solution for inverse
% % of a sum of cumulative normals
% sampP = rand(1,numRuns);
% 
% [~,minInds] = min(abs(sampP - cPrior));
% 
% priorSamps = x(minInds);
% 
% % check
% f1 = figure;
% f1.Position = [100 100 650 600];
% hold on
% histogram(priorSamps,'Normalization','pdf');
% plot(x,prior,'r','linewidth',2);
% set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20);
% xlabel('x');
% ylabel('p(v)');
% title('Distribution of prior samples');
% end


%% Find bivariate PDF p(m,v) for each stimulus contrast

nBins = 30;
nBinsrs = 70;
edges = linspace(minVal,maxVal,nBins+1);
edgesrs = linspace(minVal,maxVal,nBinsrs+1);
binCents = edges(1:end-1)+diff(edges(1:2));
binCentsrs = linspace(binCents(1),binCents(end),nBinsrs);
dBC = diff(binCentsrs(1:2));

[binCentsmx,binCentsmy] = meshgrid(binCents,binCents);
[binCentsrsmx,binCentsrsmy] = meshgrid(binCentsrs,binCentsrs);

if ~pmvAna
    
%     % Generate measurements from samples
%     % should be fine to use the same prior samples, right?
%     measSamps1 = priorSamps + sig1*randn(1,numRuns);
%     measSamps2 = priorSamps + sig2*randn(1,numRuns);
%     
%     [hc1]  = histcounts2(priorSamps,measSamps1,edges,edges);
%     [hc2]  = histcounts2(priorSamps,measSamps2,edges,edges);
%     
%     jProb1 = hc1'/sum(hc1(:));
%     jProb2 = hc2'/sum(hc2(:));
%     
%     % Linearly interpolate to denser sampling
%     jProb1 = interp2(binCentsmx,binCentsmy,jProb1,binCentsrsmx,binCentsrsmy);
%     jProb2 = interp2(binCentsmx,binCentsmy,jProb2,binCentsrsmx,binCentsrsmy);
    
else
    % Just calculate posterior via pointwise multiplication of analytically
    % derived matrices of likelihood and prior
    switch priorType
        case 'single Gaussian'
            jProb1 = normpdf(binCentsrsmx-binCentsrsmy,0,sig1).*normpdf(binCentsrsmx,muP,gam);
            jProb2 = normpdf(binCentsrsmx-binCentsrsmy,0,sig2).*normpdf(binCentsrsmx,muP,gam);
            
            postProb1 = normpdf(binCentsrsmx-binCentsrsmy,0,sig1).*normpdf(binCentsrsmy,muP,gam)./repmat(sum(jProb1,2),[1 nBinsrs]);
            postProb2 = normpdf(binCentsrsmx-binCentsrsmy,0,sig2).*normpdf(binCentsrsmy,muP,gam)./repmat(sum(jProb2,2),[1 nBinsrs]);
            
        case 'MoG'
            prior2 = zeros(1,numel(binCentsrs));
            
            for ii = 1:numel(gam)
                prior2 = prior2 + w(ii)*normpdf(binCentsrs,muP(ii),gam(ii));    
            end
            
            jProb1 = normpdf(binCentsrsmy-binCentsrsmx,0,sig1).*repmat(prior2,[nBinsrs 1]);
            jProb2 = normpdf(binCentsrsmy-binCentsrsmx,0,sig2).*repmat(prior2,[nBinsrs 1]);
%             postProb1 = normpdf(binCentsrsmy-binCentsrsmx,0,sig1).*repmat(prior2',[1 nBinsrs])./repmat(sum(jProb1,2),[1 nBinsrs]);
%             postProb2 = normpdf(binCentsrsmy-binCentsrsmx,0,sig2).*repmat(prior2',[1 nBinsrs])./repmat(sum(jProb2,2),[1 nBinsrs]);

            postProb1 = normpdf(binCentsrsmy-binCentsrsmx,0,sig1).*repmat(prior2',[1 nBinsrs]);
            postProb2 = normpdf(binCentsrsmy-binCentsrsmx,0,sig2).*repmat(prior2',[1 nBinsrs]);
            
            % normalize each column individually
            postProb1 = postProb1./(repmat(sum(postProb1,1),[nBinsrs 1])*dBC);
            postProb2 = postProb2./(repmat(sum(postProb2,1),[nBinsrs 1])*dBC);
            
    end
end

mu1 = 0.3;
mu2 = linspace((1/3)*minVal,(1/3)*maxVal,20) + mu1;


%% INSTEAD OF JOIN P(M,V), CALCULATE posterior P(X|M)
for ii = 1:nBinsrs
    
    thism = binCentsrs(ii);
    
    % for this measurement, get the modified parameters
    [alphas1,muTildes1,wTildes1] = getMogPostPars(w,gam,muP,sig1,thism);
    [alphas2,muTildes2,wTildes2] = getMogPostPars(w,gam,muP,sig2,thism);
    
    % for each component
    for jj = 1:numel(gam)
        temp1(jj,:) = wTildes1(jj)*normpdf(binCentsrs,alphas1(jj)*thism + muTildes1(jj),sqrt(alphas1(jj)*sig1^2));
        temp2(jj,:) = wTildes2(jj)*normpdf(binCentsrs,alphas2(jj)*thism + muTildes2(jj),sqrt(alphas2(jj)*sig2^2));
    end
    
    postProb1ana(ii,:) = sum(temp1,1);
    postProb2ana(ii,:) = sum(temp2,1);
    
end


% figure;
% subplot(2,2,1);
% imagesc(edgesrs,edgesrs,jProb1);
% set(gca,'fontsize',20,'xlim',[x(1) x(end)],'ylim',[x(1) x(end)],'plotboxaspectratio',[1 1 1],'ydir','normal');
% title('num posterior, stim1');
% subplot(2,2,2);
% imagesc(edgesrs,edgesrs,jProb1ana);
% set(gca,'fontsize',20,'xlim',[x(1) x(end)],'ylim',[x(1) x(end)],'plotboxaspectratio',[1 1 1],'ydir','normal');
% title('ana posterior, stim1');
% subplot(2,2,3);
% imagesc(edgesrs,edgesrs,jProb2);
% set(gca,'fontsize',20,'xlim',[x(1) x(end)],'ylim',[x(1) x(end)],'plotboxaspectratio',[1 1 1],'ydir','normal');
% title('num posterior, stim2');
% subplot(2,2,4);
% imagesc(edgesrs,edgesrs,jProb2ana);
% set(gca,'fontsize',20,'xlim',[x(1) x(end)],'ylim',[x(1) x(end)],'plotboxaspectratio',[1 1 1],'ydir','normal');
% title('ana posterior, stim2');
% 
% figure;
% subplot(1,2,1);
% imagesc(edgesrs,edgesrs,jProb1 - jProb1ana);
% set(gca,'fontsize',20,'xlim',[x(1) x(end)],'ylim',[x(1) x(end)],'plotboxaspectratio',[1 1 1],'ydir','normal');
% title('num posterior - analytic, stim1');
% subplot(1,2,2);
% imagesc(edgesrs,edgesrs,jProb2 - jProb2ana);
% set(gca,'fontsize',20,'xlim',[x(1) x(end)],'ylim',[x(1) x(end)],'plotboxaspectratio',[1 1 1],'ydir','normal');
% title('num posterior - analytic, stim2');


%% Sample from rows of p(m|v) to get p(x|m) and calculate T_bayes for each
%  pair of m1,m2

px2gtx1 = nan(nBinsrs,nBinsrs);

% integrate above unity line
[xm,ym] = meshgrid(binCentsrs,binCentsrs);

posInds = ym>xm;
uniInds = ym==xm;

for ii = 1:nBinsrs
    for jj = 1:nBinsrs

        % Get bivariate posterior for m1,m2 via outer product
        bivPost = postProb2ana(ii,:)'*postProb1ana(jj,:);
        
        % Normalize, since we're feeding this into px2gtx1 and calculating
        % when p(x2>x1) > p(x2>1) which is mutually exclusive and total prob must sum to 1 
        bivPost = bivPost/sum(bivPost(:));
        
        px2gtx1(ii,jj) = sum(bivPost(posInds),'omitnan') + sum(bivPost(uniInds),'omitnan')/2;
        
    end
end

% Find values of m1,m2 where px2 > px1 (2AFC posterior)
Tbayes = px2gtx1 > 0.5;

% Do it analytically to test Eqn 23
switch priorType
    
    case 'single Gaussian'
        alpha1 = (gam^2)/(gam^2 + sig1^2);
        alpha2 = (gam^2)/(gam^2 + sig2^2);
        
        px2gtx1ana = normcdf((alpha2*binCentsrsmy - alpha1*binCentsrsmx)/sqrt(alpha1*sig1^2 + alpha2*sig2^2));
        Tbayesana  = px2gtx1ana > 0.5;
        
    case 'MoG'
        px2gtx1ana = nan(nBinsrs,nBinsrs);
        
        for ii = 1:nBinsrs
            for jj = 1:nBinsrs
                thism1 = binCentsrs(ii);
                thism2 = binCentsrs(jj);
                px2gtx1ana(jj,ii) = get2AFCpostMoG(thism1,thism2,sig1,sig2,gam,muP,w);
            end
        end
        Tbayesana  = px2gtx1ana > 0.5;
end

[JSD] = getJSDiv(px2gtx1ana(:)',px2gtx1(:)');

% Find deviation between decision boundary found numerically and estimated
% with a single line according to our approximation

% Estimate numerical decision boundary as mean of indices spanning the
% boundary
[bndLow,~] = find(diff(Tbayes,[],1)==1);
bndHi = bndLow+1;

dBoundy = mean([binCentsrs(bndLow);binCentsrs(bndHi)]);
dBoundx = binCentsrs;

dBoundyAna = (wTildes1'*alphas1)/(wTildes1'*alphas1)*dBoundx;

figDecision = figure;
figDecision.Position = [100 100 1075 500];

subplot(1,2,1);
hold on;
plot(dBoundx,dBoundy,'k','linewidth',2);
plot(dBoundx,dBoundyAna,'r','linewidth',2);
ylabel('m_{2}');
xlabel('m_{1}');
set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1],'xlim',[-3 3]);

subplot(1,2,2);
plot(dBoundx,dBoundy-dBoundyAna,'r','linewidth',2);
ylabel('Approx. Error (Num-Ana)');
xlabel('m_{1}');
set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1],'xlim',[-3 3]);


%% Make an example psychometric function from numeric Tbayes and analytic Tbayes

% Stimulus 1 and 2 speeds
mu2 = linspace((1/3)*minVal,(1/3)*maxVal,20) + mu1;

pPsychNum = nan(numel(mu2),1);

for ii = 1:numel(mu2)
    
    % likelihoods
    thismu1 = mu1;
    thismu2 = mu2(ii);
    
    % bivariate likelihood
    pm1m2{ii} = normpdf(binCentsrs,thismu2,sig2)'*normpdf(binCentsrs,thismu1,sig1);
    
    % classify and sum to get psych function
    classLike  = pm1m2{ii}.*Tbayes;
    pPsychNum(ii) = sum(classLike(:))*dBC*dBC;
    
    classLikeAna  = pm1m2{ii}.*Tbayesana;
    pPsychAna(ii) = sum(classLikeAna(:))*dBC*dBC;
    
end


%% Calculate psych function using our APPROXIMATE analytical approach

if usePsychFxnApprox
    pPsychAna = nan(numel(mu2),1);
    
    for ii = 1:numel(mu2)
        
        switch priorType
            case 'MoG'
                pPsychAna(ii) = getPsychFxnMoG(mu1,mu2(ii),sig1,sig2,gam,muP,w);
                
            case 'single Gaussian'
                pPsychAna(ii) = getPsychFxnSG(mu1,mu2(ii),sig1,sig2,gam);
                
        end
        
    end
    
end

%% Plotting/Saving

if plotFig1
priorLikeFig = figure;
priorLikeFig.Position = [100 100 1500 440];

tab = uitable(priorLikeFig,'Data',[w',gam',muP'],'ColumnName',{'w','gam','muP'});
pos = get(subplot(1,3,3),'position');
delete(subplot(1,3,3));
set(tab,'units','normalized','position',pos);

subplot(1,3,1);
plot(binCentsrs,prior2,'k','linewidth',2);
title(['MoG prior (',num2str(numComp), ' comps)']);
xlabel('x');
ylabel('p(x)');
set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1]);

subplot(1,3,2);
hold on
plot(x,normpdf(x,mu1,sig1),'r','linewidth',2);
plot(x,normpdf(x,1.5*mu1,sig2),'b','linewidth',2);
title('Gaussian likelihoods');
xlabel('m');
ylabel('p(m|x)');
set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1]);
ax = gca;
yMax = ax.YLim(2);
xMin = ax.XLim(1);
text(0.9*xMin,0.9*yMax,['\sigma_{1} = ',num2str(round(sig1,2))],'fontsize',20);
text(0.9*xMin,0.75*yMax,['\sigma_{2} = ',num2str(round(sig2,2))],'fontsize',20);
text(0.9*xMin,0.6*yMax,['\mu_{ref} = ',num2str(round(mu1,2))],'fontsize',20);

if ~userDef
saveas(priorLikeFig,[figDir,'1-priorLike',figSuffix,'.png']);
end

end

if plotFig2
% check point estimate posteriors
f2 = figure;
f2.Position = [300 100 1625 835];
hold on
subplot(2,3,1);
imagesc(edgesrs,edgesrs,postProb1');
cbar = colorbar;
cbar.Label.String = 'p(x|m)';
set(gca,'fontsize',20,'xlim',[x(1) x(end)],'ylim',[x(1) x(end)],'plotboxaspectratio',[1 1 1],'ydir','normal');
xlabel('x');
ylabel('m');
title(['\sigma_{1} = ',num2str(sig1),' num']);

subplot(2,3,2);
imagesc(edgesrs,edgesrs,postProb1ana);
cbar = colorbar;
cbar.Label.String = 'p(x|m)';
set(gca,'fontsize',20,'xlim',[x(1) x(end)],'ylim',[x(1) x(end)],'plotboxaspectratio',[1 1 1],'ydir','normal');
xlabel('x');
ylabel('m');
title(['\sigma_{1} = ',num2str(sig1),' ana']);

subplot(2,3,4);
imagesc(edgesrs,edgesrs,postProb2');
cbar = colorbar;
cbar.Label.String = 'p(x|m)';
set(gca,'fontsize',20,'xlim',[x(1) x(end)],'ylim',[x(1) x(end)],'plotboxaspectratio',[1 1 1],'ydir','normal');
xlabel('x');
ylabel('m');
title(['\sigma_{2} = ',num2str(sig2),' num']);

subplot(2,3,5);
imagesc(edgesrs,edgesrs,postProb2ana);
cbar = colorbar;
cbar.Label.String = 'p(x|m)';
set(gca,'fontsize',20,'xlim',[x(1) x(end)],'ylim',[x(1) x(end)],'plotboxaspectratio',[1 1 1],'ydir','normal');
xlabel('x');
ylabel('m');
title(['\sigma_{2} = ',num2str(sig2),' ana']);

subplot(2,3,3);
imagesc(edgesrs,edgesrs,postProb1' - postProb1ana);
set(gca,'fontsize',20,'xlim',[x(1) x(end)],'ylim',[x(1) x(end)],'plotboxaspectratio',[1 1 1],'ydir','normal');
title('num - ana, stim1');
cbar = colorbar;
cbar.Label.String = 'p(x|m) (num-ana)';
subplot(2,3,6);
imagesc(edgesrs,edgesrs,postProb2' - postProb2ana);
set(gca,'fontsize',20,'xlim',[x(1) x(end)],'ylim',[x(1) x(end)],'plotboxaspectratio',[1 1 1],'ydir','normal');
title('num - ana, stim2');
cbar = colorbar;
cbar.Label.String = 'p(x|m) (num-ana)';
end

if plotFig3
% Plot numerical & analytical 2AFC posteriors and compare
f3 = figure;
f3.Position = [100 100 1750 415];
subplot(1,3,1);
imagesc(binCentsrs,binCentsrs,px2gtx1)
set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1],'ydir','normal',...
    'ylim',[binCentsrs(1) binCentsrs(end)],'xlim',[binCentsrs(1) binCentsrs(end)]);
title('p(x_{2}>x_{1}|m_{1},m_{2}) (Num)');

subplot(1,3,2);
imagesc(binCentsrs,binCentsrs,px2gtx1ana)
set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1],'ydir','normal',...
    'ylim',[binCentsrs(1) binCentsrs(end)],'xlim',[binCentsrs(1) binCentsrs(end)]);
title('p(x_{2}>x_{1}|m_{1},m_{2}) (Ana)');


subplot(1,3,3);
hold on
imagesc(binCentsrs,binCentsrs,px2gtx1ana-px2gtx1);
set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1],'ydir','normal',...
    'ylim',[binCentsrs(1) binCentsrs(end)],'xlim',[binCentsrs(1) binCentsrs(end)]);
% title('Ana - Num p(x_{2}>x_{1}|m_{1},m_{2})');
title(['JSD = ',num2str(JSD)]);
% ax = gca;
% yMax = ax.YLim(2);
% xMin = ax.XLim(1);
% text(0.9*xMin,09*yMax,['JSD = ',num2str(JSD)],'fontsize',15);
cbar = colorbar;
cbar.Label.String = 'Ana - Num';
if ~userDef
saveas(f3,[figDir,'3-anaNum2AFCPosts',figSuffix,'.png']);
end
end

if plotFig4
% Plot Tbayes analytical and numerical as well as example 2D likelihood
f4 = figure;
f4.Position = [100 100 1750 415];

subplot(1,3,1);
hold on
imagesc(binCentsrs,binCentsrs,Tbayes);
plot([minVal maxVal],(wTildes1'*alphas1)/(wTildes2'*alphas2)*[minVal maxVal],'--r','linewidth',5);
set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1],'ydir','normal',...
    'ylim',[binCentsrs(1) binCentsrs(end)],'xlim',[binCentsrs(1) binCentsrs(end)]);
title('T_{Bayes} Num');

subplot(1,3,2);
hold on
imagesc(binCentsrs,binCentsrs,Tbayesana);
% plot([minVal maxVal],(wTildes1'*alphas1)/(wTildes2'*alphas2)*[minVal maxVal],'--r','linewidth',5);
set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1],'ydir','normal',...
    'ylim',[binCentsrs(1) binCentsrs(end)],'xlim',[binCentsrs(1) binCentsrs(end)]);
title('T_{Bayes} Ana');    
    
    
subplot(1,3,3);
imagesc(binCentsrs,binCentsrs,pm1m2{round(end/2)});
set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,'ylim',[binCentsrs(1) binCentsrs(end)],...
    'ydir','normal','xlim',[binCentsrs(1) binCentsrs(end)]);
ylabel('p(m_{2}|x_{2})');
xlabel('p(m_{1}|x_{1})');
title(['\mu_{1} = ',num2str(mu1),', \mu_{2} = ',num2str(mu2(round(end/2)))]);
end

if plotFig5
% Compare analytical and numerical psychometric functions
f5 = figure;
f5.Position = [100 500 650 600];
hold on;

plot(mu2,pPsychNum,'k','linewidth',2);
set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,'ylim',[0 1],'xlim',[mu2(1) mu2(end)],'ydir','normal');
xlabel('x_{2}');
ylabel('p(R=2|x_{1},x_{2})');
title(['x_{1} = ',num2str(mu1),'; ','\sigma_{1} = ',num2str(sig1),', \sigma_{2} = ',num2str(sig2)]);

plot(mu2,pPsychAna,'--r','linewidth',2);
plot([mu2(1) mu1],[0.5 0.5],'--k');
plot(mu1*[1 1],[0 0.5],'--k');
legend({'Numerical','Analytical'},'Location','southeast');

if ~userDef
saveas(f5,[figDir,'5-anaNumPfxn',figSuffix,'.png']);
end
end

if numObs > 1
allVars = who;

for av = 1:numel(allVars)
    
    m1 = strcmp('userDef',allVars{av});
    m2 = strcmp('restrictMeasNoise',allVars{av});
    m3 = strcmp('usePsychFxnApprox',allVars{av});
    m4 = strcmp('figDir',allVars{av});
    m5 = strcmp('numObs',allVars{av});
    m6 = strcmp('plotFig1',allVars{av});
    m7 = strcmp('plotFig2',allVars{av});
    m8 = strcmp('plotFig3',allVars{av});
    m9 = strcmp('plotFig4',allVars{av});
    m10 = strcmp('plotFig5',allVars{av});
    m11 = strcmp('nOb',allVars{av});
    m12 = strcmp('allVars',allVars{av});
    
    matchV = m1 | m2 | m3 | m4 | m5 | m6 | m7 | m8 | m9 | m10 | m11 | m12;
    
    if ~matchV
    clear(allVars{av})
    end
end

close all
clear m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12 matchV allVars av
end

end

%% Helper functions

function pFxn = getPsychFxnMoG(mu1,mu2,sig1,sig2,gam,muP,w)

    %%%%%%%%%%%%% CHECK ME, TRY OTHER DERIVATION (WITH SUM OUTSIDE)?

    [alphas1,muTildes1,wTildes1] = getMogPostPars(w,gam,muP,sig1,mu1);
    [alphas2,muTildes2,wTildes2] = getMogPostPars(w,gam,muP,sig2,mu2);

    denom = sqrt( (sig1^2)*(wTildes1'*alphas1)^2 + (sig2^2)*(wTildes2'*alphas2)^2 );

    Z = ( sum(wTildes2.*(muTildes2 + alphas2*mu2)) - sum(wTildes1.*(muTildes1 + alphas1*mu1)) )/denom;

    pFxn = normcdf(Z);

end

function pFxn = getPsychFxnSG(mu1,mu2,sig1,sig2,gam)

    alpha1 = (gam^2)/(gam^2 + sig1^2);
    alpha2 = (gam^2)/(gam^2 + sig2^2);

    denom = sqrt( (alpha1^2)*(sig1^2) + (alpha2^2)*(sig2^2) );

    Z = (mu2*alpha2 - mu1*alpha1)/denom;

    pFxn = normcdf(Z);

end

function pFxn = get2AFCpostMoG(m1,m2,sig1,sig2,gam,muP,w)

    ind = 1;
    
    [alphas1,muTildes1,wTildes1] = getMogPostPars(w,gam,muP,sig1,m1);
    [alphas2,muTildes2,wTildes2] = getMogPostPars(w,gam,muP,sig2,m2);
            
    for ii = 1:numel(gam)
        for jj = 1:numel(gam)
            
            denom = sqrt( (sig1^2)*alphas1(jj) + (sig2^2)*alphas2(ii) );
            
            Z = ( muTildes2(ii) + alphas2(ii)*m2 - muTildes1(jj) - alphas1(jj)*m1 )/denom;
            
            pFxnInd(ind) = wTildes1(jj)*wTildes2(ii)*normcdf(Z);
            
            ind = ind + 1;
            
        end
    end
    
    pFxn = sum(pFxnInd);
    
end


function pFxn = get2AFCpostSG(mu1,mu2,sig1,sig2,gam)

    alpha1 = (gam^2)/(gam^2 + sig1^2);
    alpha2 = (gam^2)/(gam^2 + sig2^2);

    denom = sqrt( (alpha1)*(sig1^2) + (alpha2)*(sig2^2) );

    Z = (mu2*alpha2 - mu1*alpha1)/denom;

    pFxn = normcdf(Z);

end










