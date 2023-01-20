function [outmat] = getMinDistEst2(runStruct,JSDpriorsTrCnts,priorDist,dim,domain)

% Estimate the minimum distance apart along one parameter of prior shape
% two priors must be to tell them apart at 95% significance level
%
% Fit JSD for three priors with gamma function and interpolate
% distributions between them

numPriors = size(runStruct,1);
numPairs  = size(JSDpriorsTrCnts,2);
numTrCnts = size(runStruct,3);
numSampPairs = size(JSDpriorsTrCnts{1,1},1);

priorPairs = makeUniquePairs(numPriors);

aW = nan(numPriors,numTrCnts);
bW = nan(numPriors,numTrCnts);
aB = nan(numPairs,numTrCnts);
bB = nan(numPairs,numTrCnts);

switch priorDist 
    
    case 'expon'
        % should include bit here on fitting expon with mixgauss?
        priorPars = squeeze(arrayfun(@(x) x.avlog(1),runStruct));
        priorPars = priorPars(:,1);
        
        xlab = 'Exponential prior log-slope';
        
        theseinds = 1:numPriors;
        thesePairinds = 1:numPriors-1;
        
    case 'mixGauss'
        switch dim
            case 'sig'
%                 theseinds = 1:3;
                theseinds = 1:5;
                thesePairinds = 1:numPriors-1;
                
                priorPars = squeeze(arrayfun(@(x) x.sigP(1),runStruct));
                priorPars = priorPars(theseinds,1);
                xlab = 'Gaussian prior sigma';
                
            case 'mu'
%                 theseinds = [2 4 5];
                theseinds = [1 6 7 8 9];
                thesePairinds = [2 9 16:21];
                
                priorPars = squeeze(arrayfun(@(x) x.muP(1),runStruct));
                priorPars = priorPars(theseinds,1);
                
                xlab = 'Gaussian prior mu';
                
            case 'tail'
%                 theseinds = [2 6 7];
                theseinds = [2 6 7];
                
%                 priorPars = squeeze(arrayfun(@(x) x.muP(1),runStruct));
%                 priorPars = priorPars(theseinds,1);
                priorPars = 1:3;
                
                xlab = 'Gaussian prior tail';
        end
        
    case 'oneGauss'
        priorPars = squeeze(arrayfun(@(x) x.sigP(1),runStruct));
        priorPars = priorPars(:,1);
          
        xlab = 'Gaussian prior sigma';
        
        theseinds = 1:numPriors;
        thesePairinds = 1:numPriors-1;
        
    otherwise
        error('prior distribution undefined.');
        
end


%% Collect data points into single mat

% Within observer pairs
if isrow(runStruct(1,1,1).JSDdist)
    t1 = squeeze(arrayfun(@(x) x.JSDdist',runStruct,'uniformoutput',0));
else
    t1 = squeeze(arrayfun(@(x) x.JSDdist,runStruct,'uniformoutput',0));
end

% Between observer pairs (drop non-canonical pairs for now)
t2 = JSDpriorsTrCnts(:,thesePairinds)';

numResamps = 100;

fita2 = nan(numTrCnts,numResamps);
fitPrPar = nan(numTrCnts,numResamps);

for bb = 1:numResamps + 1
    
    t1 = t1(theseinds(1),:);
    
    if bb == numResamps + 1
        t1rs = t1;
        t2rs = t2;
    else
        t1rs = cellfun(@(x) datasample(x,numSampPairs),t1,'uniformoutput',0);
        t2rs = cellfun(@(x) datasample(x,numSampPairs),t2,'uniformoutput',0);
    end
    
    % Concat data into desired pairs/trial number counts
    % First make an m x n cell array (m=#priors, n=#trial counts)
    theseData = [t1rs;t2rs];
    
    % Linearize array so it's now column vector of
    % [{set of priors, trcnt1};{set of priors trcnt2};...]
    % Then convert to mat which also expands this [mxn,1] array to [mxnxp,1]
    % where p is the number of elements in JSD distribution
    theseData = cell2mat(theseData(:));

    % Keep track of the number of what inds contribute to each fit parameter
    [aW,bW] = fitGamDist(theseData,numPriors,numSampPairs,numTrCnts);
    
    aW = reshape(aW,[numPriors numTrCnts]);
    
    
    %% Numerically solve for minimum value of a2 that maintains 95% significance
    
    dx = 0.001;
    x = 0:dx:1;
    b1 = bW;
    b2 = bW;
    
    
    
    for ii = 1:numTrCnts
        
        a1 = aW(1,ii);
        gammaSepFxn = @(a2) abs( 0.95 - sum( gampdf(x,a2,b2(ii)).*gamcdf(x,a1,b1(ii)) ,'omitnan')*dx );
        
        opts      = optimset('display','off');
        initVal   = 1;
        
        fita2(ii,bb) = fminunc(gammaSepFxn,initVal,opts);
        pVals(ii) = sum( gampdf(x,fita2(ii),b2(ii)).*gamcdf(x,a1,b1(ii)) )*dx;
        
    end
    
    
    %% Relate a2 parameter back to prior distribution parameters
    
    for ii = 1:numTrCnts
        
        switch domain
            case 'log'
                fitPrPar(ii,bb) = interp1(log(aW(theseinds,ii)),priorPars,log(fita2(ii,bb)));
            case 'lin'
                fitPrPar(ii,bb) = interp1(aW(theseinds,ii),priorPars,fita2(ii,bb));
        end
        
    end
    
end

fitPrParGT = fitPrPar(:,end);
fita2GT    = fita2(:,end);

fitPrParRS = fitPrPar(:,1:end-1);
fita2RS    = fita2(:,1:end-1);

%% Check plots

if 0
colorMat = colororder;

fig = figure;
fig.Position = [100 100 810 670];
hold on;

switch domain
    case 'lin'
        ylim = [0 10];
        loca = 'northeast';
    case 'log'
        ylim = [-2 3];
        loca = 'southeast';
end

for ii = 1:numTrCnts

    switch domain
        case 'log'
            plotdata1 = log(aW(theseinds,ii));
            fitdata = log(fita2(ii));
        case 'lin'
            plotdata1 = aW(theseinds,ii);
            fitdata = fita2(ii);
    end

%     plotdata1 = aW(:,ii) - aW(1,ii) + 0.1;
%     plotdata1 = log(aW(:,ii) - aW(1,ii) + 0.1);

    p{ii} = plot(priorPars,plotdata1,'color',colorMat(ii,:),'linewidth',3);
    
    plot([priorPars(1) priorPars(end)],fitdata*[1 1],'color',colorMat(ii,:),'linestyle','--','linewidth',3);
    plot(fitPrPar(ii)*[1 1],[-3 5],'color',colorMat(ii,:),'linestyle','--','linewidth',3);
    
end

xlabel(xlab);
ylabel('JSD_{diff} gamma distribution parameter');
set(gca,'fontsize',20,'xlim',[priorPars(1) priorPars(end)],'xtick',priorPars);
legend([p{1},p{2},p{3},p{4}],{'900 trials','4,200 trials','8,400 trials','16,800 trials'},'location',loca);
% legend([p{1},p{2},p{3}],{'4,200 trials','8,400 trials','16,800 trials'},'location',loca);
end

%% Make output matrix
outmat.aWGT = aW;
outmat.bWGT = bW;
outmat.fitGamDistGT = fita2GT;
outmat.fitPrParGT = fitPrParGT;
outmat.fitGamDistRS = fita2RS;
outmat.fitPrParRS = fitPrParRS;
outmat.domain = domain;

end