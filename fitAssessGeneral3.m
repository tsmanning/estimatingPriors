% Fitting analysis

clear all
close all

%% Define Datasets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% Define exp/observer params

% model       = 'ss';
model       = 'mg';
% model       = 'mgfix';
% model       = 'mgfixneg';
% model       = 'zeroMu';
% model       = 'zeroMuWeq';

% priorF      = 'expon';
% priorF      = 'splitExp';
priorF      = 'mixgauss';

% priorInds = 1:5;
% priorInds = 1:7;
% priorInds = 1:9;
priorInds = 1:3;

% dStr = '20210404a';
% dStr = '2021-04-12';
% dStr = '2021-05-01';
% dStr = '2021-05-02';
% dStr = '2021-07-28';
% dStr = '2021-08-12';
% dStr = '2021-08-13';
dStr = '2021-09-01';

% lapseRates  = {'0','0.05'};
lapseRates  = {'0'};

% trCnts      = [5 10 20];
trCnts      = [20];

% Load/sort AND/OR do analysis and plot
loadSort = 1;
analysisPlot = 1;

% Define figure save directory
saveON = 0;
subDir = 'debug';

% Define which saved runs to use
% setDir = '';
% setDir = ['1May2021',filesep];
% setDir = ['1Aug2021_NVBLScorrect',filesep];
% setDir = ['13Aug2021_oldVarBLS',filesep];
% setDir = ['13Aug2021_newVarBLS',filesep];
% setDir = ['16Aug2021_NVBLS',filesep];
% setDir = ['16Aug2021_NVBLS_expon',filesep];
% setDir = ['17Aug2021_Expon',filesep];
setDir = ['1Sep2021_mgNV',filesep];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

domain      = 'log';

switch model
    case 'ss'
        interp      = 2;
    otherwise
        if strcmp(priorF,'expon')
            %%% should really fix this properly in instruction set
            interp  = 2;
        else
            interp  = 0;
        end
end

%%%% fix
interp = 0;

tbtON       = 1;
binScl      = 0;

% method      = 'numeric';
method      = 'closed';

% Define root simulated data directory
splPath = regexp(which('fitAssessGeneral3'),filesep,'split');
workDir = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep,'SimData',filesep,method,filesep,setDir];


%% Setup filenames and pars

numPriors = numel(priorInds);

% Set save dir for figures
figDir = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep,'Figures',filesep,subDir,filesep];   

% Define dirs and filenames to load in
ind = 1;

for i = 1:numel(lapseRates)
    for j = 1:numPriors

        simdat{ind} = ['gt_',domain,'_',priorF,'_pInd_',num2str(j),...
                         '_',num2str(numPriors),'_LR_',lapseRates{i}];
        
        dataSets{ind} = ['gt_',domain,'_',priorF,'_pInd_',num2str(j),...
                         '_',num2str(numPriors),'_LR_',lapseRates{i},...
                         '_',model,'_',dStr];
        
        ind = ind + 1;
    end
end

for k = 1:numel(trCnts)
    if tbtON
        if interp == 1
            parEstNames{k} = ['parEsts_',model,'_n',num2str(trCnts(k)),'_collectTBTinterp'];
        elseif interp == 2
            parEstNames{k} = ['parEsts_',model,'_n',num2str(trCnts(k)),'_collectTBTNN'];
        else
            parEstNames{k} = ['parEsts_',model,'_n',num2str(trCnts(k)),'_collectTBT'];
        end
    else
        if binScl
            parEstNames{k} = ['parEsts_',model,'_n',num2str(trCnts(k)),'_collectbiScl'];
        else
            parEstNames{k} = ['parEsts_',model,'_n',num2str(trCnts(k)),'_collect'];
        end
    end
end

numDatasets = numel(dataSets);
numEstRuns  = numel(parEstNames);
numLRs      = numel(lapseRates);
numTrCnts   = numel(trCnts);

% Get MATLAB default color indexing
colorMat = colororder;
   
% Hardcoding ref stim range for labeling... watch out
vStim1 = [0.5 1 2 4 8 12];
for ii = 1:6
    vStimLab{ii} = num2str(vStim1(ii));
end

%% Load it all in, put it in an index struct
    
if loadSort

prInd   = 1;
lrInd   = 1;

% Loop over simulated datasets
for ds = 1:numDatasets
    
    simDatPath = [workDir,simdat{ds},filesep];
    dsPath = [workDir,dataSets{ds},filesep];
    
    % Loop over estimate runs
    for ne = 1:numEstRuns     
       
        load([simDatPath,'simDataN',num2str(trCnts(k))]);
        load([dsPath,parEstNames{ne}]);
        
        % dim 1: priorInd
        % dim 2: lapse rate
        % dim 3: sample size
        
        switch model
            
            % Piecewise
            case 'ss'
                runStruct(prInd,lrInd,ne).avlogF     = avlogF;
                
                if strcmp(priorF,'expon')
                    runStruct(prInd,lrInd,ne).avlog  = avlog;
                else
                    runStruct(prInd,lrInd,ne).avlog  = avlog;
                    runStruct(prInd,lrInd,ne).muP    = muP;
                    runStruct(prInd,lrInd,ne).sigP   = sigP;
                    runStruct(prInd,lrInd,ne).wP     = wP;
                end
                
            % MoG
            otherwise
                runStruct(prInd,lrInd,ne).muPF   = muPF;
                runStruct(prInd,lrInd,ne).sigPF  = sigPF;
                runStruct(prInd,lrInd,ne).wPF    = wPF;
                
                if strcmp(priorF,'expon')
                    runStruct(prInd,lrInd,ne).avlog  = avlog;
                else
                    runStruct(prInd,lrInd,ne).muP    = muP;
                    runStruct(prInd,lrInd,ne).sigP   = sigP;
                    runStruct(prInd,lrInd,ne).wP     = wP;
                end
                
        end
        
        runStruct(prInd,lrInd,ne).gvlogF      = gvlogF;
        runStruct(prInd,lrInd,ne).gvlog       = gvlog;
        
        runStruct(prInd,lrInd,ne).hcF         = hcF;
        runStruct(prInd,lrInd,ne).hc          = hc;
        
        runStruct(prInd,lrInd,ne).nllF        = nllF;
        runStruct(prInd,lrInd,ne).nllGT       = nllGT(1);
        
        runStruct(prInd,lrInd,ne).dataset     = dataSets{ds};
        runStruct(prInd,lrInd,ne).filename    = parEstNames{ne};
        
        runStruct(prInd,lrInd,ne).numTrials   = numTrials;
        
        runStruct(prInd,lrInd,ne).numModelPar = size(initConds,2);
        
        temp                                  = cellfun(@(x) numel(x),ptvs_data);
        runStruct(prInd,lrInd,ne).numDataPt   = sum(temp(:));
        
        % Calculate BIC for this fit
        numPars                               = runStruct(prInd,lrInd,ne).numModelPar;
        numPts                                = runStruct(prInd,lrInd,ne).numDataPt;
        theseNll                              = runStruct(prInd,lrInd,ne).nllF;
        runStruct(prInd,lrInd,ne).bic         = getBIC(numPars,numPts,theseNll);
        
        runStruct(prInd,lrInd,ne).vStim1      = simData2.vStim1;
        runStruct(prInd,lrInd,ne).vStim2Delta = simData2.vStim2Delta;
        runStruct(prInd,lrInd,ne).cStim1      = simData2.cStim1;
        runStruct(prInd,lrInd,ne).cStim2      = simData2.cStim2;
    end
    
    if prInd < numPriors
        prInd    = prInd + 1;
    else
        prInd    = 1;
        lrInd   = lrInd + 1;
    end
    
end

save([workDir,'runStruct_',model,priorF,dStr],'runStruct');

end


%% Analyze piecewise estimation

if analysisPlot

if ~exist('runStruct')
    load([workDir,'runStruct_',model,priorF,dStr]);
end
    
% Lapse Rate
numLRs = 1;

% Num unique GT priors
numPriors = size(runStruct,1);

% Num samples per condition
numSamps  = size(runStruct(1,1,1).gvlogF,1);

% Unique pairings of samples & priors
sampPairs  = makeUniquePairs(numSamps);
priorPairs = makeUniquePairs(numPriors);

if numPriors > 5
    % Reorder so that sigma = 1 prior is the "canonical" one
%     rsTemp = runStruct;
%     runStruct(1,:,:) = rsTemp(2,:,:);
%     runStruct(2,:,:) = rsTemp(1,:,:);
    
    sigPairs  = [1 2 3 4 5];
    muPairs   = [3 6 7 8 9];
    
    allPairs  = {sigPairs,muPairs};
else
    allPairs = {[1:5]};
end

% Get matlab colors
[colormat] = colororder;

% Subplot positions
if numPriors > 5
    % Multiple dims, hard-code here since it's goofy
    spPos = [1 2 3 4 5 6 7 8 9];
    titleStr = {'\sigma = ','\sigma = ','\sigma = ','\sigma = ','\sigma = ',...
                '\mu = ','\mu = ','\mu = ','\mu = '};
    priorVars = [arrayfun(@(x) x.sigP(1),runStruct(1:5,1,1));...
                 arrayfun(@(x) x.muP(1),runStruct([1 6 7 8 9],1,1))];     
             
    for hh = 1:numel(allPairs)
        % Loop through first five, then last 4; goofy because titlestr is
        % accessed by looping through number of priors later on; should fix
        jsdStr{hh} = titleStr{1+5*(hh-1)};
        
        theseVals = priorVars([1:5] + 5*(hh-1));
        
        for ii = 1:numel(theseVals)
            jsdTix{hh,ii} = num2str(theseVals(ii));
        end
    end
else
    % Comparison along single dim
    spPos = 1:numPriors;
    switch priorF
        case 'expon'
            titleStr   = 'Log-slope = ';
            priorVars = arrayfun(@(x) x.avlog(1),runStruct(:,1,1));
        case 'mixgauss'
            titleStr   = '\sigma = ';
            priorVars = arrayfun(@(x) x.sigP(1),runStruct(:,1,1));
    end
    
    jsdStr{1} = titleStr;
    for ii = 1:numel(priorVars)
        jsdTix{ii} = num2str(priorVars(ii));
    end
    
end

% Plot a figure for each number of repeats, comparing fits
for i = 1:numTrCnts
    
    fitsFig{i} = figure;
    ftemp = fitsFig{i};
    if numPriors > 5
        ftemp.Position = [100 100+(i-1)*100 1770 660];
    else
        ftemp.Position = [100 100+(i-1)*100 1850 400];
    end
    hold on;
    
    clear theseFits thisGT
    
    for j = 1:numPriors
        
        switch model
            case 'ss'
                theseSlopes   = runStruct(j,numLRs,i).avlogF;
                
                % Recover best fit prior for each sample & plot with GT
                for k = 1:numSamps
                    
                    [thesePriors(k,:),theseVels(k,:)] = buildPiecewisePrior(theseSlopes(k,:));
                    
                end
                
                thisgtPrior = buildPiecewisePrior(runStruct(j,numLRs,i).avlog);
                
            otherwise
                % Get the fits
                theseMus   = runStruct(j,numLRs,i).muPF;
                theseSigs  = runStruct(j,numLRs,i).sigPF;
                theseWs    = runStruct(j,numLRs,i).wPF;
                
                % Get support vels that match piecewise prior
                [~,theseVels] = buildPiecewisePrior([]);
                
                % Recover best fit prior for each sample & plot with GT
                thesePriors = buildMoGPrior(theseSigs,theseMus,theseWs,theseVels);
                
                if strcmp(priorF,'expon')
                    thisgtPrior = buildPiecewisePrior(runStruct(j,numLRs,i).avlog);
                else
                    gtMu   = runStruct(j,numLRs,i).muP;
                    gtSig  = runStruct(j,numLRs,i).sigP;
                    gtW    = runStruct(j,numLRs,i).wP;
                    
                    thisgtPrior = buildMoGPrior(gtSig,gtMu,gtW,theseVels);
                end
        end
        
        % Collect sigmas
        sigArray{j} = runStruct(j,numLRs,i).hcF;
        sigGT{j}    = runStruct(j,numLRs,i).hc;
        
        % Get JSD measure of accuracy between fits and GT
        theseJSD(i,j,:) = getJSDiv(thesePriors,repmat(thisgtPrior,[numSamps 1]));
        
        % Build up distribution of JSD for all pairs in this set
        JSDdist{j}  = getJSDiv(thesePriors(sampPairs(:,1),:),thesePriors(sampPairs(:,2),:));
        
        runStruct(j,numLRs,i).JSDdist      = JSDdist{j};
        
        priors{j}   = thesePriors;
        
        % Get bias for each sample point of prior
        runStruct(j,numLRs,i).gtPrior   = thisgtPrior;
        runStruct(j,numLRs,i).fitPriors = thesePriors;
        runStruct(j,numLRs,i).velPrior  = theseVels(1,:);
        
        % Plot all the fits and the ground truth disribution
        if numPriors > 5
%             subplot(3,3,spPos(j));
            subplot(2,5,spPos(j));
        else
%             subplot(1,5,spPos(j));
            subplot(1,numPriors,spPos(j));
        end

        hold on;
        for kk = 1:numSamps
            plot(theseVels(1,:),thesePriors(kk,:),'color',colormat(i,:));
        end
        plot(theseVels(1,:),thisgtPrior,'k','linewidth',2);
        xlabel('Velocity');
        ylabel('p(V)');
        set(gca,'xtick',getLogXform(vStim1,0.3),'xticklabel',vStimLab,...
            'ylim',[1e-7 1e1],'xlim',getLogXform([0.5 12],0.3),'yscale',...
            'log','plotboxaspectratio',[1 1 1],'fontsize',15,'ytick',...
            [1e-6 1e-4 1e-2 1e0]);
        if numPriors > 5 %% hacky, but whatever.
            title([titleStr{j},num2str(priorVars(spPos(j)))]);
        else
            title([titleStr,num2str(priorVars(j))]);
        end
    end

    % Determine probability of detecting true difference between priors in
    % this set
    sigMat = zeros(numPriors,numPriors);
    pMat = zeros(numPriors,numPriors);
    
    numPairs = size(priorPairs,1);
    
    for k = 1:numPairs
    
        p1    = priors{priorPairs(k,1)};
        p2    = priors{priorPairs(k,2)};
     
        JSDpriorComp{k} = getJSDiv(p1(sampPairs(:,1),:),p2(sampPairs(:,2),:));
        
        % Percentage of estimates from different prior > same prior
        pDiff(k,i) = sum(JSDpriorComp{k} > JSDdist{priorPairs(k,1)})/size(sampPairs,1);
        sigDiff(k,i) = pDiff(k,i) > 0.95;
        
        sigMat(priorPairs(k,1),priorPairs(k,2)) = sigDiff(k,i);
        pMat(priorPairs(k,1),priorPairs(k,2)) = pDiff(k,i);
        
    end
    
    JSDpriorsTrCnts(i,1:numPairs) = JSDpriorComp;
    
    % Make Colormap
    gamma = linspace(0,1,256).^1;
    map = repmat(gamma',[1,3]);
    map = map.*repmat(colormat(i,:),[256 1]);
    
    % Create three mini comparison plots along mu/sigma/kurtosis dimensions
%     for kk = 1:numel(allPairs)
%         
%         thesePairs = allPairs{kk};
%         numThesePairs = numel(thesePairs);
%         
%         %%%% update indexing
%         pairFig{i,kk} = figure;
%         ftemp2 = pairFig{i,kk};
%         hold on;
%         xt = 1:numThesePairs;
%         yt = xt;
%         imagesc(xt,yt,pMat(thesePairs,thesePairs));
%         colormap(map);
%         set(gca,'xtickLabel',jsdTix(kk,:),'ytickLabel',jsdTix(kk,:),...
%                 'ytick',1:numThesePairs,'plotboxaspectratio',[1 1 1],...
%                 'fontsize',20,'xlim',[0.5 numThesePairs+0.5],...
%                 'ylim',[0.5 numThesePairs+0.5],'ydir','reverse');
%         xlabel(jsdStr{kk});
%         ylabel(jsdStr{kk});
%         title('p(JSD_{diff}>JSD_{same})');
%         
%         % label significant comparisons
%         [xm,ym] = meshgrid(xt,yt);
%         xms = xm(logical(sigMat(thesePairs,thesePairs)));
%         yms = ym(logical(sigMat(thesePairs,thesePairs)));
%         for kk = 1:numel(xms)
%             text(xms(kk),yms(kk),'*','fontsize',40,'fontweight','bold','horizontalalignment','center');
%         end
%         ftemp2.Position = [300 100+(i-1)*100 560 525];
%     end
end

save([workDir,'runStruct_',model,priorF,dStr],'runStruct','JSDpriorsTrCnts','priorPairs');

accFig = figure;
% offset = [-0.25 0 0.25];
offset = [0];
hold on;

for ii = 1:numTrCnts
    xpts = repmat([1:numPriors]',[1 numSamps]) + offset(ii);
    ypts = squeeze(theseJSD(ii,:,:));
    scatter(xpts(:),ypts(:),75,'linewidth',3);
end

tickLabs = cell(ii,1);
for ii = 1:numPriors
    if numPriors > 5
        tickLabs{ii} = [titleStr{ii},num2str(priorVars(spPos(ii)))];
        tickAng = 45;
        xlab = '';
    else
        tickLabs{ii} = [num2str(round(priorVars(spPos(ii)),1))];
        tickAng = 0;
        xlab = titleStr(1:end-3);
    end
end

set(gca,'xtick',1:numPriors,'xticklabels',tickLabs,'fontsize',20,'ylim',...
        [1e-5 1],'yscale','log','plotboxaspectratio',[1 1 1]);
xtickangle(gca,tickAng);
ylabel('JSD (fit vs ground truth)');
xlabel(xlab);
title('Fit Accuracy');
for ii = 1:numel(trCnts)
    legendEnts{ii} = ['Trial Reps = ',num2str(trCnts(ii))];
end
legend(legendEnts,'location','northwest');
accFig.Position = [700 200 575 535];

% Plot spread of likelihoods
likeFig = figure;
likeFig.Position = [500 400 750 550];
hold on;

conts = [0.05 0.075 0.1 0.2 0.4 0.5 0.8];

hcFlin = sigArray{j};
hcLin  = 2 - conts;
contsArr = repmat(conts,[size(hcFlin,1) 1]);

scatter(contsArr(:),hcFlin(:),75,'lineWidth',3);
plot(conts,hcLin,'--k','lineWidth',3)

xlabel('Contrast');
ylabel('Best fit hc');
set(gca,'xlim',[conts(1) conts(end)],'xtick',conts,'xscale','log','fontsize',20);

end


%% Save figures

if saveON
    
    if numel(allPairs) > 1
        pairSuffix = {'sig','mu','kurt'};
    else
        pairSuffix = {'allPairs'};
    end
    
    for i = 1:numTrCnts
        
        figID = ['N',num2str(trCnts(i))];
        figName1 = [figDir,model,priorF,figID,'_fits.svg'];
        
        saveas(fitsFig{i},figName1);
        
%         for j = 1:numel(allPairs)
%             figName2 = [figDir,model,priorF,figID,'_sigDiff_',pairSuffix{j},'.svg'];
%             saveas(pairFig{i,j},figName2);
%         end
    end
    
    figName3 = [figDir,model,priorF,'_acc.svg'];
    
    saveas(accFig,figName3);


end