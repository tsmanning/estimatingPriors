% Fitting analysis for log-gaussian prior in stocker-simoncelli framework


%% Select analysis to do

loadIn      = 1;
plotPars    = 0;
plotLRPrec  = 0;
pSDTcomp    = 0;
biasPlot    = 1;

%% Select parameters to plot

numRunsInd = 3;
% [1 - n=20, 2 - n=80, 3 - n=500, 4 - n=5000]

pSelect = 1;
% [1 - a(v)"slope", 2 - g(v), 3 - h(c)]

pNames = {'avlog','gvlog','hc'};
param = pNames{pSelect};
parFit = [param,'F'];


%% Define Datasets

workDir = which('LRComp');
workDir = [workDir(1:end-numel('LRComp.m')),'ssGauss/'];

% dataSets = {'gt_gauss_0.1_LR_0','gt_gauss_0.2_LR_0','gt_gauss_0.4_LR_0','gt_gauss_0.8_LR_0',...
%             'gt_gauss_0.1_LR_0.05','gt_gauss_0.2_LR_0.05','gt_gauss_0.4_LR_0.05','gt_gauss_0.8_LR_0.05'};
dataSets = {'gt_log_gauss_0.8_LR_0',   'gt_log_gauss_1_LR_0',   'gt_log_gauss_1.2_LR_0',   'gt_log_gauss_1.4_LR_0',...
            'gt_log_gauss_0.8_LR_0.05','gt_log_gauss_1_LR_0.05','gt_log_gauss_1.2_LR_0.05','gt_log_gauss_1.4_LR_0.05'};
% dataSets = {'gt_log_gauss_1.4_LR_0'};

% parEstNames = {'parEsts_n20_collect','parEsts_n80_collect',...
%                'parEsts_n250_collect','parEsts_n_jitter_5000'};
parEstNames = {'parEsts_ss_jitter_n20','parEsts_ss_jitter_n80','parEsts_ss_jitter_n250'};
        
trCnts      = [20 80 250];
lapseRates  = [0 0.05];
sigmas      = [0.8 1 1.2 1.4];
% trCnts      = [80];
% lapseRates  = [0];
% sigmas      = [1.4];
           
numDatasets = numel(dataSets);
numEstRuns  = numel(parEstNames);
numLRs      = numel(lapseRates);
numSigmas   = numel(sigmas);

minlikeli   = 1e-12;


%% Load it all in, put it in an index struct

if loadIn
    
aInd    = 1;
lrInd   = 1;

% Loop over simulated datasets
for i = 1:numDatasets
    
    dsPath = [workDir,dataSets{i},filesep,'simData'];
    load(dsPath);   
    
    % Loop over estimate runs
    for j = 1:numEstRuns
        
        thisEstRun = parEstNames{j};        
        
        estPath = [workDir,dataSets{i},filesep,thisEstRun];
        load(estPath);
        
        % dim 1: slope
        % dim 2: lapse rate
        % dim 3: sample size
        
        runStruct(aInd,lrInd,j).avlogF      = avlogF;
        runStruct(aInd,lrInd,j).avlog       = avlog;
        
        runStruct(aInd,lrInd,j).gvlogF      = gvlogF;
        runStruct(aInd,lrInd,j).gvlog       = gvlog;
        
        runStruct(aInd,lrInd,j).hcF         = hcF;
        runStruct(aInd,lrInd,j).hc          = hc;
        
        runStruct(aInd,lrInd,j).nllF        = nllF;
        
        runStruct(aInd,lrInd,j).dataset     = dataSets{i};
        runStruct(aInd,lrInd,j).filename    = thisEstRun;
        
        runStruct(aInd,lrInd,j).numTrials   = numTrials;
    end
    
    dsStruct(i).vStim1      = vStim1;
    dsStruct(i).vStim2Delta = vStim2Delta;
    dsStruct(i).cStim1      = cStim1;
    dsStruct(i).cStim2      = cStim2;
    
    if aInd < numSigmas
        aInd    = aInd + 1;
    else
        aInd    = 1;
        lrInd   = lrInd + 1;
    end
    
end

end


%% Compare precision & accuracy of parameter estimates

dNorm = @(x,sigma,mu) -x.*(1/(sqrt(2*pi)*sigma^3)).*exp(-0.5*(x/sigma).^2);

if plotPars

% probably institute some sort of outlier detection and culling?

f1 = figure;
f1.Position = [100 100 630 630];
hold on;

% Set subplot index
spInd = 1;

for i = 1:numLRs
    for j = 1:numSigmas
        
        subplot(numLRs,numSigmas,spInd);
        hold on;
        
        fitNLLs   = runStruct(j,i,numRunsInd).nllF;
        theseFits   = runStruct(j,i,numRunsInd).(parFit);
        thisGT      = dNorm(getLogXform(vStim1,0.3),sigmas(j),1e-14);
        
        if pSelect < 3
            thisDV      = vStim1;
        else
            thisDV      = cStim2;
        end
        
        numRuns = size(theseFits,1);
        [~,bestFit] = min(fitNLLs);

        % Plot ground truth values
        plot(thisDV,thisGT,'--k');
        
        % Plot parameter estimates from each run
%         for h = 1:numRuns
            
            scatter(thisDV,theseFits(bestFit,:),[],[0.5 0.5 1],'filled');
            
%         end

%         % Plot median PEs
%         scatter(thisDV,median(theseFits,1),[],'k','filled');
        
        % Label yer axes
        
        if pSelect == 1
            title(['LR = ',num2str(lapseRates(i)),'; ','\sigma',' = ',num2str(sigmas(j))]);
            set(gca,'ylim',[-1.5 1.5],'xlim',[vStim1(1) vStim1(end)],'xtick',vStim1,'xscale','log');
            xlabel('V_{standard}');
        elseif pSelect == 2
            title(['LR = ',num2str(lapseRates(i)),'; ',param(1),' = ',num2str(sigmas(j))]);
            set(gca,'ylim',[-2 2],'xlim',[vStim1(1) vStim1(end)],'xtick',vStim1,'xscale','log');
            xlabel('V_{standard}');
        elseif pSelect == 3
            title(['LR = ',num2str(lapseRates(i)),'; ',param(1),' = ',num2str(sigmas(j))]);
            set(gca,'ylim',[-2 4],'xlim',[cStim2(1) cStim2(end)],'xtick',cStim2,'xscale','log');
            xlabel('Contrast');
            xtickangle(45);
            f1.Position = [100 100 1200 900];
        end
        
        ylabel('Best Fit \sigma');

        
        spInd = spInd + 1;
        
    end
end

end


%% Compare precision via slope estimates & likelihood ratios

if plotLRPrec
    
f2 = figure;
f2.Position = [400 100 1075 680];
hold on;

numtrials = trCnts(numRunsInd);

numVels = numel(vStim1);

% Set subplot index
spInd = 1;

for i = 1:numLRs
    for j = 1:numSigmas
        
        disp([num2str(i),'/',num2str(numLRs),'; ',num2str(j),'/',num2str(numSigmas)]);
        
        subplot(numLRs,numSigmas,spInd);
        hold on;
        
        % Grab parameters for this run
        aP = runStruct(j,i,numRunsInd).avlogF;
        gP = runStruct(j,i,numRunsInd).gvlogF;
        cP = runStruct(j,i,numRunsInd).hcF;
        thisSigma = sigmas(j);
        testSigs = thisSigma*[0.25 0.5 1.5 2];
        numTestSigs = numel(testSigs);
        
        % number of fitting runs
        numRuns = numel(ptvsF);
        
        for h = 1:numRuns
            % THIS SHOULD BE THE SAVED DATA VALUES, NOT THE BEST FIT FOR 
            % EACH RUN
            ptvs_meas = ptvs_data;
            
            % Make new set of parameters varying slope from best fit value
            parVec_temp = [gP(h,:) cP(h,:)];
            xvals = getLogXform(vStim1,0.3);
            slopeMat = [dNorm(xvals,testSigs(1));
                        dNorm(xvals,testSigs(2));
                        dNorm(xvals,testSigs(3));
                        dNorm(xvals,testSigs(4))];
%             slopeMat = repmat(sigmas',[1 numVels]); % just use other GT vals?
            parVec = [slopeMat repmat(parVec_temp,[numTestSigs 1])];
            
            % Add actual best fit pars to bottom 
            % (parVec dims: numSlopes + 1 x numParams)
            parVec = [parVec;[dNorm(xvals,sigmas(j)) gP(h,:) cP(h,:)]];

            for k = 1:numTestSigs+1
                nll(h,k) = calculate_model_nll(parVec(k,:),ptvs_meas,vStim1,...
                                        cStim1,vStim2Delta,cStim2,numtrials);
            end
        end
        
        % Get likelihood ratio for range of sigmas (assuming constant over 
        % velocity range), including best fit slope
        
        likeRats = -2*(-nll + repmat(nll(:,end),[1 numTestSigs + 1]));
        
        % Plot Values
        for h = 1:numRuns
            
            [sortSlopes,sI] = sort([testSigs thisSigma]);
            
            plot(sortSlopes,likeRats(h,sI),'o-k');
%             scatter(mean(aP(h,:)),likeRats(h,end),[],[1 0.5 0.5],'filled');
            plot([thisSigma thisSigma],[-1.1 2.1],'--b'); 
        end
        
        % Label yer axes
        title(['LR = ',num2str(lapseRates(i)),'; a = ',num2str(sigmas(j))]);
        xlabel('Input Sigma');
        ylabel('Likelihood Ratio');
        set(gca,'xtick',sortSlopes,'xlim',[0 sortSlopes(end)],'ylim',[-1 1]);
        xtickangle(45);
        
        spInd = spInd + 1;
        
    end
end

end

%% Compare psychometric function fits

if pSDTcomp
    
    f3 = figure;
    f3.Position = [100 100 650 650];
    hold on;
    
    for i = 1:numLRs
        for j = 1:numSigmas
            
            tempStruc = runStruct(j,1,numRunsInd);
            
            % Define which condition to plot across different observers
            vStimInd = 3;
            c1Ind = 1;
            c2Ind = 2;
            
            vStim1 = dsStruct(i).vStim1(vStimInd);
            cStim1 = dsStruct(i).cStim1;
            vStim2Delta = dsStruct(i).vStim2Delta;
            cStim2 = dsStruct(i).cStim2;
            numtrials = tempStruc.numTrials;
            
            vStim2 = vStim1.*vStim2Delta;
            
            % Get ground truth psychometric function
            
            avlog = tempStruc.avlog;
            gvlog = tempStruc.gvlog;
            hc = tempStruc.hc;  
            
            [ptvsGT,~] = ...
                calculate_ptvs(vStim1,cStim1,vStim2Delta,cStim2,numtrials,avlog,gvlog,hc); 
            
            % Plot em
            
            subplot(numLRs,numSigmas,j + numSigmas*(i-1));
            hold on;
            
            plot(vStim2,squeeze(ptvsGT(c1Ind,1,c2Ind,:)),'k','linewidth',2);
            
            numRuns = numel(ptvsF);
            
            for h = 1:numRuns
                
                thisRun = ptvsF{h};
                
                plot(vStim2,squeeze(thisRun(c1Ind,vStimInd,c2Ind,:)),'--b');
                
            end
            
            xlabel('V_{test}');
            ylabel('p(v_{test} > v_{standard})');
            title(['LR = ',num2str(lapseRates(i)),'; a = ',num2str(sigmas(j))]);
            set(gca,'ylim',[0 1],'xlim',[vStim2(1) vStim2(end)],'xtick',round(vStim2,2));
            xtickangle(45);
        end
    end
    
    
end

%% Compare bias vs slope

if biasPlot
    
    f4 = figure;
    f4.Position = [100 100 1100 950];
    hold on;
    
    % Should get this from numrows
    numResamps = 5;
    
    numS1Vels = numel(vStim1);
    numS2Conts = numel(cStim2);
        
    % Index of desired contrast level into unique contrast vector - 
    % should find better way to do this 
    s1C = 6; % 0.5
    
    % Lapse rate to test
    lr2t = 1;
    vs2D = vStim2Delta;
    % dim 1: slope
    % dim 2: lapse rate
    % dim 3: sample size
    
    % Loop over vels
    for i = 1:numS1Vels
        % Loop over conts
        for j = 1:numS2Conts

            disp([num2str(i),'/',num2str(numS1Vels),'; ',num2str(j),'/',num2str(numS2Conts)]);
            
            subplot(numS1Vels,numS2Conts,j + numS2Conts*(i-1));
            hold on;
            
            % Loop over observer slopes
            for sl = 1:numSigmas
                
                avlogF = runStruct(sl,lr2t,numRunsInd).avlogF;
                gvlogF = runStruct(sl,lr2t,numRunsInd).gvlogF;
                hcF    = runStruct(sl,lr2t,numRunsInd).hcF;
                
                % Loop over samples
                for k = 1:numResamps
                    
                    % Find bias for each run
                    [~,~,biasF(k)] = ...
                        findPSENum(vStim1,vs2D,cStim2,avlogF(k,:),gvlogF(k,:),hcF(k,:),...
                        [s1C j],i,0);
                    
                end
                
                % Get mean slope (or should this just use locally
                % estimated slope??)
                slope = mean(avlogF,2);
                    
                scatter(slope,biasF,[],[1-0.75*(sl/numSigmas) 0 0],'filled');

            end
            
            xlabel('Best Fit Prior Slope');
            ylabel('Bias (PSE_{fit}:True)');
            title(['V_{s}:',num2str(vStim1(i)),...
                   '; c_{t}:',num2str(cStim2(j))]);
            set(gca,'ylim',[0.5 1.5],'xlim',[-4 1]);        
                                
        end
    end
    
end

