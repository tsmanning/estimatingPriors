clear all
close all

%% set up parameters

% Normal function (support, mean, std dev)
gF    = @(x,mu,sig) normpdf(x,mu,sig);

% support
x = -20:0.01:20;

% number of measurements to run
num_trials = 7500;

% Define grid of observer parameters
gam1 = linspace(0.5,1.5,10);
gam2 = linspace(0.5,1.5,10);
w1   = linspace(0.1,0.5,5);
sig  = linspace(0.25,1,8);
stim = [0.5 1 2 4 8 12];

% gam1 = 0.5;
% gam2 = 0.5;
% w1   = 0.5;
% sig  = 0.3;
% stim = [0.5 1 2 4 8 12];

% Define unique conditions to bootstrap
[allConds1] = makeCombos([numel(gam1) numel(gam2) numel(w1) numel(sig)]);

% Cull repeated observers
endInds = repmat([numel(gam1):numel(gam1):numel(gam1)*numel(gam2)]',[numel(w1)*numel(sig) 1]);
begInds = repmat([1:numel(gam1)]' + [0:(numel(gam1)-1)]'*numel(gam1),[numel(w1)*numel(sig) 1]);

mult = repelem(numel(gam1)*numel(gam2)*[0:(numel(w1)*numel(sig)-1)],numel(gam1))';

bI = begInds + mult;
eI = endInds + mult;

allConds = [];

for ii = 1:numel(bI)
    allConds = [allConds;allConds1(bI(ii):eI(ii),:)];
end

% Set num hist bins
numBins = 80;

% Save dir
splPath = regexp(which('errorBootstrap'),filesep,'split');
fDir    = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];
sDir    = [fDir,'SimData',filesep];

% date string ID
ds      = datestr(now,30);

%% run bootstrapping

% hardcode stimulus value for now
mInd = 2;

% initialize 
numerSampDist = nan(size(allConds,1),numBins);
newSampDist   = nan(size(allConds,1),numBins);
oldSampDist   = nan(size(allConds,1),numBins);

% Loop over unique observer parameters
if 1 
for ii = 1:size(allConds,1)
    
    if mod(ii,10) == 0
       disp(['running bootstrap ',num2str(ii),' of ',num2str(size(allConds,1))]); 
    end
    
    thisgam1 = gam1(allConds(ii,1));
    thisgam2 = gam2(allConds(ii,2));
    thisw1   = w1(allConds(ii,3));
    thisw2   = 1 - thisw1;
    thissig  = sig(allConds(ii,4));
    
    % generate likelihood for this set of runs
    ms = thissig*randn(1,num_trials) + stim(mInd);
    
    % Get tilde parameters using expected value of likelihood
    [aExp,mExp,wTExp] = getMogPostPars([thisw1 thisw2],[thisgam1 thisgam2],[0 0],thissig,stim(mInd));
    
    for t = 1:num_trials
        
        meas = ms(t);
    
%         [~,~,wT,xHat(t)] = getMogPostPars([thisw1 thisw2],[thisgam1 thisgam2],[0 0],thissig,meas);
        [~,~,wT,~] = getMogPostPars([thisw1 thisw2],[thisgam1 thisgam2],[0 0],thissig,meas);
        
        mu = [0 0];
        g = [thisgam1 thisgam2];
        
        % numerically compute likelihood and prior
        lik   = gF(x,meas,thissig);
        
        prior = zeros(1,numel(x));
        
        for kk = 1:numel(wT)
            prior = prior + wT(kk)*normpdf(x,mu(kk),g(kk));
        end
        
        %% numerically compute posterior
        post    = lik.*prior;                               % sums to ~10.4 and should match TM derivation (no normalization)
        
        xHatNum(t) = sum(post.*x)/sum(post);
        
    end
    
    % Define edge starts/ends based on +/-3SD mean of estimate distribution                    
    hstart = (wTExp(1)*(aExp(1)*stim(mInd) + mExp(1)) + ...
              wTExp(2)*(aExp(2)*stim(mInd) + mExp(2))) - ...
              3*(sqrt( (thissig^2)*(wTExp(1)*aExp(1) + wTExp(2)*aExp(2))^2 ));
    hend   = (wTExp(1)*(aExp(1)*stim(mInd) + mExp(1)) + ...
              wTExp(2)*(aExp(2)*stim(mInd) + mExp(2))) + ...
              3*(sqrt( (thissig^2)*(wTExp(1)*aExp(1) + wTExp(2)*aExp(2))^2 ));
    
    % Define edges to use for histograms
    edges = linspace(hstart,hend,numBins+1);
    centers = edges(1:end-1) + diff(edges(1:2))/2;
    
    % histogram xHat
    numerSampDist(ii,:) = histcounts(xHatNum,edges,'normalization','pdf');
    
    % New estimate distribution
    newSampDist(ii,:)  = gF(centers,wTExp(1)*(aExp(1)*stim(mInd) + mExp(1)) + ...
                                    wTExp(2)*(aExp(2)*stim(mInd) + mExp(2)),...
                            sqrt( (thissig^2)*(wTExp(1)*aExp(1) + wTExp(2)*aExp(2))^2 ));
    % Old estimate distribution
    oldSampDist(ii,:)  = gF(centers,wTExp(1)*(aExp(1)*stim(mInd) + mExp(1)) + ...
                                    wTExp(2)*(aExp(2)*stim(mInd) + mExp(2)),...
                            sqrt( (thissig^2)*(wTExp(1)*aExp(1)^2 + wTExp(2)*aExp(2)^2) ));
                        
    edgesOut(ii,:)    = edges;
    centersOut(ii,:)  = centers;
    expectedVal(ii,1) = wTExp(1)*(aExp(1)*stim(mInd) + mExp(1)) + ...
                        wTExp(2)*(aExp(2)*stim(mInd) + mExp(2));
    newSig(ii,1)      = sqrt( (thissig^2)*(wTExp(1)*aExp(1) + wTExp(2)*aExp(2))^2 );                   
    oldSig(ii,1)      = sqrt( (thissig^2)*(wTExp(1)*aExp(1)^2 + wTExp(2)*aExp(2)^2) );      
    
end

errNew = (newSampDist-numerSampDist)./numerSampDist;
errOld = (oldSampDist-numerSampDist)./numerSampDist;

save([sDir,'approxBootstraps_',ds],'numerSampDist','newSampDist','oldSampDist','errNew','errOld',...
                                   'edgesOut','centersOut','expectedVal','newSig','oldSig');
else
    
end
    
%% Compare numeric vs analytic approximations
figure;
histogram(newSampDist-oldSampDist,'edgecolor','none');
set(gca,'fontsize',20);
xlabel('New approx - old approx');
ylabel('Count');

figure;
hold on;
histogram(errNew,'edgecolor','none');
histogram(errOld,'edgecolor','none');
set(gca,'fontsize',20,'xlim',[-1 1]);
xlabel('Analytic - numerical (fractional error)');
ylabel('Count');
legend({'New approx','Old approx.'});

% Reorganize
numerSD = reshape(numerSampDist,[numel(gam1) numel(gam2) numel(w1) numel(sig) numBins]);
newSD = reshape(newSampDist,[numel(gam1) numel(gam2) numel(w1) numel(sig) numBins]);


% figure;
% hold on;
% plot(numerSampDist,'k','linewidth',2);
% plot(newSampDist,'--b','linewidth',2);
% plot(oldSampDist,'--r','linewidth',2);
% 
% figure;
% hold on;
% plot(errNew,'--b','linewidth',2);
% plot(errOld,'--r','linewidth',2);






