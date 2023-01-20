function f1 = plotter(datamat,fitStr1,fitStr2)

v1s = unique(datamat(1,:));
v1sLog = getLogXform(v1s,0.3);
c1s = unique(datamat(2,:));
c2s = unique(datamat(4,:));

numV1 = numel(v1s);
numC1 = numel(c1s);
numC2 = numel(c2s);
numPfxns = numV1*numC2;

pars1{1} = fitStr1.gvlogF;
pars1{2} = fitStr1.hcF;
pars1{3} = fitStr1.avlogF;
% pFxn_pw = calcPFxn_piece(pars1,datamat);
[~,pFxn_pw] = cNLL_pw([fitStr1.hcF fitStr1.avlogF],datamat,v1s,c2s,fitStr1.gvlogF);

pars2{1} = fitStr2.gvlogF;
pars2{2} = fitStr2.hcF;
pars2{3} = fitStr2.wF;
pars2{4} = fitStr2.sigF;
pars2{5} = fitStr2.muF;
pFxn_mog = calcPFxn_mog(pars2,datamat);

ind = 1;

% Generate psychometric functions from MoG prob vector and collect subject
% response proportions
for ii = 1:numV1
    
    for kk = 1:numC2
        
        for jj = 1:numC1
            
            % Loop through different combos of stims
            pFxns(ii,jj,kk).v1 = v1s(ii);
            pFxns(ii,jj,kk).c1 = c1s(jj);
            pFxns(ii,jj,kk).c2 = c2s(kk);
            
            % Find all columns for this pFxn
            theseInds = datamat(1,:) == v1s(ii) & ...
                        datamat(2,:) == c1s(jj) & ...
                        datamat(4,:) == c2s(kk);
            
            % Grab the set of responses
            theseResps        = datamat(9,theseInds);
            
            % Grab test stimulus velocities and indices
            [theseVels,tx,tv] = unique(datamat(3,theseInds));
            
            pResp = nan(1,numel(theseVels));
            pRespTr = nan(1,numel(theseVels));
            
            % Calculate proportion of responses
            for hh = 1:numel(theseVels)
                vInds = tv == hh;
                pResp(hh)   = sum(theseResps(vInds))/sum(vInds);
                pRespTr(hh) = sum(vInds);
            end
            
            pFxns(ii,jj,kk).v2      = theseVels;
            pFxns(ii,jj,kk).pResp   = pResp;
            pFxns(ii,jj,kk).pRespTr = pRespTr;
            
            % Collect all probability fits for this set of trials
            theseMog         = pFxn_mog(theseInds);
            % Select only unique instances (should be >=1 copies)
            ptvsF2{ii,jj,kk} = theseMog(tx);
            
            % Piecewise
            thesePW          = pFxn_pw(theseInds);
            % Select only unique instances (should be >=1 copies)
            ptvsF1{ii,jj,kk} = thesePW(tx);
            
        end
        
        pFxnPos(ind,:) = [ii,kk];
        ind = ind + 1;
        
    end
    
end

% Permute velocities since calcPtvsSSData expects it in c1,v1,c2 order,
% then permute psych function vals back to v1,c1,c2 order
% a = arrayfun(@(x) x.v2,pFxns,'uniformoutput',false);
% a = permute(a,[2 1 3]);
% ptvsF1 = calcPtvsSSData(c1s',v1s',c2s',a,...
%         fitStr1.avlogF,fitStr1.gvlogF,fitStr1.hcF,1);
% ptvsF1 = permute(ptvsF1,[2 1 3]); 

% keyboard

c1 = 1;

refContInds = [2 6];

f1 = figure;
f1.Position = [100 100 1266 1175];

for ii = 1:numPfxns
    
    thisVel    = pFxnPos(ii,1);
    thisCont   = pFxnPos(ii,2);
    theseVels  = pFxns(thisVel,c1,thisCont).v2;
    theseProps = pFxns(thisVel,c1,thisCont).pResp;
    thesePpsTr = pFxns(thisVel,c1,thisCont).pRespTr;
    
    if ~isempty(theseVels)
    subplot(numV1,numC2,ii);
    hold on
    
    % Calculate PSEs
    piecePSE = v1sLog(thisVel) + fitStr1.avlogF(thisVel)* ...
                              ( (fitStr1.gvlogF(thisVel)*fitStr1.hcF(refContInds(c1)))^2 - ...
                                (fitStr1.gvlogF(thisVel)*fitStr1.hcF(thisCont))^2  );
    if numel(fitStr2.sigF) == 1
        mogPSE   = v1sLog(thisVel) * ((fitStr2.gvlogF(thisVel)*fitStr2.hcF(thisCont))^2 + fitStr2.sigF^2)/...
                                     ((fitStr2.gvlogF(thisVel)*fitStr2.hcF(refContInds(c1)))^2 + fitStr2.sigF^2);
    else
        %%% no closed form solution for PSE with # comps > 1... 
        mogPSE = interp1(theseVels,ptvsF2{thisVel,c1,thisCont},0.5);
    end
        
    piecePSE = getLinXform(piecePSE,0.3);
    
    mogPSE   = getLinXform(mogPSE,0.3);                         
    
    % Fit weibull function
    [thisDataFit,dfPSE] = fitCumGauss(theseVels,theseProps,thesePpsTr);
                          
    scatter(theseVels,theseProps,50,'k');
    plot(theseVels,ptvsF1{thisVel,c1,thisCont},'r','linewidth',2);
    plot(piecePSE*[1 1],[0 0.5],'r','linewidth',1.5);
    plot(theseVels,ptvsF2{thisVel,c1,thisCont},'b','linewidth',2);
    plot(mogPSE*[1 1],[0 0.5],'b','linewidth',1.5);
    plot(theseVels,thisDataFit,'g','linewidth',2);
    plot(dfPSE*[1 1],[0 0.5],'g','linewidth',1.5);
    
    plot(v1s(thisVel)*[1 1],[0 1],'--k');
    
    set(gca,'fontsize',10,'ylim',[0 1],'xtick',v1s,'xlim',[theseVels(1) theseVels(end)]);
    if thisVel == numV1
    xlabel('Test velocity (\circ/s)');
    end
    if thisCont == 1
    ylabel('p("v_{2}">"v_{1}")');
    end
    else
        continue
    end
    
end

end

function [pFxn] = calcPFxn_piece(pars,dataMat)

% Calculate psych function values given a set of observer/exp parameters

% Extract pars
gvlog = pars{1};
hc    = pars{2};
avlog = pars{3};

vStim1 = dataMat(1,:);
cStim1 = dataMat(2,:);
vStim2 = dataMat(3,:);
cStim2 = dataMat(4,:);

% Get logxforms
rf_vlog   = log(1 + vStim1/0.3);
test_vlog = log(1 + vStim2/0.3);

% Index into contrast pars
[~,~,tmpCntVec] = unique([cStim1;cStim2]);
tmpCntMat       = reshape(tmpCntVec,[2 size(dataMat,2)]);
tsCont          = tmpCntMat(1,:);
rfCont          = tmpCntMat(2,:);

% Index into velocity pars; use nearest neighbor Ref vel to index Test vel
[uniqVels,~,rfVel] = unique(vStim1);
% distVals  = abs(repmat(uniqVels',[1 size(dataMat,2)]) - ...
%                 repmat(vStim2,[numel(uniqVels) 1]));
% [~,tsVel] = min(distVals);

% Just make test index equal reference like scratch script
tsVel = rfVel;
    
% Grab likelihood sigmas from exp parameters
thisRSigL    = gvlog(rfVel).*hc(rfCont);
thisTSigL    = gvlog(tsVel).*hc(tsCont);

% Use nearest neighbor Ref vel to index into fitted slopes
refSlope  = avlog(rfVel);
testSlope = avlog(tsVel);

% Calculate p(v2_hat>v1_hat) with piecewise method (using linear
% interpolation)
rfMu    = rf_vlog + refSlope.*thisRSigL.^2;
rfVar   = thisRSigL.^2;

testMu  = test_vlog + testSlope.*thisTSigL.^2;
testVar = thisTSigL.^2;

Z    = (testMu - rfMu) ./ sqrt(rfVar + testVar);

% pFxn = min(max(normcdf(Z), eps), 1 - eps);
pFxn = normcdf(Z);

end

function [pFxn] = calcPFxn_mog(pars,dataMat)

% Calculate psych function values given a set of observer/exp parameters

% Assumes dataMat is a mxn matrix where first 4 rows hold stimulus values
% and n is number of trials

numTrials = size(dataMat,2);

% Extract pars (assuming all components free version)
gvlog = pars{1};
hc    = pars{2};
wP    = pars{3};
sigP  = pars{4};
muP   = pars{5};

vStim1 = dataMat(1,:);
cStim1 = dataMat(2,:);
vStim2 = dataMat(3,:);
cStim2 = dataMat(4,:);

% Get logxforms
rf_vlog   = log(1 + vStim1/0.3);
test_vlog = log(1 + vStim2/0.3);

% Assign index to each unique contrast
[~,~,tmpCntVec] = unique([cStim1;cStim2]);
tmpCntMat       = reshape(tmpCntVec,[2 size(dataMat,2)]);
rfCont          = tmpCntMat(1,:);
tsCont          = tmpCntMat(2,:);

% Assign index to each unique velocity; use nearest neighbor Ref vel to 
% index Test vel
[uniqVels,~,rfVel] = unique(vStim1);
distVals  = abs(repmat(uniqVels',[1 size(dataMat,2)]) - ...
                repmat(vStim2,[numel(uniqVels) 1]));
[~,tsVel] = min(distVals);

% Grab likelihood sigmas from exp parameters
thisRSigL    = gvlog(rfVel).*hc(rfCont);
thisTSigL    = gvlog(tsVel).*hc(tsCont);

% Calculate p(v2_hat>v1_hat) with MoG estimation
% MoG post pars are (nxm) where n = # components, m = # likelihoods
[alphas1,muTildes1,wTildes1] = getMogPostPars(wP,sigP,muP,thisRSigL,rf_vlog);
[alphas2,muTildes2,wTildes2] = getMogPostPars(wP,sigP,muP,thisTSigL,test_vlog);

% Loop over trials (or figure out how to reorg pars, since you'd need a 3D
% mat)
pFxn = nan(1,numTrials);

for ii = 1:numTrials
    denom   = sqrt((thisTSigL(ii)^2)*(alphas2(:,ii).^2) + (thisRSigL(ii)^2)*(alphas1(:,ii).^2));
    
    Z       = (muTildes2(:,ii)  + alphas2(:,ii)*test_vlog(ii) ...
             - muTildes1(:,ii)' - alphas1(:,ii)'*rf_vlog(ii)   )./denom;

    pFxn(ii)    = wTildes2(:,ii)'*normcdf(Z)*wTildes1(:,ii);
end

end

function [nll,rho] = cNLL_pw(parVec,subjData,vStim1,cStim2,GTgvlog)

%% Likelihood of data given current model pars

% avlog = parVec(1:numel(vStim1));
% hc    = parVec(numel(vStim1) + 1:end);
avlog = parVec(numel(cStim2) + 1:end);
hc    = parVec(1:numel(cStim2));
gvlog = GTgvlog;

% moving test slope outside of for loop for speed when using
% interpolation. Actually, nothing else
% here really needs a for loop, just used it for clarity when debugging

% Define test velocities
testVelLin    = subjData(3,:);

% Convert to normalized log space
testVelLog  = log(1 + testVelLin/0.3);

% for each trial, calculate mus and vars of map sampling for ref and
% test
for x = 1:size(subjData,2)
    
    %% reference
    rfCont = subjData(2,x);
    
    % Convert to normalized vlog space
    rfVelLin(x) = subjData(1,x);
    rfVelLog(x)  = log(1 + rfVelLin(x)/0.3);
    
    % mean and var of standard posterior in tranformed velocity space
    % hc = all contrasts; assume numel(cStim2) == numel(hc) and
    % cStim1 is a subset of cStim2
    refSlope(x) = avlog(rfVelLin(x) == vStim1);
    
    refMu(x)   = rfVelLog(x) + refSlope(x)*(gvlog(rfVelLin(x) == vStim1)*hc(rfCont == cStim2))^2;
    refVar(x)  = (gvlog(rfVelLin(x) == vStim1)*hc(rfCont == cStim2))^2;
    
    %% test
    testCont = subjData(4,x);
    testSlope(x) = avlog(rfVelLin(x) == vStim1);
    
    testMu(x)      = testVelLog(x) + testSlope(x)*(gvlog(rfVelLin(x) == vStim1)*hc(testCont == cStim2))^2;
    testVar(x)     = (gvlog(rfVelLin(x) == vStim1)*hc(testCont == cStim2))^2;
    
end

% naeker way
Z = (testMu - refMu) ./ sqrt(refVar + testVar);

rho = normcdf(Z);

% Compute log-likelihood - log of bernoulli distribution
nll = -sum(subjData(9,:) .* log(rho) + ...
    (1 - subjData(9,:)) .* log(1 - rho));

end