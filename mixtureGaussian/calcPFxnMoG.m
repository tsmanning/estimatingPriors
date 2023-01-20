function [ptvs] = calcPFxnMoG(vs1,cs1,vStim2Delta,cs2,gvlog,hc,sigP,muP,wP,method,estDistEval)

% Data wrangling
if ~iscolumn(vs1)
    vs1 = vs1';
end

if ~iscolumn(vStim2Delta)
    vStim2Delta = vStim2Delta';
end

if ~isrow(gvlog)
    gvlog = gvlog';
end

if ~isrow(hc)
    hc = hc';
end

% Make vectors of trial parameters
% vs2 = vs1*vStim2Delta;

% combs = makeCombos([numel(vs1) numel(cs1) numel(vs2) numel(cs2)]);
combs = makeCombos([numel(vs1) numel(cs1) numel(vStim2Delta) numel(cs2)]);
combs = combs';

vStim1 = vs1(combs(1,:)); 
cStim1 = cs1(combs(2,:)); 
vStim2 = vStim2Delta(combs(3,:)).*vs1(combs(1,:)); 
cStim2 = cs2(combs(4,:));

numTrials = size(combs,2);

% Get logxforms
rf_vlog   = log(1 + vStim1/0.3);
test_vlog = log(1 + vStim2/0.3);

% Assign index to each unique contrast
[~,~,tmpCntVec] = unique([cStim1;cStim2]);
tmpCntMat       = reshape(tmpCntVec,[2 numTrials]);
rfCont          = tmpCntMat(1,:);
tsCont          = tmpCntMat(2,:);

% Assign index to each unique velocity; use nearest neighbor Ref vel to 
% index Test vel
[uniqVels,~,rfVel(1,:)] = unique(vStim1);
distVals  = abs(repmat(uniqVels,[1 numTrials]) - ...
                repmat(vStim2',[numel(uniqVels) 1]));
% [~,tsVel] = min(distVals);
tsVel = 1;

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
    switch method
        case 'numeric'
            
            error('Not yet implemented');
            
        case 'closed'
            switch estDistEval
                case 'original'
                    denom   = sqrt((thisTSigL(ii)^2)*(alphas2(:,ii).^2) + (thisRSigL(ii)^2)*(alphas1(:,ii).^2));
                    
                    Z       = (muTildes2(:,ii)  + alphas2(:,ii)*test_vlog(ii) ...
                        - muTildes1(:,ii)' - alphas1(:,ii)'*rf_vlog(ii)   )./denom;
                    
                    pFxn(ii) = wTildes2(:,ii)'*normcdf(Z)*wTildes1(:,ii);
                    
                case 'rederived'
                    denom   = sqrt((thisTSigL(ii).^2).*sum(wTildes2(:,ii).*alphas2(:,ii)).^2 + ...
                        (thisRSigL(ii).^2).*sum(wTildes1(:,ii).*alphas1(:,ii)).^2);
                    
                    Z       = ( sum(wTildes2(:,ii).*(muTildes2(:,ii) + alphas2(:,ii).*test_vlog(ii))) - ...
                        sum(wTildes1(:,ii).*(muTildes1(:,ii) + alphas1(:,ii).*rf_vlog(ii)  )) )./denom;
                    
                    %             db.denom(ii) = denom;
                    %             db.Z(ii)     = Z;
                    %             db.thisTSigL(ii) = thisTSigL(ii);
                    %             db.thisRSigL(ii) = thisRSigL(ii);
                    %             db.alphas2(:,ii)   = alphas2(:,ii);
                    %             db.alphas1(:,ii)   = alphas1(:,ii);
                    %             db.muTildes2(:,ii)   = muTildes2(:,ii);
                    %             db.muTildes1(:,ii)   = muTildes1(:,ii);
                    %             db.wTildes2(:,ii)   = wTildes2(:,ii);
                    %             db.wTildes1(:,ii)   = wTildes1(:,ii);
                    pFxn(ii) = normcdf(Z);
                    
            end
    end
end

% Convert to 4D mat
temp = reshape(pFxn,[numel(vs1) numel(cs1) numel(vStim2Delta) numel(cs2)]);
ptvs = permute(temp,[2 1 4 3]);

end