function [fitStr] = fitSSDatags(subjData,outDir,outFName,solver)

% Takes original datasets from SS2006 and fits them with our reconstruction
% of their model
%
% Usage: [fitStr] = fitSSDatags(subjData)

%% Rearrange data raw data mat so we can work with it

[datStr,expData,vStim2] = reorgSSData(subjData);

vStim1          = expData.uniqRefVels;
cStim1          = expData.uniqRefConts;
cStim2          = expData.uniqTestConts;

numRefVels      = expData.numRefVels;
numRefConts     = expData.numRefConts;
numTestConts    = expData.numTestConts;
numPars         = numRefVels + numTestConts + 1;


%% Setup anonymous loss function & fitting options

lossFxn = @(parVec) calculateModelNllSSDatags(parVec,datStr,vStim1,cStim1,cStim2,vStim2);

opts    = optimset('useparallel',true,'display','off','tolx',1e-13,...
                   'maxfunevals',1e4,'largescale', 'off');
lb      = [zeros(1,numRefVels) zeros(1,numTestConts) 0];
ub      = inf*[ones(1,numRefVels) ones(1,numTestConts) 1];


%% Fit data

numReps = 25;
indCnt  = 1;

sigPF       = nan(numReps,1);
gvlogF      = nan(numReps,numRefVels);
hcF         = nan(numReps,numTestConts);
ptvsF       = cell(numReps,1);
nllF        = nan(numReps,1);
initConds   = nan(numReps,numPars);
os          = cell(numReps,1);

while indCnt <= numReps
    
    disp(['Run ',num2str(indCnt),'/',num2str(numReps)]);
    
    tsl = tic;
    
    % Randomize initial parameters within reasonable window
    gvloginit       = rand(1,numRefVels)*(2-0.1)  + 0.1;
    hcinit          = rand(1,numTestConts)*(2-0.1) + 0.1;
    sigPinit        = rand(1,1)*(1-0.1) + 0.1;
    
    parVec0         = [gvloginit hcinit sigPinit];
    
    % Run fitting using selected solver
    switch solver
        case 'fminunc'
            [prshat,~,exFlag,outStruc]  = fminunc(lossFxn,parVec0,opts);
        case 'fminsearch'
            [prshat,~,exFlag,outStruc]  = fminsearch(lossFxn,parVec0,opts);
        case 'fmincon'
            [prshat,~,exFlag,outStruc]  = fmincon(lossFxn,parVec0,[],[],[],[],lb,ub,[],opts);
    end
    
    gvlogF(indCnt,:) = prshat(1:numRefVels);
    hcF(indCnt,:)    = prshat(numRefVels + 1:numRefVels + numTestConts);
    sigPF(indCnt,:)   = prshat(end);

    initConds(indCnt,:) = parVec0;
    os{indCnt}          = outStruc;
    
    % Output best fitting psychometric function and nLL
    ptvsF{indCnt}     = calcPtvsSSDatags(cStim1,vStim1,cStim2,vStim2,...
                                         sigPF(indCnt,:),gvlogF(indCnt,:),hcF(indCnt,:));
    nllF(indCnt)      = lossFxn(prshat);
    
    % Only increment if optimization didn't get stuck in a black hole
    % Fminsearch almost always hits ceiling so don't even bother.
    if (exFlag > 0) || ~(strcmp(os(indCnt),'quasi-newton'))
        indCnt = indCnt + 1;
    elseif exFlag == 0
        disp('Max function evals hit, repeat loop.');
    end
    
    t = toc(tsl);
    disp(['Time Elapsed: ',num2str(t)]);
    
end


%% Save fits
saveFPath = [outDir,filesep,outFName,'_',solver,'_gs'];

expStr.vStim2 = vStim2;
save(saveFPath,'initConds','gvlogF','hcF','sigPF','ptvsF','nllF',...
               'datStr','expStr','os','vStim1','cStim1','cStim2');

fitStr.initConds = initConds;
fitStr.gvlogF    = gvlogF;
fitStr.hcF       = hcF;
fitStr.sigPF    = sigPF;
fitStr.ptvsF     = ptvsF;
fitStr.nllF      = nllF;
fitStr.datStr    = datStr;
fitStr.expStr    = expStr;
fitStr.os        = os;


end