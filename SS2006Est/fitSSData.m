function [fitStr] = fitSSData(subjData,outDir,outFName,solver)

% Takes original datasets from SS2006 and fits them with our reconstruction
% of their model
%
% Usage: [fitStr] = fitSSData(subjData)

%% Rearrange data raw data mat so we can work with it

[datStr,expData,vStim2] = reorgSSData(subjData);

vStim1          = expData.uniqRefVels;
cStim1          = expData.uniqRefConts;
cStim2          = expData.uniqTestConts;

numRefVels      = expData.numRefVels;
numRefConts     = expData.numRefConts;
numTestConts    = expData.numTestConts;
numPars         = numRefVels + numTestConts;


%% Setup anonymous loss function & fitting options

lossFxn = @(parVec) calculateModelNllSSData(parVec,datStr,vStim1,cStim1,cStim2,vStim2);
% Try using Naecker's method
% lossFxn = @(parVec) calculateModelNllSSData_altCalc(parVec,subjData,vStim1);

opts    = optimset('useparallel',true,'display','off','tolx',1e-13,...
                   'maxfunevals',1e4,'largescale', 'off');
% lb      = [-inf*ones(1,numRefVels) zeros(1,numRefVels) zeros(1,numTestConts)];
% ub      = inf*[ones(1,numRefVels) ones(1,numRefVels) ones(1,numTestConts)];
lb = [-100*ones(1,numRefVels) 0.001*ones(1,numTestConts)];
ub = [100*ones(1,numRefVels) 100*ones(1,numTestConts)];


%% Fit data

% numReps = 25;
numReps = 1;
indCnt  = 1;

avlogF      = nan(numReps,numRefVels);
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
%     avloginit       = rand(1,numRefVels)*(5-(-5)) +  -5;
%     gvloginit       = rand(1,numRefVels)*(2-0.1)  + 0.1;
%     hcinit          = rand(1,numTestConts)*(2-0.1) + 0.1;
    avloginit       = 1*ones(1,numRefVels);
    hcinit          = 2*ones(1,numTestConts);
    
%     parVec0         = [avloginit gvloginit hcinit];
    parVec0         = [avloginit hcinit];
    
    % Run fitting using selected solver
    switch solver
        case 'fminunc'
            [prshat,~,exFlag,outStruc]  = fminunc(lossFxn,parVec0,opts);
        case 'fminsearch'
            [prshat,~,exFlag,outStruc]  = fminsearch(lossFxn,parVec0,opts);
        case 'fmincon'
            [prshat,~,exFlag,outStruc]  = fmincon(lossFxn,parVec0,[],[],[],[],lb,ub,[],opts);
    end
    
%     avlogF(indCnt,:) = prshat(1:numRefVels);
%     gvlogF(indCnt,:) = prshat(numRefVels + 1:2*numRefVels);
%     hcF(indCnt,:)    = prshat(2*numRefVels + 1:numPars);
    avlogF(indCnt,:) = prshat(1:numRefVels);
    gvlogF(indCnt,:) = 0.2*ones(1,numRefVels);
    hcF(indCnt,:)    = prshat(numRefVels + 1:numPars);
    
    initConds(indCnt,:) = parVec0;
    os{indCnt}          = outStruc;
    
    % Output best fitting psychometric function and nLL
    ptvsF{indCnt}     = calcPtvsSSData(cStim1,vStim1,cStim2,vStim2,...
                                       avlogF(indCnt,:),gvlogF(indCnt,:),hcF(indCnt,:));
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
saveFPath = [outDir,filesep,outFName,'_',solver];

expStr.vStim2 = vStim2;
save(saveFPath,'initConds','gvlogF','hcF','avlogF','ptvsF','nllF',...
               'datStr','expStr','os','vStim1','cStim1','cStim2');

fitStr.initConds = initConds;
fitStr.gvlogF    = gvlogF;
fitStr.hcF       = hcF;
fitStr.avlogF    = avlogF;
fitStr.ptvsF     = ptvsF;
fitStr.nllF      = nllF;
fitStr.datStr    = datStr;
fitStr.expStr    = expStr;
fitStr.os        = os;


end