function [simData,s1Conts,s2Conts,s1Vels,s2Vels_relative,ptvs_CF] = autoGenGroundTruth(sDir,fName,parStr)
% Generate multiple ground truth datasets without UI
%
% Usage: [simData,s1Conts,s2Conts,s1Vels,s2Vels_relative,ptvs_CF] = autoGenGroundTruth(sDir,fName,parStr)


%% Define/create output folder
saveDir = [sDir,filesep,fName];

if exist(saveDir,'dir') == 0
    mkdir(saveDir);
end

%% Define stimulus parameters

s1Vels          = parStr.s1Vels;
s1Conts         = parStr.s1Conts;
s2Vels_relative = parStr.s2VelsDelta;
s2Conts         = parStr.s2Conts;
trPoolSize      = parStr.trPoolSize;

s2Vels = s1Vels'*s2Vels_relative;

% Get unique contrast values (differs fom Tyler's indexing) - this is
% considered to be == s2Conts later in the code & that is true as currently
% written
[allConts,~,~] = unique([s1Conts s2Conts]);

%% Define Observer parameters

priorF     = parStr.priorF;
priorStruc = parStr.priorStruc;
contScF    = parStr.contScF;
domain     = parStr.domain;

prInd      = parStr.priorInd;
LRind      = parStr.LRind;

lapseRate  = parStr.lapseRate(LRind);

switch domain
    case 'linear'
        stimDom = 'lin';
    case 'log'
        stimDom = 'log';
end

% Save ground truth for later & determine likelihood sigmas
switch priorF
    
    case 'gauss'
        
        % Gaussian width
        sigP = priorStruc.sigP(prInd);
        
        priorPar(1) = priorStruc.muP(prInd);
        priorPar(2) = sigP;
        
    case 'expon'  
        
        % Store slopes for fitting
        slope = priorStruc.avlog(prInd);
        avlog = repmat(slope,1,numel(s1Vels));
        
        priorPar(1) = slope;
        priorPar(2) = priorStruc.bP(prInd);
        
    case 'mixgauss'
                
        % Gaussian widths
        gammaP = priorStruc.gammaP{prInd};
        
        % Gaussian means
        muP    = priorStruc.muP{prInd};
        
        % Component weights
        wP     = priorStruc.wP{prInd};
        
        priorPar{1} = [gammaP;muP;wP];
        
    case 'splitExp'
        
        % Store slopes for fitting
        aP    = priorStruc.aP{prInd};
        avlog = repelem(aP,1,numel(s1Vels)/numel(aP));
        
        priorPar{1} = aP;
        priorPar{2} = priorStruc.bP{prInd};
        
    case 'genGauss'
        
        % Store generalized Gaussian pars
        uP = priorStruc.uP;
        aP = priorStruc.aP;
        bP = priorStruc.bP;
        
        priorPar{1} = [uP;aP;bP];
        
end

% Speed-dependent scale factor for likelihood width (flat for now)
gvlog = repmat(0.2,1,numel(s1Vels));

% Contrast-dependent scale factor for likelihood width (2 for SS)
hc = -allConts + contScF;

% Get likelihood SDs (reference matrix for s1 and s2)
for v = 1:numel(s1Vels)
    for c = 1:numel(allConts)
        sSD(v,c) = gvlog(v)*hc(c);
    end
end


%% Generate dataset

numS1Vels = numel(s1Vels);
numS1Conts = numel(s1Conts);
numS2Conts = numel(s2Conts);

indCnt = 1;

% Loop over reference stimulus contrasts
for sc = 1:numS1Conts

    % Get index into allConts that corresponds to this contrast
    % level of the standard
    thisContrastS1 = (allConts == s1Conts(sc));
    
    % Loop over reference stimulus velocities
    for sv = 1:numS1Vels
        
        % Loop over test stimulus contrasts
        for tc = 1:numS2Conts
            
            % Get index into allConts that corresponds to this contrast
            % level of the test
            thisContrastS2 = (allConts == s2Conts(tc));
            
            % Grab parameters for each likelihood
            like1Par = [s1Vels(sv) sSD(sv,thisContrastS1)];
            like2Par = {s2Vels(sv,:),sSD(sv,thisContrastS2)};
            
            % Generate ground truth data for these parameters
            [tempStruct] = ...
                genGroundTruthData(priorF,priorPar,like1Par,like2Par,...
                                   stimDom,trPoolSize,lapseRate);
            
            simData(indCnt) = tempStruct;
            
            % With lapse rate
            ptvs_GT(sc,sv,tc,:) = simData(indCnt).simProps;
            
            % Calculated with ss/gs model
            ptvs_CF(sc,sv,tc,:) = simData(indCnt).cf_pSDT;
            
            indCnt = indCnt + 1;
        end
    end
end

% Package up GT pars into vector
switch priorF
    
    case 'expon'
        
        parVec = [gvlog hc avlog];
        
    case 'gauss'
        
        parVec = [gvlog hc sigP];
        
    case 'mixgauss'
        
        parVec = [gvlog hc gammaP muP wP];
        
    case 'splitExp'
        
        parVec = [gvlog hc avlog];
        
    case 'genGauss'
        
        parVec = [gvlog hc uP aP bP];
        
end

% Do some relabeling
vStim1      = s1Vels;
cStim1      = s1Conts;
vStim2Delta = s2Vels_relative;
cStim2      = s2Conts;

%% Save dataset
simDataPath = [saveDir,filesep,'simData'];
save(simDataPath,'simData','parVec','vStim1','cStim1','vStim2Delta',...
    'cStim2','trPoolSize','ptvs_GT','ptvs_CF');


end
