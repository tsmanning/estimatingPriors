function instrStr = instr_mgExpon

% Set instructions for:

% Mixture of Gaussians fit
% Mixture of Gaussians Prior
% Bernoulli estimation

%% Ground Truth Observer Options

% Set domain in which velocity is encoded ('log'/'linear')
domain      = 'log';

% Set large number of trials that approaches "ground truth" precision for
% numerical calculations
trPoolSize  = 5000;

% Set widths of observer prior (either slope of exponential or sigma of
% gaussian)
% priorF      = 'gauss';
priorF      = 'expon';
% priorF      = 'splitExp';
% priorF      = 'mixgauss';

switch priorF
    case 'gauss'
        % Gaussian
        priorStruc.sigP      = linspace(0.5,1.5,5);
        priorStruc.numPriors = numel(priorStruc.sigP);
        
    case 'expon'
        % Global exponential
        priorStruc.avlog     = linspace(-8,-1.5,5);
%         priorStruc.avlog     = linspace(-8,-1.5,1);
        priorStruc.bP        = zeros(numel(priorStruc.avlog),1);
        priorStruc.numPriors = numel(priorStruc.avlog);
        
    case 'splitExp'
        % Split Exp
        priorStruc.aP        = {[-3 -1.5],[-1.5 -3]};
        priorStruc.bP        = {0,0};
        priorStruc.numPriors = numel(priorStruc.aP);
        priorStruc.splitLoc  = 3;
        
    case 'mixgauss'
        % Mixture Gaussian [sigmas; mus; weights] - otta find a way to package this
        % up to allow multiple different priors and still work with other types of
        % priors - each prior is one 3xN matrix
        priorStruc.gammaP    = {[0.5 0.5 0.5],[1 1 1],[1.5 1.5 1.5],...
                                [0.5 1 1.5],...
                                [1 1 1]};
        priorStruc.muP       = {[0 0 0],[0 0 0],[0 0 0],...
                                [0 0 0],...
                                [0 0.25 0.5]};
        priorStruc.wP        = {[1 1 1]/3,[1 1 1]/3,[1 1 1]/3,...
                                [1 1 1]/3,...
                                [1 1 1]/3};
        priorStruc.numPriors = numel(priorStruc.gammaP);
end

% Set observer lapse rates to simulate
% lapseRate = [0 0.05];
lapseRate   = [0];

%% Experiment Options

% Define Reference Stimuli
s1Vels      = [0.5 1 2 4 8 12]; % velocity (deg/sec)
% s1Conts     = [0.075 0.5];      % contrast (% max)
s1Conts     = 0.5;

% Define Test Stimuli
s2VelsDelta = linspace(0.5,2,10);
s2Conts     = [0.05 0.075 0.1 0.2 0.4 0.5 0.8];
% s2Conts     = [0.05 0.1 0.2 0.5 0.8];


%% Estimation Options

% Set encoding/decoding model
% model     = 'gs';
% model     = 'ss';
model     = 'mg';

% Set number of trials in sample
% numTrials   = [5 10 20];
numTrials   = [1];

% Resample how many times?
numRuns     = 25;
% numRuns     = 1;

% How many randomized parameter vectors to initialize fmincon with?
numReps     = 15;
% numReps     = 2;

% Use trial-by-trial Bernoulli calculation for NLL instead of Binomial?
tbtOn       = 1;

% Choose fit method: numerical or closed form
% fitmethod = 'numeric';
fitmethod   = 'closed';

% Interpolate slopes (0: off, 1:interpolation, 2:nearest neightbor)
interpOn    = 2;

% Scale binomial NLL calc to match bernoulli?
binScale    = 0;

% Set contrast-dependent scale factor (keep)
contScF     = 2;
% contScF     = 1.2;

% Set number of Gaussians you want to estimate prior with
numComps    = 3;

%% Package 'em up
instrStr.domain         = domain;
instrStr.trPoolSize     = trPoolSize;
instrStr.priorF         = priorF;
instrStr.priorStruc     = priorStruc;
instrStr.lapseRate      = lapseRate;
instrStr.s1Vels         = s1Vels;
instrStr.s1Conts        = s1Conts;
instrStr.s2VelsDelta    = s2VelsDelta;
instrStr.s2Conts        = s2Conts;
instrStr.numTrials      = numTrials;

% Fitting model
instrStr.model          = model;
instrStr.numRuns        = numRuns;
instrStr.numReps        = numReps;
instrStr.tbtOn          = tbtOn;
instrStr.fitmethod      = fitmethod;
instrStr.interpOn       = interpOn;
instrStr.binScale       = binScale;
instrStr.contScF        = contScF;
instrStr.numComps       = numComps;

end