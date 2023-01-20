function [ptvs_data,numSamps,ptvs_trData] = choiceResamplerGauss(cStim1,cStim2,vStim1,vStim2Delta,numTrialsFull,simData,numTrials)

% Takes a random subsample of n choices from full dataset (Gaussian sampling
% from set of test speeds)

%% Set up sampling 
gauF = @(x,mu,sig) (1/(sig*sqrt(2*pi)))*exp(-0.5*((x-mu)/sig).^2);
inds2sig = linspace(-2.5,2.5,10);
dist = gauF(inds2sig,0,1);
numSamps = round(numTrials * (dist/sum(dist)));

% Make sure each test index has at least one sample (be careful this
% doesn't cause total numTrials to deviate - maybe decrease sigma)
numSamps(numSamps == 0) = 1;

%% Set up vars
numS1Conts = numel(cStim1);
numS2Conts = numel(cStim2);
numS1Vels = numel(vStim1);
numS2Vels = numel(vStim2Delta);

ptvs_data = nan(numS1Conts,numS1Vels,numS2Conts,numS2Vels);

indCnt = 1;

for sc = 1:numS1Conts
    
    for sv = 1:numS1Vels
        
        for tc = 1:numS2Conts
            
            randTrInds = datasample(1:numTrialsFull,max(numSamps),'Replace',false); 
            
            % Grab a block of trials
            ptvsTemp = simData(indCnt).simChoices(randTrInds,:);
            
            % Loop over number of test vels
            for tv = 1:numS2Vels
                
                % Grab number of trials according to gaussian profile
                testTempInds = datasample(1:max(numSamps),numSamps(tv),'Replace',false);
                
                ptvs_trData{sc,sv,tc,tv} = ptvsTemp(testTempInds,tv);
                ptvs_data(sc,sv,tc,tv) = mean(ptvsTemp(testTempInds,tv));
                
            end
            
            indCnt = indCnt + 1;
            
        end
    end
end

end