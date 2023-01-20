function [ptvs_data,ptvs_rawdata,likeMeans1,likeMeans2] = choiceResampler(cStim1,cStim2,vStim1,vStim2Delta,numTrialsFull,simData,numTrials)

% Takes a random subsample of n choices from full dataset (evenly sampling
% from each test speed)

numS1Conts = numel(cStim1);
numS2Conts = numel(cStim2);
numS1Vels = numel(vStim1);
numS2Vels = numel(vStim2Delta);

ptvs_data = nan(numS1Conts,numS1Vels,numS2Conts,numS2Vels);

nftinit = numTrialsFull;

indCnt = 1;

for rc = 1:numS1Conts
    
    for rv = 1:numS1Vels
        
        for tc = 1:numS2Conts
            
            % hacky bit for debugging piecewise
            if isfield(simData(indCnt),'incInds1')
                
                % Sample from included inds
                incChoices    = simData(indCnt).incInds1 & simData(indCnt).incInds2;
                numTrialsFull = sum(incChoices);
                
                if ~any(numTrialsFull)
                    error('No trials fit inclusion criteria.');
                end
                
                for ii = 1:numS2Vels
                    
                    if numTrials > numTrialsFull(ii)
                        error('Not enough trials fit inclusion criteria.');
                    end
                    
                    indsoinds = datasample(1:numTrialsFull(ii),numTrials,'Replace',false); 
                    incInds   = find(incChoices(:,ii));
                    
                    randTrInds(:,ii) = sub2ind([nftinit,numS2Vels],incInds(indsoinds),ii*ones(numTrials,1));
                    
                end
                
            else
                
                randTrInds    = datasample(1:numTrialsFull,numTrials,'Replace',false);
                
            end
            
            % Store noisy likelihood means for each trial for debugging
%             likeMeans1{rc,rv} = simData(indCnt).likeMeans1(randTrInds,1);
%             lm2Mat            = simData(indCnt).likeMeans2(randTrInds,:);
            likeMeans1{rc,rv} = [];            
            
            % Get matrix of choices for each test velocity for this triplet
            %%%% get rid of test below after debugging
            if ~isrow(randTrInds)
                
                tvMat         = simData(indCnt).simChoices(randTrInds);
                
            else
                tvMat         = simData(indCnt).simChoices(randTrInds,:);
            end
           
            % Subsample and mean to convert to proportion of trials s2 > s1
            ptvs_data(rc,rv,tc,:) = mean(tvMat);
            
            % Get all repeats for this quadruplet
            for tv = 1:numS2Vels
                ptvs_rawdata{rc,rv,tc,tv}   = tvMat(:,tv);
                
%                 likeMeans2{rc,rv,tc,tv}     = lm2Mat(:,tv);
                likeMeans2{rc,rv,tc,tv} = [];
            end
            
            indCnt = indCnt + 1;
            
        end
    end
end

end