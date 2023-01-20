function [outputMat,expDataMat,vStim2] = reorgSSData(inputMat)

% Rearranges/sorts Stocker & Simoncelli 2006 Data for fitting in our
% version of the analysis pipeline

% Our data arrangement:
% p(VTest>VRef) is stored in 4D matrix with 
% dims (Ref Cont, Ref Vel, Test Cont, Test Vel)

% Their data arrangement:
% 11 x n mat where n is the number of trials
% Row 1: Ref Speed (deg/s)
% Row 2: Ref Cont (%)
% Row 3: Test Speed (deg/s)
% Row 4: Test Cont (%)
% Row 9: Response (Test>Ref? - 0/1)

% Issue arises since they used a Bayesian adaptive staircase to determine
% test vels while we simulated an exp with method of constant stimuli

%% Sorting
% Get unique vals for each stimulus parameter
uniqRefVels     = unique(inputMat(1,:));
uniqRefConts    = unique(inputMat(2,:));
uniqTestConts   = unique(inputMat(4,:));

numRefVels      = numel(uniqRefVels);
numRefConts  	= numel(uniqRefConts);
numTestConts    = numel(uniqTestConts);

% Get vectors of logical inds
for rv = 1:numRefVels
    refVelInds(:,rv)   = [inputMat(1,:) == uniqRefVels(rv)]';
end

for rc = 1:numRefConts
    refContInds(:,rc)  = [inputMat(2,:) == uniqRefConts(rc)]';
end

for tc = 1:numTestConts
    testContInds(:,tc) = [inputMat(4,:) == uniqTestConts(tc)]';
end

% Organize into 3D structure for now with test vels and responses as
% subfields

data   = struct;
vStim2 = cell(numRefConts,numRefVels,numTestConts);

for rc = 1:numRefConts
    for rv = 1:numRefVels
        for tc = 1:numTestConts
            
            % Logical index for unique triplets of (rc,rv,tc)
            theseInds = refVelInds(:,rv) & ...
                        refContInds(:,rc) & ...
                        testContInds(:,tc);
            
            % All 80 trials for this triplet
            data(rc,rv,tc).testVels      = inputMat(3,theseInds);
            data(rc,rv,tc).responses     = inputMat(9,theseInds);
            
            % Find unique test velocities
            data(rc,rv,tc).uniqTestVels  = unique(inputMat(3,theseInds));
            vStim2{rc,rv,tc} = data(rc,rv,tc).uniqTestVels;
            
            % Find proportion of test>ref for each unique vel
            for tv = 1:numel(data(rc,rv,tc).uniqTestVels)
                
                vTemp     = data(rc,rv,tc).uniqTestVels(tv);
                theseVels = data(rc,rv,tc).testVels == vTemp;
                
                rTemp = data(rc,rv,tc).responses;
                
                data(rc,rv,tc).numCh(tv)   = sum(rTemp(theseVels));
                data(rc,rv,tc).numTr(tv)  = numel(rTemp(theseVels));
                 
            end
            
        end
    end
end
            
outputMat = data;            

expDataMat.uniqRefVels      = uniqRefVels;
expDataMat.uniqRefConts     = uniqRefConts;
expDataMat.uniqTestConts    = uniqTestConts;
expDataMat.numRefVels       = numRefVels;
expDataMat.numRefConts      = numRefConts;
expDataMat.numTestConts     = numTestConts;


end