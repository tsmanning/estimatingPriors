function [ptvs] = calcPtvsSSDatags(cStim1,vStim1,cStim2,vStim2,sigP,gvlog,hc)

% Calculate p(s2>s1) for the S&S 2006 dataset

% cStim1: n x 1 double
% vStim1: m x 1 double
% cStim2: p x 1 double
% vStim2: n x m x p cell


%% Initialize

expCnts = [numel(cStim1) numel(vStim1) numel(cStim2)];
ptvs    = cell(expCnts(1),expCnts(2),expCnts(3));

%% Calculate psychometric functions for input parameters

% For each reference contrast
for rc = 1:expCnts(1)
    
    rfCont = cStim1(rc);
    
    % For each reference speed
    for rv = 1:expCnts(2)
        
        % Convert to normalized vlog space
        rfVel  = log(1 + vStim1(rv)/0.3);        
        
        % mean and var of standard posterior in tranformed velocity space
        % hc = all contrasts; assume numel(cStim2) == numel(hc) and 
        % cStim1 is a subset of cStim2
        thisRSigL     = gvlog(rv)*hc(rfCont == cStim2);
        [refMu,rfStd] = findGaussProd(rfVel,thisRSigL,0,sigP);
        refVar        = rfStd^2;
        
        % For each test contrast
        for tc = 1:expCnts(3)
            
            % Define test velocities
            testVels    = vStim2{rc,rv,tc};
            
            % Convert to normalized log space
            testVelLog  = log(1 + testVels/0.3);
            
            thisTSigL        = gvlog(rv)*hc(tc);
            [testMu,testStd] = findGaussProd(testVelLog,thisTSigL,0,sigP);
            testVar          = testStd^2;
            
            ptvsTemp = nan(numel(testMu),1);
            
            % For each test velocity
            for ts = 1:numel(testMu)
            
                % evaluate p(test > standard) for this test velocity
                v_range = linspace(...
                    min([testMu(ts) refMu]) - 3*max([testVar refVar]),...
                    max([testMu(ts) refMu]) + 3*max([testVar refVar]),...
                    100);
                
                ptvsTemp(ts) = ...
                    sum( ...
                        (normpdf(v_range,testMu(ts),sqrt(testVar))...
                    /sum(normpdf(v_range,testMu(ts),sqrt(testVar)))).*...
                         cumsum(normpdf(v_range,refMu,sqrt(refVar))...
                           /sum(normpdf(v_range,refMu,sqrt(refVar)))) );
            
            end
            
            ptvs{rc,rv,tc} = ptvsTemp;
            
            clear testVelLog testMu testVar
                  
        end
        
    end
    
end


end