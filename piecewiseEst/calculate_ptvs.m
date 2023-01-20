function [ptvs,valMat] = calculate_ptvs(vStim1,cStim1,vStim2Delta,cStim2,avlog,gvlog,hc,method,interpOn)

% Assumptions:
% Exponential prior; Gaussian likelihoods centered on veridical speed
% Velocity encoded in normalized log domain
% Reference contrasts are a subset of test contrasts

% Initialize vars
inds     = [numel(cStim1) numel(vStim1) numel(cStim2) numel(vStim2Delta)];
rf_mu    = nan(inds(1),inds(2));
rf_var   = nan(inds(1),inds(2));
test_mu  = nan(inds(1),inds(2),inds(3),inds(4));
test_var = nan(inds(1),inds(2),inds(3),inds(4));
ptvs     = nan(inds(1),inds(2),inds(3),inds(4));

if interpOn == 2
    % Get indices of nearest neighbor slope estimate for each vStim2
    [avlogInds] = getNNslope(vStim1,vStim2Delta);
end

% for each reference contrast
for rc = 1:inds(1)
    
    rf_cont = cStim1(rc);
    
    % for each reference speed
    for rv = 1:inds(2)
        
        % convert to normalized vlog space
        rf_vlog = log(1 + vStim1(rv)/0.3);
        
        % mean and var of standard posterior in tranformed velocity space
        % hc = all contrasts; assume numel(cStim2) == numel(hc) and 
        % cStim1 is a subset of cStim2
        rf_mu(rc,rv)   = rf_vlog + avlog(rv)*(gvlog(rv)*hc(rf_cont == cStim2))^2;
        rf_var(rc,rv)  = (gvlog(rv)*hc(rf_cont == cStim2))^2;
        
        % define test vels and convert to normalized vlog space
        test_vlog = log(1 + (vStim1(rv)*vStim2Delta)/0.3);
        
        if interpOn == 1
            
            % Interpolate slopes across reference speed for speed range
            % defined by test_vlog
            theseSlopes = interp1(log(1 + vStim1/0.3), avlog, ...
                test_vlog, 'linear', 'extrap');
            
        elseif interpOn == 2
            
            % Use nearest neighbor estimate of slope
            theseSlopes = avlog(avlogInds(rv,:));
            
        end
        
        % for each test contrast
        for tc = 1:inds(3)
        
            % for each test speed
            for ts = 1:inds(4)
                
                if interpOn ~= 0
                    
                    testSlope = theseSlopes(ts);
                    %%%% for debugging
                    tsTest(rc,rv,tc,ts) = testSlope;
                    
                else
                    
                    testSlope = avlog(rv);
                    
                end
               
                test_mu(rc,rv,tc,ts)   = test_vlog(ts) + testSlope*(gvlog(rv)*hc(tc))^2;
                test_var(rc,rv,tc,ts)  = (gvlog(rv)*hc(tc))^2;
                
                % Evaluate p(test > standard) for this test velocity
                
                switch method
                    case 'numeric'
                        
                        v_range = linspace(...
                            min([test_mu(rc,rv,tc,ts) rf_mu(rc,rv)]) - ...
                        3*max(sqrt([test_var(rc,rv,tc,ts) rf_var(rc,rv)])),...
                            max([test_mu(rc,rv,tc,ts) rf_mu(rc,rv)]) + ...
                        3*max(sqrt([test_var(rc,rv,tc,ts) rf_var(rc,rv)])),...
                            5000);
                
                        ptvs(rc,rv,tc,ts) = ...
                            sum( ...
                            (normpdf(v_range,test_mu(rc,rv,tc,ts),sqrt(test_var(rc,rv,tc,ts)))...
                            /sum(normpdf(v_range,test_mu(rc,rv,tc,ts),sqrt(test_var(rc,rv,tc,ts))))).*...
                             cumsum(normpdf(v_range,rf_mu(rc,rv),sqrt(rf_var(rc,rv)))...
                            /sum(normpdf(v_range,rf_mu(rc,rv),sqrt(rf_var(rc,rv))))) );
                       
                    case 'closed'
                        
                        % Get z-score AKA d'/sqrt(2):
                        Z = (test_mu(rc,rv,tc,ts) - rf_mu(rc,rv)) ./ ...
                             sqrt(rf_var(rc,rv) + test_var(rc,rv,tc,ts));
                        
                        % Get p(v2Hat > v1Hat) with phi function:
                        rho = min(max(normcdf(Z), eps), 1 - eps);
                        
                        ptvs(rc,rv,tc,ts) = rho;
                        
                end
                
            end
            
        end
        
    end
    
end

valMat.rf_mu = rf_mu;
valMat.rf_var = rf_var;
valMat.test_mu = test_mu;
valMat.test_var = test_var;
if interpOn
valMat.tsTest = tsTest;
end