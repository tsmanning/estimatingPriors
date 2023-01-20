function [ptvs] = calculate_ptvs_gauss(vStim1,cStim1,vStim2Delta,cStim2,gvlog,hc,sigP,method)

% Assumptions:
% Gaussian prior and likelihoods centered on veridical speed
% Velocity encoded in normalized log domain
% Reference contrasts are a subset of test contrasts

%%%% fix here, need to precompute mu/var for each, THEN calculate 

% Initialize vars
inds     = [numel(cStim1) numel(vStim1) numel(cStim2) numel(vStim2Delta)];
rf_mu    = nan(inds(1),inds(2));
rf_var   = nan(inds(1),inds(2));
test_mu  = nan(inds(1),inds(2),inds(3),inds(4));
test_var = nan(inds(1),inds(2),inds(3),inds(4));
ptvs     = nan(inds(1),inds(2),inds(3),inds(4));

% Loop over reference contrasts
for rc = 1:inds(1)
    
    rf_cont = cStim1(rc);
    
    % Loop over reference speeds
    for rv = 1:inds(2)
        
        % convert to normalized vlog space
        rf_vlog = log(1 + vStim1(rv)/0.3);
        
        % MAP sampling distribution is product of likelihood EV and prior
        %%% index into sigL according to contrasts
        thisRSigL             = gvlog(rv)*hc(rf_cont == cStim2);
%         [rf_mu(rc,rv),~]      = findGaussProd(rf_vlog,thisRSigL,0,sigP);
        scaleFac_ref          = (sigP^2)/(sigP^2 + thisRSigL^2);
        rf_var(rc,rv)         = (scaleFac_ref^2)*(thisRSigL^2);
        rf_mu(rc,rv)          = scaleFac_ref*rf_vlog;
                
        % Loop over test contrasts
        for tc = 1:inds(3)
                    
            % Loop over test speeds
            for ts = 1:inds(4)
                
                test_velDelta = vStim2Delta(ts);
                
                % define test vels and convert to normalized vlog space
                test_vlog     = log(1 + (vStim1(rv)*test_velDelta)/0.3);
                
                thisTSigL                       = gvlog(rv)*hc(tc);
%                 [test_mu(rc,rv,tc,ts),~]        = findGaussProd(test_vlog,thisTSigL,0,sigP);
                scaleFac_test                   = (sigP^2)/(sigP^2 + thisTSigL^2);
                test_var(rc,rv,tc,ts)           = (scaleFac_test^2)*(thisTSigL^2);
                test_mu(rc,rv,tc,ts)            = scaleFac_test*test_vlog;
                
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
                        
                        Z = (test_mu(rc,rv,tc,ts) - rf_mu(rc,rv)) ./ ...
                             sqrt(rf_var(rc,rv) + test_var(rc,rv,tc,ts));
                         
                        rho = min(max(normcdf(Z), eps), 1 - eps);
                         
                        ptvs(rc,rv,tc,ts) = rho;
                        
                end
                
            end
            
        end
        
    end
    
end