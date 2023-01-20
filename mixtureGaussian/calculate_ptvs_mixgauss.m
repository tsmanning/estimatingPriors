function [ptvs,db] = calculate_ptvs_mixgauss(vStim1,cStim1,vStim2Delta,cStim2,gvlog,hc,sigP,muP,wP,method,estDistEval)

% Assumptions:
% Prior is a mixture of n Gaussian components with indivdual means and
% variances; velocity encoded in normalized log domain
%
% Usage: [ptvs] = calculate_ptvs_mixgauss(vStim1,cStim1,vStim2Delta,cStim2,gvlog,hc,sigP,muP,wP)

% Initialize vars
inds     = [numel(cStim1) numel(vStim1) numel(cStim2) numel(vStim2Delta)];
ptvs     = nan(inds(1),inds(2),inds(3),inds(4));

% Loop over reference contrasts
for rc = 1:inds(1)
    
    rf_cont = cStim1(rc);
    
    % Loop over reference speeds
    for rv = 1:inds(2)
        
        % convert to normalized vlog space
        rf_vlog   = log(1 + vStim1(rv)/0.3);
        
        % MAP sampling distribution is product of likelihood EV and prior
        %%% index into sigL according to contrasts
        thisRSigL = gvlog(rv)*hc(rf_cont == cStim2);
        
        [alphas1,muTildes1,wTildes1] = getMogPostPars(wP,sigP,muP,thisRSigL,rf_vlog);
        db.alphas1{rc,rv} = alphas1;
        db.muTildes1{rc,rv} = muTildes1;
        db.wTildes1{rc,rv} = wTildes1;
        
        % Loop over test contrasts
        for tc = 1:inds(3)
            
            % Loop over test speeds
            for ts = 1:inds(4)
                
                test_velDelta = vStim2Delta(ts);
                
                % define test vels and convert to normalized vlog space
                test_vlog     = log(1 + (vStim1(rv)*test_velDelta)/0.3);
                
                thisTSigL     = gvlog(rv)*hc(tc);
                
                [alphas2,muTildes2,wTildes2] = getMogPostPars(wP,sigP,muP,thisTSigL,test_vlog);
                db.alphas2{rc,rv,tc,ts} = alphas2;
                db.muTildes2{rc,rv,tc,ts} = muTildes2;
                db.wTildes2{rc,rv,tc,ts} = wTildes2;
                
                % Evaluate p(test > standard) for this test velocity
                
%                 estDistEval = 'original';
%                 estDistEval = 'rederived';
                
                switch method
                    
                    case 'numeric'
                        
                        mu1  = [];
                        mu2  = [];
                        var1 = [];
                        var2 = [];
                        
                        % Define support
                         vRange = linspace(...
                            min([mu1 mu2]) - ...
                        3*max(sqrt([test_var(rc,rv,tc,ts) rf_var(rc,rv)])),...
                            max([test_mu(rc,rv,tc,ts) rf_mu(rc,rv)]) + ...
                        3*max(sqrt([test_var(rc,rv,tc,ts) rf_var(rc,rv)])),...
                            5000);
                        
                        % Define estimate distribution stim 1
                        
                        sPost1 = zeros(1,numel(vRange));
                        
                        for ii = 1:numel(alphas1)
                            
                            sPost1 = sPost1 + ...
                                     wTildes2(ii)*normpdf(vRange,...
                                                          muTildes1(ii) + alphas1(ii)*ref_vlog,...
                                                          alphas1(ii)*thisRSigL);
                            
                        end
                        
                        sEst1 = [];
                            
                        % Define estimate distribution stim 2
                        
                        sPost2 = zeros(1,numel(vRange));
                        
                        for ii = 1:numel(alphas2)
                        
                            sPost2 = sPost2 + [];
                            
                        end
                        
                        sEst2 = [];
                        
                        % Sum over distributions
                        ptvs(rc,rv,tc,ts) = sum(sEst1.*cumsum(sEst2));
                
                    case 'closed'
                        
                        % This would all be better handled outside of a
                        % loop...
                        
                        switch estDistEval
                            case 'original'
                                denom   = sqrt((thisTSigL^2)*(alphas2'.^2) + (thisRSigL^2)*(alphas1'.^2));
                                
                                db.denom{rc,rv,tc,ts} = denom;
                                
                                % Careful, this is a matrix: Test dim 1, ref dim 2
                                Z       = (muTildes2' + alphas2'*test_vlog - ...
                                           muTildes1 - alphas1*rf_vlog)./denom;
                                
                                db.Z{rc,rv,tc,ts} = Z;
                                
                                % Again careful, order matters
                                ptvs(rc,rv,tc,ts)    = wTildes2'*normcdf(Z)*wTildes1;
                                
                            case 'rederived'
                                denom   = sqrt((thisTSigL^2)*sum(wTildes2.*alphas2)^2 + ...
                                               (thisRSigL^2)*sum(wTildes1.*alphas1)^2);
                                
                                Z       = ( sum(wTildes2.*(muTildes2 + alphas2*test_vlog)) - ...
                                            sum(wTildes1.*(muTildes1 + alphas1*rf_vlog  )) )/denom;
                                
                                db.denom{rc,rv,tc,ts} = denom;
                                db.Z{rc,rv,tc,ts}    = Z;
                                ptvs(rc,rv,tc,ts)    = normcdf(Z);
                                
                        end
                        
                end
                
            end
            
        end
        
    end
    
end

end
