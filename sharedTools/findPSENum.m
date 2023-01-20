function [PSE,bias,biasF,pPSE] = findPSENum(vStim1,vStim2Delts,cStim2,avlog,gvlog,hc,hcInds,velInd,plotOn,interpOn)

% Numerically estimate PSE
%%%%%% with closed form solution, can probably just do this analytically. 
%
% Usage: [PSE,bias,biasF] = findPSENum(vStim1,vStim2Delts,cStim2,avlog,gvlog,hc,hcInds,velInd,plotOn)

% Define loss function as absolute difference between 0.5 and ptvs
fixPars.vStim1  = vStim1;
fixPars.avlog   = avlog;
fixPars.gvlog   = gvlog;
fixPars.hc      = hc;
fixPars.hcInds  = hcInds;
fixPars.velInd  = velInd;

lossF = @(vStim2Delta) abs(calcPSE(fixPars,vStim2Delta,interpOn) - 0.5);
opts  = optimset('display','off','tolx',1e-13,'maxfunevals',1e4,'largescale', 'off');
lb = eps;
ub = inf;

% Find PSE (fails sometimes for init = 1?)
try
vDel0 = 1;
vDelHat = fmincon(lossF,vDel0,[],[],[],[],lb,ub,[],opts);
catch
vDel0 = 2;
vDelHat = fmincon(lossF,vDel0,[],[],[],[],lb,ub,[],opts);   
end

PSE = vDelHat*vStim1(fixPars.velInd);
pPSE = calcPSE(fixPars,vDelHat,interpOn);

% Get rest of psychometric FXN
ptvs = calcPSE(fixPars,vStim2Delts,interpOn);
vStim2 = vStim2Delts*vStim1(fixPars.velInd);

% Calculate raw bias and percentage of veridical
bias = PSE-vStim1(fixPars.velInd);
biasF = PSE/vStim1(fixPars.velInd);

% Plot
if plotOn
    
    figure;
    hold on;
    
    plot(vStim2,ptvs);
    scatter(PSE,pPSE,[],[1 0 0],'filled');
    plot([vStim1(fixPars.velInd) vStim1(fixPars.velInd)],[0 1],'--k');
    set(gca,'ylim',[0 1],'xlim',[vStim2(1) vStim2(end)]);
    xlabel('V_{test}');
    ylabel('p(V_{test} > V_{ref})');
    title(['V_{ref} = ',num2str(vStim1(fixPars.velInd)),...
           '; c_{ref} = ',num2str(cStim2(hcInds(1))),...
           '; c_{test} = ',num2str(cStim2(hcInds(2)))]);
end

end

function [pSDT] = calcPSE(fixPars,vStim2Delta,interpOn)
    
    vStim1      = fixPars.vStim1;
    avlog       = fixPars.avlog;
    gvlog       = fixPars.gvlog;
    hc          = fixPars.hc;
    hcInds      = fixPars.hcInds;
    velInds     = fixPars.velInd;
    
    % Define standard pars
    stdV = vStim1(velInds);
    
    st_vlog = log(1 + stdV/0.3);
    st_mu   = st_vlog + avlog(velInds)*(gvlog(velInds)*hc(hcInds(1)))^2;
    st_var  = (gvlog(velInds)*hc(hcInds(1)))^2;
    
    % Define test vels/interpolated slopes if desired
    testV = vStim2Delta*vStim1(velInds);
    test_vlog = log(1 + (testV)/0.3);
    
    if interpOn
        theseSlopes = interp1(log(1 + vStim1/0.3), avlog, ...
                            test_vlog, 'linear', 'extrap');      
    end
    
    for ts = 1:numel(vStim2Delta)
        
        if interpOn
            testSlope = theseSlopes(ts);
        else
            testSlope = avlog(velInds);
        end
        
        % Define test pars
        
        test_mu  = test_vlog(ts) + testSlope*(gvlog(velInds)*hc(hcInds(2)))^2;
        test_var = (gvlog(velInds)*hc(hcInds(2)))^2;

        % Evaluate p(test > standard) for this test velocity
%         v_range = linspace(...
%             min([test_mu st_mu]) - 3*max([test_var st_var]),...
%             max([test_mu st_mu]) + 3*max([test_var st_var]),...
%             100);
%         
%         pSDT(ts) = ...
%             sum( ...
%             (    normpdf(v_range,test_mu,sqrt(test_var))...
%             /sum(normpdf(v_range,test_mu,sqrt(test_var)))).*...
%             cumsum(normpdf(v_range,st_mu,sqrt(st_var))...
%             /sum(       normpdf(v_range,st_mu,sqrt(st_var)))) );
%         
        Z = (test_mu - st_mu) ./ sqrt(st_var + test_var);
                         
        rho = min(max(normcdf(Z), eps), 1 - eps);
                         
        pSDT(ts) = rho;
    end
end