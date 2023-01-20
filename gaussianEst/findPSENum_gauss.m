function [PSE,bias,biasF,pPSE] = findPSENum_gauss(vStim1,vStim2Delts,cStim2,sigP,gvlog,hc,hcInds,velInd,plotOn)

% Numerically estimate PSE
%
% Usage: [PSE,bias,biasF] = findPSENum_gauss(vStim1,vStim2Delts,cStim2,sigP,gvlog,hc,hcInds,velInd,plotOn)
% - hcInds: 1x2 vector - [RefContrastInd TestContrastInd]
%           INDS OF UNIQUE CONTRAST VECTOR
% - velInd: scalar - RefVelInd

% Define loss function as absolute difference between 0.5 and ptvs
fixPars.vStim1  = vStim1;
fixPars.sigP    = sigP;
fixPars.gvlog   = gvlog;
fixPars.hc      = hc;
fixPars.hcInds  = hcInds;
fixPars.velInd  = velInd;

lossF = @(vStim2Delta) abs(calcPSE(fixPars,vStim2Delta) - 0.5);
opts  = optimset('display','off','tolx',1e-13,'maxfunevals',1e4,'largescale', 'off');

% Find PSE (fails sometimes for init = 1?)
try
vDel0 = 1;
vDelHat = fminunc(lossF,vDel0,opts);
catch
vDel0 = 2;
vDelHat = fminunc(lossF,vDel0,opts);   
end

PSE = vDelHat*vStim1(fixPars.velInd);
pPSE = calcPSE(fixPars,vDelHat);

% Get rest of psychometric FXN
ptvs = calcPSE(fixPars,vStim2Delts);
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

function [pSDT] = calcPSE(fixPars,vStim2Delta)
    
    vStim1      = fixPars.vStim1;
    sigP        = fixPars.sigP;
    gvlog       = fixPars.gvlog;
    hc          = fixPars.hc;
    hcInds      = fixPars.hcInds;
    velInds     = fixPars.velInd;
    
    % Define standard pars
    stdV = vStim1(velInds);
    
    st_vlog        = log(1 + stdV/0.3);
    [st_mu,st_std] = findGaussProd(st_vlog,gvlog(velInds)*hc(hcInds(1)),0,sigP);
    st_var         = st_std^2;
    
    for ts = 1:numel(vStim2Delta)
        
        testV = vStim2Delta(ts)*vStim1(velInds);
        
        % Define test pars
        test_vlog          = log(1 + (testV)/0.3);
        [test_mu,test_std] = findGaussProd(test_vlog,gvlog(velInds)*hc(hcInds(2)),0,sigP);
        test_var           = test_std^2;

        % Evaluate p(test > standard) for this test velocity
        v_range = linspace(...
            min([test_mu st_mu]) - 3*max([test_var st_var]),...
            max([test_mu st_mu]) + 3*max([test_var st_var]),...
            100);
        
        pSDT(ts) = ...
            sum( ...
            (    normpdf(v_range,test_mu,sqrt(test_var))...
            /sum(normpdf(v_range,test_mu,sqrt(test_var)))).*...
            cumsum(normpdf(v_range,st_mu,sqrt(st_var))...
            /sum(       normpdf(v_range,st_mu,sqrt(st_var)))) );
    end
end