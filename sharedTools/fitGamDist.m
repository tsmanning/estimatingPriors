function [aFit,bFit] = fitGamDist(data,numPriors,numSampPairs,numTrCnts)

    % Fit gamma distribution to a set of data with only the a parameter
    
    numPars = numPriors*numTrCnts;
    
    %% Find best vals of a1-ax & b via MLE
    opts    = optimset('display','off','useparallel',0,'tolx',1e-13,...
                   'maxfunevals',1e4,'largescale', 'off');
    
    % gampdf: data, a parameter, b parameter
    % all inputs must be same number of elements: #priors x #trcounts x
    % #pairs
    lossFxn = @(pars) sum(-log( gampdf(data(:),...
                                       repelem(pars(1:end-numTrCnts),numSampPairs),...
                                       repelem(pars(end-(numTrCnts-1):end),numSampPairs*numPriors) ...
                                       )),'omitnan');
    
    parInit = [ones(numPars,1);0.02*ones(numTrCnts,1)];
    
    parFit  = fminunc(lossFxn,parInit,opts);

    aFit = parFit(1:end-numTrCnts);
    bFit = parFit(end-(numTrCnts-1):end);
    
end