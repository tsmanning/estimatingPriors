clear all

%% Collect and sort fits across resampling runs

% Define analyses to combine
% lapseRates  = [0 0.05];
lapseRates  = [0];
domain      = 'log';

% model       = 'ss';
model       = 'mg';
% model       = 'mgfix';
% model       = 'mgfixneg';
% model       = 'zeroMu';
% model       = 'zeroMuWeq';

% priorF      = 'expon';
% priorF      = 'splitExp';
priorF      = 'mixgauss';

switch priorF
    case 'expon'
%         priorList      = 1:5;
%         ds             = '2021-05-02';
        priorList      = 1:3;
        ds             = '2021-08-17';
    case 'splitExp'
        priorList      = {[-3 -1.5],[-1.5 -3]};
        ds             = '2021-04-12';
    case 'mixgauss' 
        priorList      = 1:3;
%         priorList      = 1:5;
%         ds             = '2021-04-29';
        ds             = '2021-09-01';
%         ds             = '2021-05-01';
end

numResamps = 25;

% numSamps = [5 10 20];
numSamps = 20;

% subDir = '';
subDir = ['1Sep2021_mgNV',filesep];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tbt         = 1;
binscal     = 0;

% method      = 'numeric';
method      = 'closed';

switch model
    case 'ss'
        interp      = 2;
    otherwise
        %%%% another hack to fix
        if strcmp(priorF,'expon')
            interp      = 2;
        else
            interp      = 0;
        end
end

interp = 0;

if tbt
    if interp == 1
        prefix = 'tbtInterp';
    elseif interp == 2
        prefix = 'tbtNN';
    else
        prefix = 'tbt';
    end
else
    if binscal
        prefix = 'bsc';
    else
        prefix = '';
    end
end

%% Define directories
splPath = regexp(which('dataCollector'),filesep,'split');
workDir = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep,'SimData',filesep,method,filesep,subDir];

ind = 1;

numPriors = numel(priorList);

for i = 1:numel(lapseRates)
    for j = 1:numel(priorList)
        
        switch priorF
                
            case 'splitExp'
%                 pwTemp = priorList{j};
%                 pString = [num2str(pwTemp(1)),num2str(pwTemp(2))];
%                 dataSets{ind} = ['gt_',model,'_',domain,'_',priorF,'_',...
%                                  pString,'_LR_',num2str(lapseRates(i))];
                error('havent worked this out yet');
                
            otherwise
                dataSets{ind} = ['gt_',domain,'_',priorF,'_pInd_',num2str(j),...
                                 '_',num2str(numPriors),'_LR_',num2str(lapseRates(i)),...
                                 '_',model,'_',ds];
                
        end
        
        ind = ind + 1;
    end
end

% Loop over datasets
for ds = 1:numel(dataSets)
    
    sDir = [workDir,dataSets{ds},filesep];
    
    % Loop over different num trial repeats
    for ns = 1:numel(numSamps)
        
        % Loop over different samples for this exp condition
        for rs = 1:numResamps
            
            % No suffix for first run
            suffix = ['_',num2str(rs)];
            
            thisFN = [sDir,'parEsts_',model,prefix,'_jitter_n',num2str(numSamps(ns)),suffix,'.mat'];
            
            if rs == 1
                % Only load GT stuff, for first sample
                load(thisFN);
            else
                switch model
                    case 'ss'
                        load(thisFN,'avlogF','gvlogF','hcF','initConds','nllF','ptvsF')
                    otherwise
                        load(thisFN,'sigPF','gvlogF','hcF','wPF','muPF','initConds','nllF','ptvsF')
                end
            end
            
            %%% hack for now to get rid of very small nll error runs
            nllF(nllF<1) = nan;
            
            % Find best fitting run
            [nllFt(rs),minInd] = min(nllF);
            
            % Select a different run if an obviously bad fit occurs (i.e.
            % alternative minima with increasing sigma with contrast)
            if strcmp(model,'ss')
                if hcF(minInd,end)>hcF(minInd,1)
                    
                    excInds = hcF(:,end) > hcF(:,1);
                    nllTemp = nllF;
                    nllTemp(excInds) = nan;
                    [nllFt(rs),minInd] = min(nllTemp);
                    
                end
            end

            % Put best fit into resampling mats
            switch model
                case 'ss'
                    avFt(rs,:) = avlogF(minInd,:);
                otherwise
                    spFt(rs,:) = sigPF(minInd,:);
                    muFt(rs,:) = muPF(minInd,:);
                    wFt(rs,:)  = wPF(minInd,:);
            end
            
            gvFt(rs,:) = gvlogF(minInd,:);
            hcFt(rs,:) = hcF(minInd,:);
            ict(rs,:)  = initConds(minInd,:);
            ptvsFt{rs} = ptvsF{minInd};
            oFN{rs}    = thisFN;
            nllGTt     = nllGT;
            
        end
        
        clear avlogF sigPF gvlogF hcF wPF muPF initConds nllF ptvsF
        
        % Reassign so these mats now carry fits from different SAMPLES, not
        % from different INITIALIZATIONS from same data sample
        switch model
            case 'ss'
                avlogF = avFt;
            otherwise
                sigPF  = spFt;
                muPF   = muFt;
                wPF    = wFt;
        end
        
        gvlogF    = gvFt;
        hcF       = hcFt;
        initConds = ict;
        ptvsF     = ptvsFt;
        nllF      = nllFt;
        nllGT     = nllGTt;
        
        if tbt
            if interp == 1
                saveFName = [sDir,'parEsts_',model,'_n',num2str(numSamps(ns)),'_collectTBTinterp.mat'];
            elseif interp == 2
                saveFName = [sDir,'parEsts_',model,'_n',num2str(numSamps(ns)),'_collectTBTNN.mat'];
            else
                saveFName = [sDir,'parEsts_',model,'_n',num2str(numSamps(ns)),'_collectTBT.mat'];
            end
        else
            if binscal
                saveFName = [sDir,'parEsts_',model,'_n',num2str(numSamps(ns)),'_collectbiScl.mat'];
            else
                saveFName = [sDir,'parEsts_',model,'_n',num2str(numSamps(ns)),'_collect.mat'];
            end
        end
        
        switch model
            case 'ss'
                if strcmp(priorF,'expon')
                    save(saveFName,'avlog','avlogF','gvlog','gvlogF','hc','hcF',...
                        'initConds','nllF','nllGT','numTrials','ptvs_data','ptvsF');
                else
                    save(saveFName,'avlog','avlogF','gvlog','gvlogF','hc','hcF',...
                        'sigP','muP','wP','initConds','nllF','nllGT','numTrials','ptvs_data','ptvsF');
                end
                
            otherwise
                if strcmp(priorF,'mixgauss')
                    save(saveFName,'sigP','sigPF','gvlog','gvlogF','hc','hcF',...
                        'wP','wPF','muP','muPF','initConds','nllF','nllGT',...
                        'numTrials','ptvs_data','ptvsF');
                else
                    %%%% ONLY FOR RUN ON 2021-04-12
                    avlog = avlog(end-5:end);
                    
                    save(saveFName,'avlog','sigPF','gvlog','gvlogF','hc','hcF',...
                        'wPF','muPF','initConds','nllF','nllGT',...
                        'numTrials','ptvs_data','ptvsF');
                end
                
        end
        
    end
    
end
