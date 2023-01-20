clear all

% Collect and sort fits across resampling runs

workDir = which('dataCollector');
% workDir = [workDir(1:end-numel('dataCollector.m')),'ssGauss/'];
workDir = workDir(1:end-numel('dataCollector.m'));

dataSets = {'gt_log_gauss_0.8_LR_0','gt_log_gauss_1_LR_0','gt_log_gauss_1.2_LR_0','gt_log_gauss_1.4_LR_0','gt_log_gauss_1.6_LR_0',...
            'gt_log_gauss_0.8_LR_0.05','gt_log_gauss_1_LR_0.05','gt_log_gauss_1.2_LR_0.05','gt_log_gauss_1.4_LR_0.05','gt_log_gauss_1.6_LR_0.05'};

numResamps = 3;

numSamps = [20 80 250];

% Loop over datasets
for ds = 1:numel(dataSets)
    
    sDir = [workDir,dataSets{ds},filesep];
    
    for ns = 1:numel(numSamps)
        
        for rs = 1:numResamps
            
            if rs == 1
                suffix = '';
            else
                suffix = ['_',num2str(rs)];
            end
            
            thisFN = [sDir,'parEsts_gs_jitter_n',num2str(numSamps(ns)),suffix,'.mat'];
            
            if rs == 1
                load(thisFN);
            else
                load(thisFN,'sigPF','gvlogF','hcF','initConds','nllF','ptvsF')
            end
            
            % Find best fitting run
            [nllFt(rs),minInd] = min(nllF);
            
            % Put best fit into resampling mats
            spFt(rs,:) = sigPF(minInd,:);
            gvFt(rs,:) = gvlogF(minInd,:);
            hcFt(rs,:) = hcF(minInd,:);
            ict(rs,:)  = initConds(minInd,:);
            ptvsFt{rs} = ptvsF{minInd};
            oFN{rs}    = thisFN;
            
        end
        
        clear sigPF gvlogF hcF initConds nllF ptvsF
        
        % Reassign and save
        sigPF    = spFt;
        gvlogF    = gvFt;
        hcF       = hcFt;
        initConds = ict;
        ptvsF     = ptvsFt;
        nllF      = nllFt;
        
        saveFName = [sDir,'parEsts_n',num2str(numSamps(ns)),'_collectG.mat'];
        
        save(saveFName,'sigP','sigPF','gvlog','gvlogF','hc','hcF',...
                       'initConds','nllF','numTrials','ptvs_data','ptvsF');
        
    end
    
end
