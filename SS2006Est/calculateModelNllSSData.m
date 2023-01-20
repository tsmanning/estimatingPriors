function nll = calculateModelNllSSData(parVec,datStr,vStim1,cStim1,cStim2,vStim2)

%% Likelihood of data given current model pars

% avlog = parVec(1:numel(vStim1));
% gvlog = parVec(numel(vStim1) + 1:2*numel(vStim1));
% hc    = parVec(2*numel(vStim1) + 1:end);
avlog = parVec(1:numel(vStim1));
gvlog = 0.2*ones(1,numel(vStim1));
hc    = parVec(numel(vStim1) + 1:end);

ptvs = calcPtvsSSData(cStim1,vStim1,cStim2,vStim2,avlog,gvlog,hc);

expCnts = [numel(cStim1) numel(vStim1) numel(cStim2)];

% Don't like likelihood go to zero
minlikli = 1e-12;

ind = 1;
nllChunk = nan(prod(expCnts),1);

% For each reference contrast
for rc = 1:expCnts(1)
    
    % For each reference speed
    for rv = 1:expCnts(2)
        
        % For each test contrast
        for tc = 1:expCnts(3)
            
            thisPfxn      = ptvs{rc,rv,tc}';
            thisTrCnt     = datStr(rc,rv,tc).numTr;
            thisChCnt     = datStr(rc,rv,tc).numCh;
            
            if ~isempty(thisPfxn)
                % Find likelihood of this chunk of data
                this_like = binopdf(thisChCnt,thisTrCnt,thisPfxn);
            
                nllChunk(ind) = sum(log(max(minlikli,this_like)));
            else
                % Dataset is missing comparison of 7.5% with 50% contrast
                nllChunk(ind) = 0; 
            end
                
            ind = ind + 1;
            
        end
        
    end
    
end

nll = -sum(nllChunk);

end