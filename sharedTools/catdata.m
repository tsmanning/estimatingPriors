function [runStruct,JSDpriorsTrCnts,priorPairs] = catdata(matfile1,matfile2)

datDir = '/media/tyler/Data/MATLAB/cooperLab/2-Modeling_Simulations/BayesModelComp/SimData/closed/';

load([datDir,matfile1]);
runStructA       = runStruct;
JSDpriorsTrCntsA = JSDpriorsTrCnts;

load([datDir,matfile2]);
runStructB       = runStruct;
JSDpriorsTrCntsB = JSDpriorsTrCnts;

runStruct       = cat(3,runStructA,runStructB);
JSDpriorsTrCnts = cat(1,JSDpriorsTrCntsA,JSDpriorsTrCntsB);

end