% Plot relative likelihoods of data

splPath = regexp(which('compareLogProbData'),filesep,'split');
fDir    = [filesep,fullfile(splPath{1:numel(splPath)-2}),filesep];
sDir    = [fDir,'SS2006Data',filesep];

ds = '2021-07-27';
rsMethod = 'triplet';

mog1 = load([sDir,'fitData_',ds,'_1comps_nn_',rsMethod,'_mog']);
mog2 = load([sDir,'fitData_',ds,'_2comps_nn_',rsMethod,'_mog']);
mog3 = load([sDir,'fitData_',ds,'_3comps_nn_',rsMethod,'_mog']);
load([sDir,'fitData_',ds,'_',rsMethod,'_pw']);

load([sDir,'s1']);
load([sDir,'s2']);

m1Nlls1 = mog1.mogDat.bestNllMoG_s1;
m1Nlls2 = mog1.mogDat.bestNllMoG_s2;
m2Nlls1 = mog2.mogDat.bestNllMoG_s1;
m2Nlls2 = mog2.mogDat.bestNllMoG_s2;
m3Nlls1 = mog3.mogDat.bestNllMoG_s1;
m3Nlls2 = mog3.mogDat.bestNllMoG_s2;
pNlls1 = pwDat.bestNllPW_s1;
pNlls2 = pwDat.bestNllPW_s2;

nllMat = [pNlls1 m1Nlls1 m2Nlls1 m3Nlls1;...
          pNlls2 m1Nlls2 m2Nlls2 m3Nlls2];
      
[coinFlipS1,cumGaussS1] = getModelFitBounds(s1);
[coinFlipS2,cumGaussS2] = getModelFitBounds(s2);     
      
nllMat = [(nllMat(1,:)-coinFlipS1)/(cumGaussS1-coinFlipS1);...
          (nllMat(2,:)-coinFlipS2)/(cumGaussS2-coinFlipS2)];

f7 = figure;
f7.Position = [2050 380 1100 600];

bH = bar(nllMat);

colorMat = [0.8*ones(1,3);...
            0.6*ones(1,3);...
            0.4*ones(1,3);...
            0.2*ones(1,3)];

legData = {'Piecewise','MoG (1 Comp)','MoG (2 Comp)','MoG (3 Comp)'};
        
for ii = 1:size(nllMat,2)
   bH(ii).FaceColor = colorMat(ii,:);
   bH(ii).EdgeColor = colorMat(ii,:);
   
   a = annotation(f7,'textbox',[0.75 0.84-0.05*(ii-1) 0.2 0.1],'string',legData{ii},'edgecolor','none');
   a.Color = colorMat(ii,:);
   a.FontSize = 20;
end

set(gca,'fontsize',20,'ylim',[0 1],'ytick',linspace(0,1,6),...
    'ytickLabel',{'Coin Flip','','','','','Cumulative Gaussian'},...
    'plotboxaspectratio',[1 1 1]);
ylabel('Log-probability of data');
xlabel('Subject');

