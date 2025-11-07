clc;
clear;
%% write loading matrix
%addpath(genpath('path-to-{CollaborativeBrainDecomposition}'));
for i=18
  candidateLstFile = strcat(('/HeLabData2/cxpang/CBIRD_FOR_DIDA_NGSR/result/init/fmri_volboostrap_comp'),num2str(i),'/init.txt'); % the path list for every bootstrap init.mat
  outDir =strcat( '/HeLabData2/cxpang/CBIRD_FOR_DIDA_NGSR/result/init_robust/',num2str(i));
   initV = selRobustInit(candidateLstFile,i,outDir,i);

end
