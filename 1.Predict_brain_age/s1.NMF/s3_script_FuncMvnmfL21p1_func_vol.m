

%addpath(genpath('path-to-{CollaborativeBrainDecomposition}'));

%dida1
maskFile ='/HeLabData2/cxpang/fmri/mask/GMMask_3mm.nii';
prepDataFile = '/HeLabData2/cxpang/fmri/test_CreatePrepData.mat';%neighbor
outDir = '/HeLabData2/cxpang/DIDA/NGSR/stage1_old';%for dida dataset
load didas1
resId = 'fmri_vol';
path='/HeLabData2/cxpang/DIDA/NGSR/stage1_old/'; 
%path='/HeLabData2/cxpang/treatment_NGSR/PKU/FunImgARWSDCF/';%for PKU dataset

file_name='/func.nii';%for dida dataset
%file_name='/xbcds4D.nii';%for PKU dataset
K = 18;%to do 
initName = strcat('/HeLabData2/cxpang/CBIRD_FOR_DIDA_NGSR/result/init_robust/',num2str(K),'/init.mat');%k to do
alphaS21 = 2;% same with li
alphaL = 10;
vxI = 1;
spaR = 1;
ard = 0;%
eta = 1;
iterNum = 30;% same with cui
calcGrp = 0;
parforOn = 0;
deployFuncMvnmfL21p1_func_vol(path,file_name,list_cell,maskFile,prepDataFile,outDir,resId,initName,K,alphaS21,alphaL,vxI,spaR,ard,eta,iterNum,calcGrp,parforOn);

%dida2
outDir = '/HeLabData2/cxpang/DIDA/NGSR/stage2_old';%for dida dataset
%outDir = '/HeLabData2/cxpang/DIDA/GSR/'
load dida2
%load didas1
resId = 'fmri_vol';
path='/HeLabData2/cxpang/DIDA/NGSR/stage2_old/'; 
%path='/HeLabData2/cxpang/treatment_NGSR/PKU/FunImgARWSDCF/';%for PKU dataset

file_name='/func.nii';%for dida dataset
%file_name='/xbcds4D.nii';%for PKU dataset
K = 18;%to do 
initName = strcat('/HeLabData2/cxpang/CBIRD_FOR_DIDA_NGSR/result/init_robust/',num2str(K),'/init.mat');%k to do
alphaS21 = 2;% same with li
alphaL = 10;
vxI = 1;
spaR = 1;
ard = 0;%
eta = 1;
iterNum = 30;% same with cui
calcGrp = 0;
parforOn = 0;
sbj_num=1410;
deployFuncMvnmfL21p1_func_vol(path,file_name,sbj_num,maskFile,prepDataFile,outDir,resId,initName,K,alphaS21,alphaL,vxI,spaR,ard,eta,iterNum,calcGrp,parforOn);
