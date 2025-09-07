clc;
clear;
data_folder='/HeLabData2/cxpang/';
maskFile = [data_folder '/fmri/mask/GMMask_3mm.nii'];
prepDataFile = [data_folder 'test_CreatePrepData.mat'];
outDir = [data_folder '/CBIRD_FOR_DIDA_NGSR/result/init'];

path=[data_folder '/CBIRD_FOR_DIDA_NGSR/FunImgARWSDCFB/N'];
filename='/bfcdswrarest.nii';

spaR = 1;
vxI = 1;
ard = 0;
iterNum = 2000;
tNum = 190;
alpha = 2;
beta = 10;
resId = 'fmri_vol';

load bnulist %the group atlas used data list
sbjNum = length(list_cell);
%random
random_number=50;
every_time_sub=117;%146*0.8
%this step is for not produce same random number at same time point in
%parallel
rand_num=zeros(random_number,every_time_sub);
for i=1:random_number
    i
     a=randperm(sbjNum);
     b=a(1:every_time_sub);
     rand_num(i,:)=b;
     fprintf('  sbj%d',b);
     fprintf('\n')
end
rand_num=sort(rand_num,2);

deployFuncInit_vol(random_number,rand_num,every_time_sub,path,filename,list_cell,maskFile,prepDataFile,outDir,spaR,vxI,ard,iterNum,tNum,alpha,beta,resId);%,rand_num);