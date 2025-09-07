clc;
clear;
data_folder='/HeLabData2/cxpang'
% for volumetric data
% add subcortical
maskFile = [data_folder '/fmri/mask/GMMask_3mm.nii'];
maskNii = load_untouch_nii(maskFile);

gNb = createPrepData('volumetric', maskNii.img, 1);

% save gNb into file for later use
prepDataName = [data_folder '/fmri/test_CreatePrepData.mat'];
save(prepDataName, 'gNb');
