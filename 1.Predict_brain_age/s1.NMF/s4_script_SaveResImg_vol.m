clc;
clear;
%% soft atlas
PATH_tmp =  '/HeLabData2/cxpang/CBIRD_FOR_DIDA_NGSR/result/init_robust/18/init.mat';
outDir ='/HeLabData2/cxpang/CBIRD_FOR_DIDA_NGSR/result/init_robust/18';
maskName ='/HeLabData2/cxpang/fmri/bnu/mask/GMMask_3mm.nii';
saveFig = 0;
refNiiName ='/HeLabData2/cxpang/fmri/bnu/mask/GMMask_3mm.nii';
func_saveVolRes2Nii(PATH_tmp,maskName,outDir,saveFig,refNiiName);
%% hard atlas
res = load(PATH_tmp);
initV = res.initV;
[m,clust] = max(initV,[],2);
maskName ='/HeLabData2/cxpang/fmri/bnu/mask/GMMask_3mm.nii';
saveFig = 0;
refNiiName ='/HeLabData2/cxpang/fmri/bnu/mask/GMMask_3mm.nii.gz';

maskNii = load_untouch_nii(maskName);

if ~exist(outDir,'dir')
    mkdir(outDir);
end

kNii = maskNii;
kNii.img(maskNii.img~=0) = clust;
outName = [outDir,filesep,'hard.nii.gz'];
save_untouch_nii(kNii,outName);