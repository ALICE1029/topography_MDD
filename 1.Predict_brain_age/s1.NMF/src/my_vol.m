function dataMat = my_vol( filename ,maskMat,id)
%MY_VOL Summary of this function goes here
%   Detailed explanation goes here
    sbjNii = load_untouch_nii(filename);
    vxNum = sum(maskMat(:)~=0);
    tNum = size(sbjNii.img,4);
    dataMat = zeros(tNum,vxNum,'single');
    for ti=1:tNum
        tImg = sbjNii.img(:,:,:,ti);
        dataMat(ti,:) = tImg(maskMat);
    end
    save(strcat('/HeLabData2/cxpang/fmri/bnu/func_mat/func_',(id),'.mat'),'dataMat','-v7.3')
end

