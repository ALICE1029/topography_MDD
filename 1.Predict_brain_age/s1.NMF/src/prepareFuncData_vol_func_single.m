function sbjData = prepareFuncData_vol_func_single(path,file_name,maskMat,si)
    filename=strcat(path,si,file_name);%%change
    sbjNii = load_untouch_nii((filename));
    vxNum = sum(maskMat(:)~=0);
    tNum = size(sbjNii.img,4);
    dataMat = zeros(tNum,vxNum,'single');
    for ti=1:tNum
        tImg = sbjNii.img(:,:,:,ti);
        dataMat(ti,:) = tImg(find(maskMat));
    end
    %save(cell2mat(strcat('/HeLabData2/cxpang/fmri/bnu/func_mat/func_',list_cell(si),'.mat')),'dataMat','-v7.3')
    %save(cell2mat(strcat('/HeLabData2/cxpang/fmri/DIDA/func_mat/func_',list_cell(si),'.mat')),'dataMat','-v7.3')
    %load(cell2mat(strcat('/HeLabData2/cxpang/fmri/bnu/func_mat/func_',list_cell(si),'.mat')))  
    sbjData= dataMat;
    rowsToKeep = ~all(isnan(sbjData), 2);
    sbjData = sbjData(rowsToKeep, :);
end

