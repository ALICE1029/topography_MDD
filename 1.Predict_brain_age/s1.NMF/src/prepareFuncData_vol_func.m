function sbjData = prepareFuncData_vol_func(path,filename,list_cell,maskName)
% fileList -- path of functional files (.nii): each row for one image
% maskName -- path of the brain mask nii
%
% output
% sbjData contains a cell structure sbjData (cell(sbjNum,1)), in which each
% entry sbjData{i} is a matrix of size t x v (# of time points by # of voxels)
%
maskNii = load_untouch_nii(maskName);
maskMat = maskNii.img~=0;



sbjNum = length(list_cell);
sbjData = cell(sbjNum,1);
disp('Read images...');

for si=1:sbjNum%%
    
    
    
    file_name=strcat(path,list_cell(si),filename);
    sbjNii = load_untouch_nii(cell2mat(file_name));
    vxNum = sum(maskMat(:)~=0);
    tNum = size(sbjNii.img,4);
    dataMat = zeros(tNum,vxNum,'single');
    for ti=1:tNum
        tImg = sbjNii.img(:,:,:,ti);
        dataMat(ti,:) = tImg(maskMat);
    end
    sbjData{si} = dataMat;
end


