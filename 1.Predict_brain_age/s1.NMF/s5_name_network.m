%%
%determine all
load('/HeLabData2/cxpang/CBIRD_FOR_DIDA_NGSR/result/init_robust/18/init.mat')
[m,p] = max(initV,[],2);
load bnulist
%% calculate the correlation matrix
con_sub=cell(length(list_cell),1);

maskName ='/HeLabData2/cxpang/fmri/bnu/mask/GMMask_3mm.nii';
maskNii = load_untouch_nii(maskName);
maskMat = maskNii.img~=0;

path='/HeLabData2/cxpang/CBIRD_FOR_DIDA_NGSR/FunImgARWSDCFB/N';
filename='/bfcdswrarest.nii';

sbjNum = length(list_cell);
sbjData = cell(sbjNum,1);
disp('Read images...');

for i=1:sbjNum%%
    i
    file_name=strcat(path,list_cell(i),filename);
    sbjNii = load_untouch_nii(cell2mat(file_name));
    vxNum = sum(maskMat(:)~=0);
    tNum = size(sbjNii.img,4);
    dataMat = zeros(tNum,vxNum,'single');
    for ti=1:tNum
        tImg = sbjNii.img(:,:,:,ti);
        dataMat(ti,:) = tImg(maskMat);
    end
     cor=corr(dataMat) ;
     %save(strcat('/HeLabData2/cxpang/fmri/bnu/silhouse/',list_cell{i},'cor'),'cor','-v7.3')
    con_sub{i}=zeros(18,18);
   % load(strcat('/HeLabData2/cxpang/CBIRD_FOR_DIDA_NGSR/result/silhouse/',list_cell{i},'cor.mat'))
    
    con=zeros(18,18);
    num=zeros(18,18);
    for s=1:18
        s
        pos=find(p ==s);
        pos_no=find(p~=s);
        
        con(s,s)=1;
        for j=1:length(pos)
            for k=1:length(pos_no)
                con(s,p(pos_no(k)))=con(s,p(pos_no(k)))+cor(pos(j),pos_no(k));
                num(s,p(pos_no(k)))=num(s,p(pos_no(k)))+1;
            end
        end
        for t=1:18
            if(num(s,t)~=0)
                con(s,t)=con(s,t)/num(s,t);
            end
        end
    end
    con_sub{i}=con;
    clear cor
end
con_ave=zeros(18,18);
for i=1:length(list_cell)
    con_ave=con_ave+con_sub{i};
end
con_ave=con_ave/length(list_cell);
%%
%method2:overlay with yeo
maskNii = load_untouch_nii('/HeLabData2/cxpang/fmri/bnu/yeo/yeo_18_new.nii.gz');
%:soft atlas
yeo_18_mask=maskNii.img(maskNii.img~=0);
pos=cell(18,1);
for i=1:18
    pos{i}=find(yeo_18_mask==i);
end

load_my=zeros(18,18);
for i=1:18
    tmp=initV(:,i);
    for j=1:18
        load_my(i,j)=sum(tmp(pos{j}));
    end
    
end
for i=1:18
    max_load=max(load_my(i,:));
    load_pos(i)=find(load_my(i,:)==max_load);
end
%hard atlas
[m,p] = max(initV,[],2);
for i=1:18
    pos_new=find(p==i);
    for j=1:18
        dice_new(i,j)=mydice(pos_new,pos{j});
    end
end
for i=1:18
    max_dice(i)=max(dice_new(i,:));
    dice_pos(i)=find(dice_new(i,:)==max_dice(i));
end