%%
%determine all
load('/HeLabData2/cxpang/CBIRD_FOR_DIDA_NGSR/result/init_robust/18/init.mat')
[m,p] = max(initV,[],2);
load bnulist

%overlay with yeo
maskNii = load_untouch_nii('/HeLabData2/cxpang/fmri/bnu/yeo/yeo_18_new.nii.gz');
%soft atlas
yeo_18_mask=maskNii.img(maskNii.img~=0);
pos=cell(18,1);
for i=1:18
    pos{i}=find(yeo_18_mask==i);
end

load_yeo=zeros(18,18);
for i=1:18
    tmp=initV(:,i);
    for j=1:18
        load_yeo(i,j)=sum(tmp(pos{j}));
    end
    
end
for i=1:18
    max_load=max(load_yeo(i,:));
    index_load(i)=find(load_yeo(i,:)==max_load);
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
    index_dice(i)=find(dice_new(i,:)==max_dice(i));
end
