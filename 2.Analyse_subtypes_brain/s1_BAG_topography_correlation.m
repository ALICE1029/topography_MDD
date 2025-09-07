%%
load dida
voxel_num=45892;
network_num=18;
sbj_num=2170;
method='NGSR';
load ('/HeLabData2/cxpang/DIDA/using_mat/glmvariable_remove65.mat')%age,group,mfd,sex
load('/home/cxpang/matlab/code/8.brain_age/index_all_remove65sex.mat','Idx')
addpath(genpath('/home/cxpang/matlab/code/8.brain_age'));
addpath(genpath('/home/cxpang/matlab/Collaborative_Brain_Decomposition-master/lib/NIfTI_20140122'))
%%
load('/home/cxpang/matlab/code/8.brain_age/adjusted_mdd_all_gap_remove65sex.mat','adjusted_mdd_all_gap')
hc=find(group==1);
mdd=find(group==2);
sub_gap=zeros(sbj_num,1);
sub_gap(mdd)=adjusted_mdd_all_gap;
sub1_index=mdd(find(Idx==1));
sub2_index=mdd(find(Idx==2));
sex_sub1=(sex(sub1_index));
mfd_sub1=mfd(sub1_index);
age_sub1=age(sub1_index);
sex_sub2=(sex(sub2_index));
mfd_sub2=mfd(sub2_index);
age_sub2=age(sub2_index);
%% lm
load(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/y_all_remove65.mat'),'y_all')
loading_hc=cell(18,1);
for i=1:18
    loading_hc{i}=y_all{i}(hc,:);
end
loading_sub1=cell(18,1);
for i=1:18
    loading_sub1{i}=y_all{i}(sub1_index,:);
end
loading_sub2=cell(18,1);
for i=1:18
    loading_sub2{i}=y_all{i}(sub2_index,:);
end
save('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/loading_subtype2_remove65.mat','loading_sub1','loading_sub2');
%% sub1
[n_voxel,n_sub] = size(loading_sub1{1}');
T_value = zeros(network_num,n_voxel);
p_value = zeros(network_num,n_voxel);
res_map = cell(network_num,1);
% regress out sub_gap
des = [mfd_sub1];
sub1_gap=sub_gap(sub1_index);
sub1_gap=abs(sub1_gap);
[sub1_gap, b, stats] = regress_out(sub1_gap, des);
des = [sub1_gap,sex_sub1,age_sub1,mfd_sub1];
for i = 1:network_num
    i
    res_map{i} = zeros(n_sub,n_voxel);
    for j = 1:n_voxel
        stat_result = regstats(loading_sub1{i}(:,j),des,'linear',{'tstat','r'});
        T_value(i,j) = stat_result.tstat.t(2);
        p_value(i,j) = stat_result.tstat.pval(2);
        res_map{i}(:,j) = stat_result.r;
        % dfe(i,j)=stat_result.tstat.dfe;
    end
end
dfe=stat_result.tstat.dfe;
save(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/betweengroup_difference/corr_sub2_1_remove65sex.mat'),'T_value','p_value','res_map','dfe');
T_value (isnan (T_value))=0;
p_value (isnan (p_value))=1;
% get z value
Z_value=zeros(network_num,voxel_num);
for i=1:network_num
    i
    tmp=T_value(i,:);
    Z_value(i,:) = spm_t2z(tmp,dfe);
end
vox = [3 3 3];% voxel size
voxel_p = 0.001;
cluster_p = 0.05;
tail = 2;
zthrd = norminv(1 - voxel_p/tail);
Z_vol=cell(network_num,1);
t_vol=cell(network_num,1);
z_grf=cell(network_num,1);
t_grf=cell(network_num,1);
cluster_size=zeros(network_num,1);
my_size=cell(18,1);
coor_x=cell(18,1);
coor_y=cell(18,1);
coor_z=cell(18,1);
cd=cell(18,1);
cd_min=cell(18,1);
br=cell(18,1);
Z=cell(18,1);
maskName ='/home/cxpang/matlab/code/3.Betweengroup_Difference/using_mat/rbrodmann.nii';
brodnii = load_untouch_nii(maskName);
brod=brodnii.img;
for i =1:18% 1:network_num
    i
    hdr_mask = spm_vol('/HeLabData2/cxpang/fmri/mask/GMMask_3mm.nii');% 2024-8-3 mask
    vol_mask = spm_read_vols(hdr_mask);
    ind = find(vol_mask);
    R_volume = zeros([hdr_mask.dim,n_sub]);
    for j = 1:n_sub
        tmp_vol = zeros(hdr_mask.dim);
        tmp=res_map{i}(j,:);
        tmp_vol(ind) = tmp;%r value
        R_volume(:,:,:,j) = tmp_vol;
    end
    j=1;
    if(i<10)
        gunzip(strcat ('/HeLabData2/cxpang/CBIRD_FOR_DIDA_NGSR/result/init_robust/18/icn_00',num2str(i),'.nii.gz'))
        hdr_mask2 = spm_vol(strcat ('/HeLabData2/cxpang/CBIRD_FOR_DIDA_NGSR/result/init_robust/18/icn_00',num2str(i),'.nii'));
    else
        gunzip(strcat ('/HeLabData2/cxpang/CBIRD_FOR_DIDA_NGSR/result/init_robust/18/icn_0',num2str(i),'.nii.gz'))
        hdr_mask2 = spm_vol(strcat ('/HeLabData2/cxpang/CBIRD_FOR_DIDA_NGSR/result/init_robust/18/icn_0',num2str(i),'.nii'));
    end
    vol_mask2 = spm_read_vols(hdr_mask2);
    [cluster_size(i),dlh,fwhm] = x_GRF(R_volume,dfe,vol_mask2,vox,voxel_p,cluster_p,tail);
    Z_vol{i} = zeros(hdr_mask.dim);
    tmp=Z_value(i,:);
    Z_vol{i}(ind) = tmp;% z value nii
    Z_vol{i} (Z_vol{i}  < zthrd & Z_vol {i} > -zthrd) = 0;
    t_vol{i}=zeros(hdr_mask.dim);
    tmp_t=T_value(i,:);
    t_vol{i}(ind) = tmp_t;
    t_vol{i} (Z_vol{i}  < zthrd & Z_vol {i} > -zthrd) = 0;
    [L, num] = bwlabeln(Z_vol{i},26);
    n = 0;
    for x = 1:num
        theCurrentCluster = L == x;     
        if length(find(theCurrentCluster)) <= cluster_size(i)
            n = n + 1;
            Z_vol{i} (logical(theCurrentCluster)) = 0;
            t_vol{i} (logical(theCurrentCluster)) = 0;
        else
            my_size{i}(j)=length(find(theCurrentCluster));
            tmp2=Z_vol{i}(find(theCurrentCluster));
            max_positive = max(tmp2(tmp2 > 0));
            min_negative = min(tmp2(tmp2 < 0));
            if ~isempty(max_positive)
                Z{i}(j) = max_positive;
            else
                Z{i}(j) = min_negative;
            end
            [vol_mask,corr]=spm_read_vols(hdr_mask);
            pos= find(Z_vol{i}==Z{i}(j));
            coor_x{i}(j)=corr(1,pos);
            coor_y{i}(j)=corr(2,pos);
            coor_z{i}(j)=corr(3,pos);
            br{i}(j)=brod(pos);
            j=j+1;
        end
    end
    z_grf{i}=Z_vol{i}(ind);
    t_grf{i}=t_vol{i}(ind);
    outDir = strcat ('/HeLabData2/cxpang/DIDA/',method,'/result_xy/corr_sub2_1_remove65sexmask');
    saveFig = 0;
    maskName ='/HeLabData2/cxpang//fmri/mask/GMMask_3mm.nii';
    maskNii = load_untouch_nii(maskName);
    if ~exist(outDir,'dir')
        mkdir(outDir);
    end
    kNii = maskNii;
    kNii.img(maskNii.img~=0) = t_grf{i};
    outName = [outDir,filesep,strcat('corr_',num2str(i),'.nii.gz')];
    save_untouch_nii(kNii,outName);
end
save(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/betweengroup_difference/corr_sub2_1_sta_remove65sexmask.mat'),'coor_x','coor_y','coor_z','my_size','Z','br')
save(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/betweengroup_difference/corr_sub2_1_grf_remove65sexmask.mat'),'z_grf','t_grf')
%% result
load(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/betweengroup_difference/corr_sub2_1_grf_remove65sexmask.mat'))
load(strcat('/HeLabData2/cxpang/DIDA/using_mat/grf_mask_',method,'.mat'))
num=zeros(2,network_num);
z_map=z_grf;
pos_sig=cell(network_num,1);
pos_nosig=cell(network_num,1);
for i=1:network_num
    pos_pn=find(z_map{i}>0);%positive voxel
    num(1,i)=length(pos_pn);
    pos_n=find(z_map{i}<0);
    num(2,i)=length(pos_n);
end
save(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/betweengroup_difference/cres_sub2_1_remove65sexmask.mat'),'num')
%% sub2
[n_voxel,n_sub] = size(loading_sub2{1}');
T_value = zeros(network_num,n_voxel);
p_value = zeros(network_num,n_voxel);
res_map = cell(network_num,1);
% regress out sub_gap
des = [mfd_sub2];
sub2_gap=sub_gap(sub2_index);
[sub2_gap, b, stats] = regress_out(sub2_gap, des);
des = [sub2_gap,sex_sub2,age_sub2,mfd_sub2];
for i = 1:network_num
    i
    res_map{i} = zeros(n_sub,n_voxel);
    for j = 1:n_voxel
        stat_result = regstats(loading_sub2{i}(:,j),des,'linear',{'tstat','r'});
        T_value(i,j) = stat_result.tstat.t(2);
        p_value(i,j) = stat_result.tstat.pval(2);
        res_map{i}(:,j) = stat_result.r;
        % dfe(i,j)=stat_result.tstat.dfe;
    end
end


dfe=stat_result.tstat.dfe;
load(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/betweengroup_difference/corr_sub2_2_remove65sex.mat'));

T_value (isnan (T_value))=0;
p_value (isnan (p_value))=1;
% get z value
Z_value=zeros(network_num,voxel_num);
for i=1:network_num
    i
    tmp=T_value(i,:);
    Z_value(i,:) = spm_t2z(tmp,dfe);
end

vox = [3 3 3];% voxel size
voxel_p = 0.001;
cluster_p = 0.05;
tail = 2;
zthrd = norminv(1 - voxel_p/tail);
Z_vol=cell(network_num,1);
t_vol=cell(network_num,1);
z_grf=cell(network_num,1);
t_grf=cell(network_num,1);
cluster_size=zeros(network_num,1);
my_size=cell(18,1);
coor_x=cell(18,1);
coor_y=cell(18,1);
coor_z=cell(18,1);
cd=cell(18,1);
cd_min=cell(18,1);
br=cell(18,1);
Z=cell(18,1);
maskName ='/home/cxpang/matlab/code/3.Betweengroup_Difference/using_mat/rbrodmann.nii';
brodnii = load_untouch_nii(maskName);
brod=brodnii.img;
for i = 1:network_num
    i
    
    hdr_mask = spm_vol('/HeLabData2/cxpang/fmri/mask/GMMask_3mm.nii');
    vol_mask = spm_read_vols(hdr_mask);
    ind = find(vol_mask);
    R_volume = zeros([hdr_mask.dim,n_sub]);
    for j = 1:n_sub
        % pos=find(mask(i,:)==0);
        tmp_vol = zeros(hdr_mask.dim);
        tmp=res_map{i}(j,:);
        %tmp(pos)=[];
        tmp_vol(ind) = tmp;%r value
        R_volume(:,:,:,j) = tmp_vol;
    end% r value nii
    j=1;
    if(i<10)
        gunzip(strcat ('/HeLabData2/cxpang/CBIRD_FOR_DIDA_NGSR/result/init_robust/18/icn_00',num2str(i),'.nii.gz'))
        hdr_mask2 = spm_vol(strcat ('/HeLabData2/cxpang/CBIRD_FOR_DIDA_NGSR/result/init_robust/18/icn_00',num2str(i),'.nii'));
    else
        gunzip(strcat ('/HeLabData2/cxpang/CBIRD_FOR_DIDA_NGSR/result/init_robust/18/icn_0',num2str(i),'.nii.gz'))
        hdr_mask2 = spm_vol(strcat ('/HeLabData2/cxpang/CBIRD_FOR_DIDA_NGSR/result/init_robust/18/icn_0',num2str(i),'.nii'));
    end
    vol_mask2 = spm_read_vols(hdr_mask2);
    [cluster_size(i),dlh,fwhm] = x_GRF(R_volume,dfe,vol_mask2,vox,voxel_p,cluster_p,tail);
    Z_vol{i} = zeros(hdr_mask.dim);
    
    tmp=Z_value(i,:);
    %tmp(pos)=[];
    Z_vol{i}(ind) = tmp;% z value nii
    Z_vol{i} (Z_vol{i}  < zthrd & Z_vol {i} > -zthrd) = 0;
    
    t_vol{i}=zeros(hdr_mask.dim);
    tmp_t=T_value(i,:);
    %tmp_t(pos)=[];
    t_vol{i}(ind) = tmp_t;
    t_vol{i} (Z_vol{i}  < zthrd & Z_vol {i} > -zthrd) = 0;
    
    [L, num] = bwlabeln(Z_vol{i},26);
    n = 0;
    for x = 1:num
        theCurrentCluster = L == x;
        
        if length(find(theCurrentCluster)) <= cluster_size(i)
            n = n + 1;
            Z_vol{i} (logical(theCurrentCluster)) = 0;
            t_vol{i} (logical(theCurrentCluster)) = 0;
        else
            my_size{i}(j)=length(find(theCurrentCluster));
            tmp2=Z_vol{i}(find(theCurrentCluster));
            %  Z{i}(j)= max(abs(Z_vol{i}(find(theCurrentCluster)))) ;
            max_positive = max(tmp2(tmp2 > 0));
            min_negative = min(tmp2(tmp2 < 0));
            
            if ~isempty(max_positive)
                Z{i}(j) = max_positive;
            else
                Z{i}(j) = min_negative;
            end
            [vol_mask,corr]=spm_read_vols(hdr_mask);
            pos= find(Z_vol{i}==Z{i}(j));
            coor_x{i}(j)=corr(1,pos);
            coor_y{i}(j)=corr(2,pos);
            coor_z{i}(j)=corr(3,pos);
            br{i}(j)=brod(pos);
            j=j+1;
        end
    end
    z_grf{i}=Z_vol{i}(ind);
    t_grf{i}=t_vol{i}(ind);
    outDir = strcat ('/HeLabData2/cxpang/DIDA/',method,'/result_xy/corr_sub2_2_remove65sexmask');
    saveFig = 0;
    maskName ='/HeLabData2/cxpang//fmri/mask/GMMask_3mm.nii';
    maskNii = load_untouch_nii(maskName);
    if ~exist(outDir,'dir')
        mkdir(outDir);
    end
    kNii = maskNii;
    kNii.img(maskNii.img~=0) = t_grf{i};
    outName = [outDir,filesep,strcat('corr_',num2str(i),'.nii.gz')];
    save_untouch_nii(kNii,outName);
end
save(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/betweengroup_difference/corr_sub2_2_sta_remove65sexmask.mat'),'coor_x','coor_y','coor_z','my_size','Z','br')
save(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/betweengroup_difference/corr_sub2_2_grf_remove65sexmask.mat'),'z_grf','t_grf')
%% result
load(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/betweengroup_difference/corr_sub2_2_grf_remove65sex.mat'))
num=zeros(2,network_num);
%percentage=zeros(2,network_num);
z_map=z_grf;
pos_sig=cell(network_num,1);
pos_nosig=cell(network_num,1);
for i=1:network_num
    %  pos=find(mask(i,:)~=0);
    pos_pn=find(z_map{i}>0);%positive voxel
    num(1,i)=length(pos_pn);
    pos_n=find(z_map{i}<0);
    num(2,i)=length(pos_n);
    %     percentage(1,i)=  num(1,i)/(length(pos));
    %     percentage(2,i)=num(2,i)/(length(pos));
end
save(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/betweengroup_difference/cres_sub2_2_remove65sexmask.mat'),'num')
%% combine sub 1
load(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/betweengroup_difference/corr_sub2_1_grf_remove65sexmask.mat'),'z_grf')
t_combine=zeros(voxel_num,1);
t=zeros(network_num,voxel_num);
for i = 1:network_num
    i
    outDir = strcat ('/HeLabData2/cxpang/DIDA/',method,'/result_xy/corr_sub2_1_remove65sexmask');
    outName = [outDir,filesep,strcat('corr_',num2str(i),'.nii.gz')];
    nii = load_untouch_nii(outName);
    grf=nii.img;
    maskName ='/HeLabData2/cxpang//fmri/mask/GMMask_3mm.nii';
    maskNii = load_untouch_nii(maskName);
    t(i,:)=grf(maskNii.img~=0);
    
end
for i=1:network_num
    tmp=t(i,:);
    tmp(find(tmp<0))=0;
    t_pos(i,:)=tmp;
    
    tmp=t(i,:);
    tmp(find(tmp>0))=0;
    t_neg(i,:)=tmp;
end

outDir = strcat ('/HeLabData2/cxpang/DIDA/',method,'/result_xy/corr_sub2_1_remove65sexmask');
saveFig = 0;
maskName ='/HeLabData2/cxpang//fmri/mask/GMMask_3mm.nii';
maskNii = load_untouch_nii(maskName);
if ~exist(outDir,'dir')
    mkdir(outDir);
end
kNii = maskNii;
kNii.img(maskNii.img~=0) = mean(t_pos,1);
outName = [outDir,filesep,strcat('t_combine_pos_grf.nii.gz')];
save_untouch_nii(kNii,outName);

kNii = maskNii;
kNii.img(maskNii.img~=0) = mean(t_neg,1);
outName = [outDir,filesep,strcat('t_combine_neg_grf.nii.gz')];
save_untouch_nii(kNii,outName);
%% combine sub 2
load(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/betweengroup_difference/corr_sub2_2_grf_remove65sexmask.mat'),'t_grf')
t_combine=zeros(voxel_num,1);
t=zeros(network_num,voxel_num);
for i = 1:network_num
    i
    outDir = strcat ('/HeLabData2/cxpang/DIDA/',method,'/result_xy/corr_sub2_2_remove65sexmask');
    outName = [outDir,filesep,strcat('corr_',num2str(i),'.nii.gz')];
    nii = load_untouch_nii(outName);
    grf=nii.img;
    maskName ='/HeLabData2/cxpang//fmri/mask/GMMask_3mm.nii';
    maskNii = load_untouch_nii(maskName);
    t(i,:)=grf(maskNii.img~=0);
end
for i=1:network_num
    tmp=t(i,:);
    tmp(find(tmp<0))=0;
    t_pos(i,:)=tmp;
    
    tmp=t(i,:);
    tmp(find(tmp>0))=0;
    t_neg(i,:)=tmp;
end

outDir = strcat ('/HeLabData2/cxpang/DIDA/',method,'/result_xy/corr_sub2_2_remove65sexmask');
saveFig = 0;
maskName ='/HeLabData2/cxpang//fmri/mask/GMMask_3mm.nii';
maskNii = load_untouch_nii(maskName);
if ~exist(outDir,'dir')
    mkdir(outDir);
end
kNii = maskNii;
kNii.img(maskNii.img~=0) = mean(t_pos,1);
outName = [outDir,filesep,strcat('t_combine_pos_grf.nii.gz')];
save_untouch_nii(kNii,outName);

kNii = maskNii;
kNii.img(maskNii.img~=0) = mean(t_neg,1);
outName = [outDir,filesep,strcat('t_combine_neg_grf.nii.gz')];
save_untouch_nii(kNii,outName);