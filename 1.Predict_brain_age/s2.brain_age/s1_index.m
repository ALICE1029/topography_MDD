load dida
voxel_num=45892;
network_num=18;
method='NGSR';
addpath(genpath('/home/cxpang/matlab/combat'))
%% hc index_generate
load ('/HeLabData2/cxpang/DIDA/using_mat/glmvariable.mat')
load ('/HeLabData2/cxpang/DIDA/using_mat/site.mat')
new_index=find(age>=65);
age(new_index)=[];
mfd(new_index)=[];
sex(new_index)=[];
group(new_index)=[];
site(new_index)=[];
save('/HeLabData2/cxpang/DIDA/using_mat/glmvariable_remove65.mat','age','sex','mfd','group')
save('/HeLabData2/cxpang/DIDA/using_mat/site_remove65.mat','site')
%% load loading matrix
x=zeros(voxel_num*network_num,length(list_cell));
%prepare for combat
for i=1:length(list_cell)
    i
    if(i<=1431)
        name=cell2mat(strcat('/HeLabData2/cxpang/DIDA/',method,'/stage1/fmri_vol_sbj_s1_comp',num2str(network_num),'_alphaS21_1_alphaL10_vxInfo1_ard0_eta1/',list_cell(i),'final_UV.mat'));
        load(name);
        x(:,i)=V(:);%v:45892*18
        clear V
        clear U
    else
        name=cell2mat(strcat('/HeLabData2/cxpang/DIDA/',method,'/stage2/fmri_vol_sbj_s1_comp',num2str(network_num),'_alphaS21_1_alphaL10_vxInfo1_ard0_eta1/',list_cell(i),'final_UV.mat'));
        load(name);
        x(:,i)=V(:);
        clear V
        clear U
    end
end
save(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/before_combat/regress_variable.mat'),'x')
%% feature 1:loading matrix 
load(strcat('/HeLabData2/cxpang/DIDA/',method,'/result/before_combat/regress_variable.mat'))%(voxel*network)*sub
pos=find(all(x==0,2));% because combat can't accept variable that is same,so find the position that everyone's loading is zero
save('/HeLabData2/cxpang/DIDA/using_mat/pos_zero.mat','pos')
x (pos,:) = [];
x(:,new_index)=[];
sex_d=dummyvar(sex);%dummy
group_d=dummyvar(group);
mod=[age sex_d(:,2) group_d(:,2) mfd];
data_harmonized = combat(x,(site),mod,1);%combat
x_tmp=zeros(voxel_num*network_num,length(site));
for i=1:voxel_num*network_num
    pos_tmp(i)=i;
end
pos_new=setdiff(pos_tmp,pos);
save('/HeLabData2/cxpang/DIDA/using_mat/pos_nozero.mat','pos_new')
x_tmp(pos_new,:)=data_harmonized;
clear data_harmonized
save(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/after_combat/variable_combat_remove65.mat'),'x_tmp')
%% another mod:18 cell (sbj_num*45892)
load(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/after_combat/variable_combat_remove65.mat'),'x_tmp')
sbj_num=size(x_tmp,2);
y=x_tmp;
clear x_tmp
y_net=zeros(sbj_num,voxel_num);
y_all=cell(network_num,1);
for j=1:network_num
    for i=1:sbj_num
        v=y(:,i);
        v=reshape(v,[voxel_num,network_num]);
        y_net(i,:)=v(:,j);
        clear v
    end
    y_all{j,1}=y_net;
end
clear y
clear y_net
save(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/y_all_remove65.mat'),'y_all')
%% feature :loading matrix hc
load(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/after_combat/variable_combat_remove65.mat'),'x_tmp')
hc=find(group==1);
site=site(hc);
age=age(hc);
mfd=mfd(hc);
sex=sex(hc);
save(strcat('/home/cxpang/matlab/code//8.brain_age/using_mat/covariate_hc_remove65.mat'),'site','age','mfd','sex')
%load loading matrix for all subjects
x_tmp=x_tmp(:,hc);
x_tmp (pos,:) = [];
loading=x_tmp';%sub*()
save(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/brain_age/loading_hc_remove65.mat'),'loading')

%% feature :loading matrix mdd
load('/HeLabData2/cxpang/DIDA/using_mat/glmvariable_remove65.mat','age','sex','mfd','group')
mdd=find(group==2);
age=age(mdd);
mfd=mfd(mdd);
sex=sex(mdd);
save(strcat('/home/cxpang/matlab/code//8.brain_age/using_mat/covariate_mdd_remove65.mat'),'age','mfd','sex')
%load loading matrix for all subjects
load(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/after_combat/variable_combat_remove65.mat'),'x_tmp')
x_tmp=x_tmp(:,mdd);
x_tmp (pos,:) = [];
loading=x_tmp';
save(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/brain_age/loading_mdd_remove65.mat'),'loading')

