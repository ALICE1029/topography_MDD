addpath(genpath('/home/cxpang/matlab/combat'))
%%
voxel_num=45892;
network_num=18;
method='NGSR';
load('/home/cxpang/matlab/code/8.brain_age/using_mat/index_all_remove65sex.mat')
load ('/HeLabData2/cxpang/DIDA/using_mat/glmvariable_remove65.mat')%age,group,mfd,sex
load ('/HeLabData2/cxpang/DIDA/using_mat/site_remove65.mat')%age,group,mfd,sex
hc=find(group==1);
mdd=find(group==2);
sub1_index=mdd(find(Idx==1));
sub2_index=mdd(find(Idx==2));
%% 
mdd_hamd=importdata('/home/cxpang/matlab/code/4.HAMD_correlation/using_mat/mdd_hamd_num_sub.txt');% index of patients have hamd score
sub1_hamd=intersect(mdd_hamd,sub1_index);
load ('/HeLabData2/cxpang/DIDA/using_mat/glmvariable_remove65.mat')%age,group,mfd,sex
age=age(sub1_hamd);
sex=sex(sub1_hamd);
mfd=mfd(sub1_hamd);
site=site(sub1_hamd);
save('/home/cxpang/matlab/code/8.brain_age/using_mat/covariate_hamd_sub_sub1_remove65.mat','age','sex','mfd','site')
for i=1:length(sub1_hamd)
sub1_hamd_index(i)=find(mdd_hamd==sub1_hamd(i));
end
load(strcat('/home/cxpang/matlab/code/4.HAMD_correlation/using_mat/hamd_sub.mat'))
hamd=hamd(sub1_hamd_index,:);
save('/home/cxpang/matlab/code/8.brain_age/using_mat/hamd_sub_sub1_remove65.mat','hamd')

sub2_hamd=intersect(mdd_hamd,sub2_index);
load ('/HeLabData2/cxpang/DIDA/using_mat/glmvariable_remove65.mat')%age,group,mfd,sex
load ('/HeLabData2/cxpang/DIDA/using_mat/site_remove65.mat')%age,group,mfd,sex
age=age(sub2_hamd);
sex=sex(sub2_hamd);
mfd=mfd(sub2_hamd);
site=site(sub2_hamd);
save('/home/cxpang/matlab/code/8.brain_age/using_mat/covariate_hamd_sub_sub2_remove65.mat','age','sex','mfd','site')
for i=1:length(sub2_hamd)
sub2_hamd_index(i)=find(mdd_hamd==sub2_hamd(i));
end
load(strcat('/home/cxpang/matlab/code/4.HAMD_correlation/using_mat/hamd_sub.mat'))
hamd=hamd(sub2_hamd_index,:);
save('/home/cxpang/matlab/code/8.brain_age/using_mat/hamd_sub_sub2_remove65.mat','hamd')

%% feature 1:loading matrix 
load ('/HeLabData2/cxpang/DIDA/using_mat/glmvariable.mat')
new_index=find(age>=65);
load(strcat('/HeLabData2/cxpang/DIDA/',method,'/result/before_combat/regress_variable.mat'))
x(:,new_index)=[];
x=x(:,mdd_hamd);
load(strcat('/home/cxpang/matlab/code/4.HAMD_correlation/using_mat/hamd_sub.mat'))
mdd_hamd=importdata('/home/cxpang/matlab/code/4.HAMD_correlation/using_mat/mdd_hamd_num_sub.txt');
sbj_num=length(mdd_hamd);
load ('/HeLabData2/cxpang/DIDA/using_mat/glmvariable_remove65.mat')%age,group,mfd,sex
load ('/HeLabData2/cxpang/DIDA/using_mat/site_remove65.mat')%
site=site(mdd_hamd);
pos=find(all(x==0,2));% because combat can't accept variable that is same,so find the position that everyone's loading is zero
save(strcat('/HeLabData2/cxpang/DIDA/',method,'/hamd/loading_pos_sub.mat'),'pos')% to further analyse loading weight
sex(find(sex==-1))=2;
sex_d=dummyvar(sex);%dummy
sex=sex_d(:,2);
mod=[age(mdd_hamd) hamd sex(mdd_hamd) mfd(mdd_hamd)];% combat: hamd as protect variable
x (pos,:) = [];
data_harmonized = combat(x,(site),mod,1);%combat
loading=data_harmonized';
save(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/hamd/loading_sub_remove65.mat'),'loading')

loading_sub1=loading(sub1_hamd_index,:);
loading_sub2=loading(sub2_hamd_index,:);

save('/HeLabData2/cxpang/DIDA/NGSR/result_xy/hamd/loading_sub_sub1_remove65.mat','loading_sub1')
save('/HeLabData2/cxpang/DIDA/NGSR/result_xy/hamd/loading_sub_sub2_remove65.mat','loading_sub2')

%% step2_pls_sub_topography.py
%% see the correlation result
fold=10;
for  k=1
    draw=cell(2,1);
    for j=0:fold-1
        load(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/hamd/cca/AtlasLoading/loading_sub_sub1_select0.001/Time_0/Fold_',num2str(j),'comp_uv.mat'))
        tmp1=u(:,k);
        tmp2=v(:,k);
        draw{1}=[draw{1};tmp1];
        draw{2}=[draw{2};tmp2];
    end
    my_corr(k)=corr( draw{1}, draw{2});
end
%% hamd subscale permut
corr_permut_all=zeros(1000,1);
bug=[];
for k=1:1000
    k
    draw=cell(2,1);
    ResultantFolder = [strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/hamd/cca/AtlasLoading/loading_sub_sub1_select0.001/Permutation/',num2str(k)),'/Time_0/'];
 
        for  s=1
             draw=cell(2,1);
            for i=0:9%to avoid the condition that predict value same,correlation nan
                ResultantFolder = [strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/hamd/cca/AtlasLoading/loading_sub_sub1_select0.001/Permutation/',num2str(k)),'/Time_0/'];
                load(strcat(ResultantFolder,'Fold_',num2str(i),'comp_uv.mat'))
                tmp1=u(:,s);
                tmp2=v(:,s);
                draw{1}=[draw{1};tmp1];
                draw{2}=[ draw{2};tmp2];
            end
            
            corr_permut_all(k,s)=corr(draw{1},draw{2});
        end
 
end
for k=1
p(k)=length(find(corr_permut_all(:,k)>=my_corr(k)))/1000;
end
h=histfit(corr_permut_all);
xline(my_corr,'r','LineWidth',2)
set(gca,'Box','off');
