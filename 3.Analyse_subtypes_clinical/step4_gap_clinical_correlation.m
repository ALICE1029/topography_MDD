%%
load dida
voxel_num=45892;
network_num=18;
sbj_num=2170;
method='NGSR';
load ('/HeLabData2/cxpang/DIDA/using_mat/glmvariable_remove65.mat')%age,group,mfd,sex
load ('/HeLabData2/cxpang/DIDA/using_mat/site_remove65.mat')%age,group,mfd,sexv
load('/home/cxpang/matlab/code/8.brain_age/index_all_remove65sex.mat','Idx')
addpath(genpath('/home/cxpang/matlab/code/8.brain_age'));
%% 
load('/home/cxpang/matlab/code/8.brain_age/adjusted_mdd_all_gap_remove65sex.mat','adjusted_mdd_all_gap')
hc=find(group==1);
mdd=find(group==2);
sub_gap=zeros(sbj_num,1);
sub_gap(mdd)=adjusted_mdd_all_gap;
sub1_index=mdd(find(Idx==1));
sub2_index=mdd(find(Idx==2));
%%  hamd
mdd_hamd=importdata('/home/cxpang/matlab/code/4.HAMD_correlation/using_mat/mdd_hamd_num_hamd.txt');
load(strcat('/home/cxpang/matlab/code/4.HAMD_correlation/using_mat/hamd.mat'))
load ('/HeLabData2/cxpang/DIDA/using_mat/glmvariable_remove65.mat')%age,group,mfd,sex
hamd_all=[mdd_hamd,hamd];
sub_hamd_1_index=intersect(mdd_hamd,sub1_index);
sub_hamd_2_index=intersect(mdd_hamd,sub2_index);
sex_sub1=sex(sub_hamd_1_index);
sex_sub2=sex(sub_hamd_2_index);
age_sub1=age(sub_hamd_1_index);
age_sub2=age(sub_hamd_2_index);
mfd_sub1=mfd(sub_hamd_1_index);
mfd_sub2=mfd(sub_hamd_2_index);
for i=1:length(sub_hamd_1_index)
    hamd_sub1(i)=hamd_all(find(hamd_all(:,1)==sub_hamd_1_index(i)),2);
end
for i=1:length(sub_hamd_2_index)
    hamd_sub2(i)=hamd_all(find(hamd_all(:,1)==sub_hamd_2_index(i)),2);
end
gap_1=sub_gap(sub_hamd_1_index);
gap_2=sub_gap(sub_hamd_2_index);

hamd_sub1=hamd_sub1';
hamd_sub2=hamd_sub2';
%sub 1
des = [age_sub1,sex_sub1];
[hamd_sub1, b, stats] = regress_out(hamd_sub1, des);
des = [mfd_sub1];
[gap_1, b, stats] = regress_out(gap_1, des);
[c1,p1]=corr(gap_1,hamd_sub1,'Type','Pearson')
%sub 2
des = [age_sub2,sex_sub2];
[hamd_sub2, b, stats] = regress_out(hamd_sub2, des);
des = [mfd_sub2];
[gap_2, b, stats] = regress_out(gap_2, des);
[c2,p2]=corr(gap_2,hamd_sub2,'Type','Pearson')
%% 4. subscale
mdd_hamd=importdata('/home/cxpang/matlab/code/4.HAMD_correlation/using_mat/mdd_hamd_num_sub.txt');
load(strcat('/home/cxpang/matlab/code/4.HAMD_correlation/using_mat/hamd_sub.mat'))
%try site z
mdd_hamd=importdata('/home/cxpang/matlab/code/4.HAMD_correlation/using_mat/mdd_hamd_num_sub.txt');
load(strcat('/home/cxpang/matlab/code/4.HAMD_correlation/using_mat/hamd_sub.mat'))
load ('/HeLabData2/cxpang/DIDA/using_mat/glmvariable_remove65.mat')%age,group,mfd,sexv

p1=zeros(17,1);
p2=zeros(17,1);
for j=1:17
    j
    hamd_all=[mdd_hamd,hamd(:,j)];
    sub_hamd_1_index=intersect(sub1_index,mdd_hamd);
    sub_hamd_2_index=intersect(mdd_hamd,sub2_index);
    
    sex_sub1=sex(sub_hamd_1_index);
    sex_sub2=sex(sub_hamd_2_index);
    
    age_sub1=age(sub_hamd_1_index);
    age_sub2=age(sub_hamd_2_index);
    
    mfd_sub1=mfd(sub_hamd_1_index);
    mfd_sub2=mfd(sub_hamd_2_index);
     
    site_sub1=site(sub_hamd_1_index);
    site_sub2=site(sub_hamd_2_index);
    hamd_sub_sub1=zeros(length(sub_hamd_1_index),1);
    for i=1:length(sub_hamd_1_index)
        hamd_sub_sub1(i)=hamd_all(find(hamd_all(:,1)==sub_hamd_1_index(i)),2);
    end
    hamd_sub_sub2=zeros(length(sub_hamd_2_index),1);
    for i=1:length(sub_hamd_2_index)
        hamd_sub_sub2(i)=hamd_all(find(hamd_all(:,1)==sub_hamd_2_index(i)),2);
    end
    %
    
    des1 = [mfd_sub1];
    sub1_gap=sub_gap(sub_hamd_1_index);
     
     [sub1_gap, b, stats] = regress_out(sub1_gap, des1);
    %
    des=[age_sub1,sex_sub1];
    [hamd_sub_sub1, b, stats] = regress_out(hamd_sub_sub1, des);
    
    [c1(j),p1(j)]=corr(sub1_gap,hamd_sub_sub1,'Type','Spearman');
    
    des2 = [mfd_sub2];
    sub2_gap=sub_gap(sub_hamd_2_index);
    [sub2_gap, b, stats] = regress_out(sub2_gap, des2);
%     
    des=[age_sub2,sex_sub2 ];
    [hamd_sub_sub2, b, stats] = regress_out(hamd_sub_sub2, des);
    
    [c2(j),p2(j)]=corr(sub2_gap,hamd_sub_sub2,'Type','Spearman');
    
end
