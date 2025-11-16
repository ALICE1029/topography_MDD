load('/home/cxpang/matlab/code/8.brain_age/adjusted_mdd_all_gap_remove65sex.mat','adjusted_mdd_all_gap')
Idx=[];
Idx(find(adjusted_mdd_all_gap<0))=1;
Idx(find(adjusted_mdd_all_gap>0))=2;
save('/home/cxpang/matlab/code/8.brain_age/index_all_remove65sex.mat','Idx')
%%
addpath(genpath('/home/cxpang/matlab/code/8.brain_age'));
addpath(genpath('/home/cxpang/matlab/Collaborative_Brain_Decomposition-master/lib/NIfTI_20140122'))
load dida
voxel_num=45892;
network_num=18;
method='NGSR';
load ('/HeLabData2/cxpang/DIDA/using_mat/glmvariable_remove65.mat')%age,group,mfd,sex
load('/home/cxpang/matlab/code/8.brain_age/index_all_remove65sex.mat','Idx')
hc=find(group==1);
mdd=find(group==2);
sub1_index=mdd(find(Idx==1));
sub2_index=mdd(find(Idx==2));

sex=sex-1;
hc_sex=sex(hc);
hc_mfd=mfd(hc);
hc_age=age(hc);

mdd1_sex=(sex(sub1_index));
mdd1_mfd=mfd(sub1_index);
mdd1_age=age(sub1_index);
mdd2_sex=(sex(sub2_index));
mdd2_mfd=mfd(sub2_index);
mdd2_age=age(sub2_index);
%% BAG-related topography trajectory
t1=zeros(network_num,voxel_num);
for i = 1:network_num
    i
    outDir = strcat ('/HeLabData2/cxpang/DIDA/',method,'/result_xy/corr_sub2_1_remove65sexmask');
    outName = [outDir,filesep,strcat('corr_',num2str(i),'.nii.gz')];
    nii = load_untouch_nii(outName);
    grf=nii.img;
    maskName ='/HeLabData2/cxpang//fmri/mask/GMMask_3mm.nii';
    maskNii = load_untouch_nii(maskName);
    t1(i,:)=grf(maskNii.img~=0);
end
t1=t1';
pos_sub1=find(t1(:)>0);
neg_sub1=find(t1(:)<0);
% t1: (45892*18)
load(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/after_combat/variable_combat_remove65.mat'),'x_tmp')
lo=x_tmp;
lo_pos=lo(pos_sub1,:);
lo_neg=lo(neg_sub1,:);

weight_p_pos=pca(lo_pos');
weight_pca_pos = lo_pos'* weight_p_pos(:,:);
weight_pca_pos=weight_pca_pos(:,1);

weight_p_neg=pca(lo_neg');
weight_pca_neg = lo_neg'* weight_p_neg(:,:);
weight_pca_neg=weight_pca_neg(:,1);

%% normalized to 1-2
if(min(weight_pca_pos)<0)
    weight_pca_pos=weight_pca_pos+abs(min(weight_pca_pos))+0.1;
else
    weight_pca_pos=weight_pca_pos;
end
if(min(weight_pca_neg)<0)
    weight_pca_neg=weight_pca_neg+abs(min(weight_pca_neg))+0.1;
else
    weight_pca_neg=weight_pca_neg;
end
hc_weight=weight_pca_pos(hc);
sub1_weight=weight_pca_pos(sub1_index);

prepare_vars_for_norm_models('all','t1psighc',1, hc_weight,hc_sex, hc_age, hc_mfd)
prepare_vars_for_norm_models('all','t1psigsub1',1, sub1_weight,mdd1_sex, mdd1_age, mdd1_mfd)

hc_weight=weight_pca_neg(hc);
sub1_weight=weight_pca_neg(sub1_index);

prepare_vars_for_norm_models('all','t1nsighc',1, hc_weight,hc_sex, hc_age, hc_mfd)
prepare_vars_for_norm_models('all','t1nsigsub1',1, sub1_weight,mdd1_sex, mdd1_age, mdd1_mfd)

%% permut sub1
for i=1:6
    m1_age_group{i}=find(mdd1_age>=i*10&mdd1_age<i*10+10);
    hc_age_group{i}=find(hc_age>=i*10&hc_age<i*10+10);
end
m1_sex_group{1}=find(mdd1_sex==1);
m1_sex_group{2}=find(mdd1_sex==0);
hc_sex_group{1}=find(hc_sex==1);
hc_sex_group{2}=find(hc_sex==0);

for j=1:1000
    j
    hc_index_permut=[];
    sub1_index_permut=[];
    for i=1:6
        hc_1=intersect(hc_age_group{i},hc_sex_group{1});
        hc_2=intersect(hc_age_group{i},hc_sex_group{2});
        
        m1_1=intersect(m1_age_group{i},m1_sex_group{1});
        m1_2=intersect(m1_age_group{i},m1_sex_group{2});
        tmp_hc=hc( hc_1);
        tmp_1=sub1_index(m1_1);
        tmp  =[tmp_hc;tmp_1];
        randd=datasample(tmp, length(hc_1),'Replace',false);
        hc_index_permut=[hc_index_permut;randd];
        sub1_index_permut=[sub1_index_permut;setdiff(tmp,randd)];
        
        tmp_hc=hc( hc_2);
        tmp_1=sub1_index( m1_2);
        tmp  =[tmp_hc;tmp_1];
        randd=datasample(tmp, length(hc_2),'Replace',false);
        hc_index_permut=[hc_index_permut;randd];
        sub1_index_permut=[sub1_index_permut;setdiff(tmp,randd)];
        
    end
    hc_index_permut_all(j,:)=hc_index_permut;
    sub1_index_permut_all(j,:)=sub1_index_permut;
end
for i=1:1000
    %
    hc_sex=sex(hc_index_permut_all(i,:));
    hc_age=age(hc_index_permut_all(i,:));
    hc_mfd=mfd(hc_index_permut_all(i,:));
    mdd1_sex=sex(sub1_index_permut_all(i,:));
    mdd1_age=age(sub1_index_permut_all(i,:));
    mdd1_mfd=mfd(sub1_index_permut_all(i,:));
    
    hc_weight_abs=weight_pca_abs(hc_index_permut_all(i,:));
    sub1_weight_abs=weight_pca_abs(sub1_index_permut_all(i,:));
    vars_permut('all','t1asighc',hc_weight_abs,hc_sex, hc_age, hc_mfd,i)
    vars_permut('all','t1asigsub1',sub1_weight_abs,mdd1_sex, mdd1_age, mdd1_mfd,i)
end
%% sub 2
t2=zeros(network_num,voxel_num);
for i = 1:network_num
    i
    outDir = strcat ('/HeLabData2/cxpang/DIDA/',method,'/result_xy/corr_sub2_2_remove65sexmask');
    outName = [outDir,filesep,strcat('corr_',num2str(i),'.nii.gz')];
    nii = load_untouch_nii(outName);
    grf=nii.img;
    maskName ='/HeLabData2/cxpang//fmri/mask/GMMask_3mm.nii';
    maskNii = load_untouch_nii(maskName);
    t2(i,:)=grf(maskNii.img~=0);
end
t2=t2';

pos_sub2=find(t2(:)>0);
neg_sub2=find(t2(:)<0);

lo=x_tmp;
lo_pos=lo(pos_sub2,:);
lo_neg=lo(neg_sub2,:);

weight_p_pos2=pca(lo_pos');
[coeff, score, latent, ~, explained, ~] = pca(lo_pos');
explained_variance = latent / sum(latent);
cumulative_explained_variance = cumsum(explained_variance);
weight_pca_pos2 = lo_pos'* weight_p_pos2(:,:);
weight_pca_pos2=weight_pca_pos2(:,1);

weight_p_neg2=pca(lo_neg');
[coeff, score, latent, ~, explained, ~] = pca(lo_neg');
explained_variance = latent / sum(latent);
cumulative_explained_variance = cumsum(explained_variance);
weight_pca_neg2 = lo_neg'* weight_p_neg2(:,:);
weight_pca_neg2=weight_pca_neg2(:,1);

%% normalized to 1-2
if(min(weight_pca_pos2)<0)
    weight_pca_pos2=weight_pca_pos2+abs(min(weight_pca_pos2))+0.1;
else
    weight_pca_pos2=weight_pca_pos2;
end
if(min(weight_pca_neg2)<0)
    weight_pca_neg2=weight_pca_pos2+abs(min(weight_pca_neg2))+0.1;
else
    weight_pca_neg2=weight_pca_neg2;
end

hc_weight=weight_pca_pos2(hc);
sub2_weight=weight_pca_pos2(sub2_index);
prepare_vars_for_norm_models('all','t2psighc',1, hc_weight,hc_sex, hc_age, hc_mfd)
prepare_vars_for_norm_models('all','t2psigsub2',1, sub2_weight,mdd2_sex, mdd2_age, mdd2_mfd)

hc_weight=weight_pca_neg2(hc);
sub2_weight=weight_pca_neg2(sub2_index);
prepare_vars_for_norm_models('all','t2nsighc',1, hc_weight,hc_sex, hc_age, hc_mfd)
prepare_vars_for_norm_models('all','t2nsigsub2',1, sub2_weight,mdd2_sex, mdd2_age, mdd2_mfd)

%% permut sub2
hc_index_permut_all=[];
sub2_index_permut_all=[];
for i=1:6
    m2_age_group{i}=find(mdd2_age>=i*10&mdd2_age<i*10+10);
    hc_age_group{i}=find(hc_age>=i*10&hc_age<i*10+10);
end
m2_sex_group{1}=find(mdd2_sex==1);
m2_sex_group{2}=find(mdd2_sex==0);
hc_sex_group{1}=find(hc_sex==1);
hc_sex_group{2}=find(hc_sex==0);
for j=1:1000
    j
    hc_index_permut=[];
    sub2_index_permut=[];
    for i=1:6
        hc_1=intersect(hc_age_group{i},hc_sex_group{1});
        
        hc_2=intersect(hc_age_group{i},hc_sex_group{2});
        
        m2_1=intersect(m2_age_group{i},m2_sex_group{1});
        
        m2_2=intersect(m2_age_group{i},m2_sex_group{2});
        tmp_hc=hc( hc_1);
        tmp_1=sub2_index( m2_1);
        tmp  =[tmp_hc;tmp_1];
        randd=datasample(tmp, length(hc_1),'Replace',false);
        hc_index_permut=[hc_index_permut;randd];
        sub2_index_permut=[sub2_index_permut;setdiff(tmp,randd)];
        
        tmp_hc=hc( hc_2);
        tmp_1=sub2_index( m2_2);
        tmp  =[tmp_hc;tmp_1];
        randd=datasample(tmp, length(hc_2),'Replace',false);
        hc_index_permut=[hc_index_permut;randd];
        sub2_index_permut=[sub2_index_permut;setdiff(tmp,randd)];
        
    end
    hc_index_permut_all(j,:)=hc_index_permut;
    sub2_index_permut_all(j,:)=sub2_index_permut;
    
end
for i=1:1000
    %
    i
    hc_sex=sex(hc_index_permut_all(i,:));
    hc_age=age(hc_index_permut_all(i,:));
    hc_mfd=mfd(hc_index_permut_all(i,:));
    mdd2_sex=sex(sub2_index_permut_all(i,:));
    mdd2_age=age(sub2_index_permut_all(i,:));
    mdd2_mfd=mfd(sub2_index_permut_all(i,:));
    
    hc_weight_abs=weight_pca_abs2(hc_index_permut_all(i,:));
    sub2_weight_abs=weight_pca_abs2(sub2_index_permut_all(i,:));
    
    vars_permut('all','t2asighc',hc_weight_abs,hc_sex, hc_age, hc_mfd,i)
    vars_permut('all','t2asigsub2',sub2_weight_abs,mdd2_sex, mdd2_age, mdd2_mfd,i)
end
%% bootstrap
for i=1:6
    m1_age_group{i}=find(mdd1_age>=i*10&mdd1_age<i*10+10);
    m2_age_group{i}=find(mdd2_age>=i*10&mdd2_age<i*10+10);
    hc_age_group{i}=find(hc_age>=i*10&hc_age<i*10+10);
end
m1_sex_group{1}=find(mdd1_sex==1);
m1_sex_group{2}=find(mdd1_sex==0);
m2_sex_group{1}=find(mdd2_sex==1);
m2_sex_group{2}=find(mdd2_sex==0);
hc_sex_group{1}=find(hc_sex==1);
hc_sex_group{2}=find(hc_sex==0);

for j=1:1000
    j
    hc_index_boot=[];
    sub1_index_boot=[];
    sub2_index_boot=[];
    for i=1:6
        hc_1=intersect(hc_age_group{i},hc_sex_group{1});
        selectedIndices = randi(length(hc( hc_1)), length(hc( hc_1)),1);%存在重复值
        tmp=hc( hc_1);
        hc_index_boot=[hc_index_boot;tmp(selectedIndices)];
        hc_2=intersect(hc_age_group{i},hc_sex_group{2});
        selectedIndices = randi(length(hc( hc_2)), length(hc( hc_2)),1);
        tmp=hc( hc_2);
        hc_index_boot=[hc_index_boot;tmp(selectedIndices)];
        
        m1_1=intersect(m1_age_group{i},m1_sex_group{1});
        selectedIndices = randi(length(sub1_index( m1_1)), length(sub1_index( m1_1)),1);
        tmp=sub1_index( m1_1);
        sub1_index_boot=[sub1_index_boot;tmp(selectedIndices)];
        m1_2=intersect(m1_age_group{i},m1_sex_group{2});
        selectedIndices = randi(length(sub1_index( m1_2)), length(sub1_index( m1_2)),1);
        tmp=sub1_index( m1_2);
        sub1_index_boot=[sub1_index_boot;tmp(selectedIndices)];
        
        m2_1=intersect(m2_age_group{i},m2_sex_group{1});
        selectedIndices = randi(length(sub2_index( m2_1)), length(sub2_index( m2_1)),1);
        tmp=sub2_index( m2_1);
        sub2_index_boot=[sub2_index_boot;tmp(selectedIndices)];
        m2_2=intersect(m2_age_group{i},m2_sex_group{2});
        selectedIndices = randi(length(sub2_index( m2_2)), length(sub2_index( m2_2)),1);
        tmp=sub2_index( m2_2);
        sub2_index_boot=[sub2_index_boot;tmp(selectedIndices)];
    end
    hc_index_boot_all(j,:)=hc_index_boot;
    sub1_index_boot_all(j,:)=sub1_index_boot;
    sub2_index_boot_all(j,:)=sub2_index_boot;
end
for i=1:1000
    i
    
    hc_sex=sex(hc_index_boot_all(i,:));
    hc_age=age(hc_index_boot_all(i,:));
    hc_mfd=mfd(hc_index_boot_all(i,:));
    mdd1_sex=sex(sub1_index_boot_all(i,:));
    mdd1_age=age(sub1_index_boot_all(i,:));
    mdd1_mfd=mfd(sub1_index_boot_all(i,:));
    mdd2_sex=sex(sub2_index_boot_all(i,:));
    mdd2_age=age(sub2_index_boot_all(i,:));
    mdd2_mfd=mfd(sub2_index_boot_all(i,:));
    
    hc_weight=weight_pca_abs(hc_index_boot_all(i,:));
    sub1_weight=weight_pca_abs(sub1_index_boot_all(i,:))
    
    vars_boot('all','t1asighc', hc_weight,hc_sex, hc_age, hc_mfd,i)
    vars_boot('all','t1asigsub1', sub1_weight,mdd1_sex, mdd1_age, mdd1_mfd,i)
    
    
    hc_weight=weight_pca_abs2(hc_index_boot_all(i,:));
    sub2_weight=weight_pca_abs2(sub2_index_boot_all(i,:));
    
    vars_boot('all','t2asighc', hc_weight,hc_sex, hc_age, hc_mfd,i)
    vars_boot('all','t2asigsub2', sub2_weight,mdd2_sex, mdd2_age, mdd2_mfd,i)
end
%%
% s2_gamlss_test.R
%% 
clear all
close all

inpath1=('D:\study\sub2\brain_age\variables_for_normative_modeling\all\');
inpath2=('D:\study\sub2\brain_age\variables_for_normative_modeling\output_normative_modeling\all\');
inpath3=('D:\study\sub2\brain_age\variables_for_normative_modeling\output_normative_modeling\boot_abs\');
for i=1:1000
    i
    yfitF=nan(1000,5);
    yfitM=nan(1000,5);
    load([inpath3,...
        'GAMLSS_','t1asighc_',num2str(i),'_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
    faf_age=age;
    load([inpath3,...
        'GAMLSS_','t1asighc_',num2str(i),'_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
    fam_age=age;
    yfit=(yfitF+yfitM)/2;
    fake_age_boot(:,i)=(faf_age+fam_age)/2;
    fit_boot(:,i)=yfit(:,3);
    % sub 1
    yfitF=nan(1000,5);
    yfitM=nan(1000,5);
    load([inpath3,...
        'GAMLSS_','t1asigsub1_',num2str(i),'_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
    faf_age=age;
    load([inpath3,...
        'GAMLSS_','t1asigsub1_',num2str(i),'_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
    fam_age=age;
    yfit=(yfitF+yfitM)/2;
    fake_age_sub1_boot(:,i)=(faf_age+fam_age)/2;
    fit_sub1_boot(:,i)=yfit(:,3);
end
% 95%ci
for i=1:1000
    mean_normal = mean(fit_boot(i,:));
    std_normal = std(fit_boot(i,:));
    
    confidence_level = 0.95; % 95%
    
    df = (length(fit_boot(i,:))-1);
    critical_value = tinv((1 + confidence_level) / 2, df);
    
    ci_normal_hc(i,:) = [mean_normal+ - critical_value * std_normal , ...
        mean_normal + critical_value * std_normal];
    
    
    mean_normal = mean(fit_sub1_boot(i,:));
    std_normal = std(fit_sub1_boot(i,:));
    
    confidence_level = 0.95; % 95%
    
    df = (length(fit_sub1_boot(i,:))-1);
    critical_value = tinv((1 + confidence_level) / 2, df);
    
    ci_normal_sub1(i,:) = [mean_normal+ - critical_value * std_normal , ...
        mean_normal + critical_value * std_normal];
end

yfitM=nan(1000,5);
load([inpath2,...
    'GAMLSS_','t1asighc_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','t1asighc_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit=(yfitF+yfitM)/2;
fake_age=(faf_age+fam_age)/2;

figure; % Ensure it's a new figure
set(gcf, 'Position', [100, 100, 800, 400]);
%gap<0
load([inpath2,...
    'GAMLSS_','t1asigsub1_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','t1asigsub1_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit_neg=(yfitF+yfitM)/2;
fake_age_neg=(faf_age+fam_age)/2;
%
x = [fake_age', fliplr(fake_age')];
y1=ci_normal_hc(:,1);
y2=ci_normal_hc(:,2);
y = [(y1)', fliplr(y2')];
fill(x, y, [0.6 0.6 0.6],'EdgeColor', 'none','facealpha',0.5);
hold on
plot(fake_age,(yfit(:,3)),'k','linewidth',2.5); %plot 50th centile

hold on
x = [fake_age', fliplr(fake_age')];
y = [ci_normal_sub1(:,1)', fliplr(ci_normal_sub1(:,2)')];
fill(x, y, [0.7, 0.7, 1], 'EdgeColor', 'none','facealpha',0.5);
hold on

box off
set(gca, 'FontSize', 14);

%% sub 2
clear all
close all

inpath1=('D:\study\sub2\brain_age\variables_for_normative_modeling\all\');
inpath2=('D:\study\sub2\brain_age\variables_for_normative_modeling\output_normative_modeling\all\');
inpath3=('D:\study\sub2\brain_age\variables_for_normative_modeling\output_normative_modeling\boot_abs\');
for i=1:1000%i boot times
    
    yfitF=nan(1000,5);
    yfitM=nan(1000,5);
    load([inpath3,...
        'GAMLSS_','t2asighc_',num2str(i),'_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
    faf_age=age;
    load([inpath3,...
        'GAMLSS_','t2asighc_',num2str(i),'_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
    fam_age=age;
    yfit=(yfitF+yfitM)/2;
    fake_age_boot(:,i)=(faf_age+fam_age)/2;
    fit_boot(:,i)=yfit(:,3);
    % sub 1
    yfitF=nan(1000,5);
    yfitM=nan(1000,5);
    load([inpath3,...
        'GAMLSS_','t2asigsub2_',num2str(i),'_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
    faf_age=age;
    load([inpath3,...
        'GAMLSS_','t2asigsub2_',num2str(i),'_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
    fam_age=age;
    yfit=(yfitF+yfitM)/2;
    fake_age_sub1_boot(:,i)=(faf_age+fam_age)/2;
    fit_sub1_boot(:,i)=yfit(:,3);
end
% 95%ci
for i=1:1000
    mean_normal = mean(fit_boot(i,:));
    std_normal = std(fit_boot(i,:));
    
    confidence_level = 0.95; % 95%
    
    df = (length(fit_boot(i,:))-1);
    critical_value = tinv((1 + confidence_level) / 2, df);
    
    ci_normal_hc(i,:) = [mean_normal+ - critical_value * std_normal , ...
        mean_normal + critical_value * std_normal];
    
    
    mean_normal = mean(fit_sub1_boot(i,:));
    std_normal = std(fit_sub1_boot(i,:));
    
    confidence_level = 0.95; % 95%
    
    df = (length(fit_sub1_boot(i,:))-1);
    critical_value = tinv((1 + confidence_level) / 2, df);
    
    ci_normal_sub1(i,:) = [mean_normal+ - critical_value * std_normal , ...
        mean_normal + critical_value * std_normal];
end

yfitM=nan(1000,5);
load([inpath2,...
    'GAMLSS_','t2asighc_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','t2asighc_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit=(yfitF+yfitM)/2;
fake_age=(faf_age+fam_age)/2;

% plot(fake_age,(yfit(:,3)),'k','linewidth',2.5); %plot 50th centile
% hold on
%gap<0
load([inpath2,...
    'GAMLSS_','t2asigsub2_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','t2asigsub2_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit_neg=(yfitF+yfitM)/2;
fake_age_neg=(faf_age+fam_age)/2;
%
figure; % Ensure it's a new figure
set(gcf, 'Position', [100, 100, 800, 400]);
x = [fake_age', fliplr(fake_age')];
y = [ci_normal_sub1(:,1)', fliplr(ci_normal_sub1(:,2)')];
fill(x, y, [1, 0.7, 0.7], 'EdgeColor', 'none','facealpha',0.5); % 灰色阴影
hold on
x = [fake_age', fliplr(fake_age')];
y1=ci_normal_hc(:,1);
y2=ci_normal_hc(:,2);
y = [(y1)', fliplr(y2')];%y的转置在里面，之前犯了先flip再转置的错误
fill(x, y, [0.5 0.5 0.5],'EdgeColor', 'none','facealpha',0.5);
hold on
plot(fake_age,(yfit(:,3)),'k','linewidth',2.5); %plot 50th centile
plot(fake_age_neg,(yfit_neg(:,3)),'r','linewidth',2.5); %plot 50th centile
hold on
% plot(fake_age,ci_normal_sub1(:,1),'b--','linewidth',1.5); %plot 50th centile
% hold on
% plot(fake_age,ci_normal_sub1(:,2),'b--','linewidth',1.5); %plot 50th centile

box off
set(gca, 'FontSize', 14);

%% sex
inpath3=('D:\study\sub2\brain_age\variables_for_normative_modeling\output_normative_modeling\boot_abs\');
for i=1:1000
    i
    yfitF=nan(1000,5);
    load([inpath3,...
        'GAMLSS_','t1asighc_',num2str(i),'_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
    faf_age=age;
    yfit=yfitF;
    fake_age_boot(:,i)=faf_age;
    fit_boot(:,i)=yfit(:,3);
    % sub 1
    yfitF=nan(1000,5);
    load([inpath3,...
        'GAMLSS_','t1asigsub1_',num2str(i),'_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
    faf_age=age;
    yfit=(yfitF);
    fake_age_sub1_boot(:,i)=(faf_age);
    fit_sub1_boot(:,i)=yfit(:,3);
end
% 95%ci
for i=1:1000
    mean_normal = mean(fit_boot(i,:));
    std_normal = std(fit_boot(i,:));
    confidence_level = 0.95; % 95%
    df = (length(fit_boot(i,:))-1);
    critical_value = tinv((1 + confidence_level) / 2, df);
    
    ci_normal_hc(i,:) = [mean_normal+ - critical_value * std_normal , ...
        mean_normal + critical_value * std_normal];
    
    
    mean_normal = mean(fit_sub1_boot(i,:));
    std_normal = std(fit_sub1_boot(i,:));
    
    confidence_level = 0.95; % 95%
    
    df = (length(fit_sub1_boot(i,:))-1);
    critical_value = tinv((1 + confidence_level) / 2, df);
    
    ci_normal_sub1(i,:) = [mean_normal+ - critical_value * std_normal , ...
        mean_normal + critical_value * std_normal];
end


load([inpath2,...
    'GAMLSS_','t1asighc_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
yfit=(yfitF);
fake_age=(faf_age);
plot(fake_age,(yfit(:,3)),'k','linewidth',2.5); %plot 50th centile
hold on
%gap<0
load([inpath2,...
    'GAMLSS_','t1asigsub1_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;

yfit_neg=(yfitF);
fake_age_neg=(faf_age);

x = [fake_age', fliplr(fake_age')];
y1=ci_normal_hc(:,1);
y2=ci_normal_hc(:,2);
y = [(y1)', fliplr(y2')];
fill(x, y, [0.6 0.6 0.6],'EdgeColor', 'none','facealpha',0.5);
hold on
plot(fake_age,(yfit(:,3)),'k','linewidth',2.5); %plot 50th centile
hold on

x = [fake_age', fliplr(fake_age')];
y = [ci_normal_sub1(:,1)', fliplr(ci_normal_sub1(:,2)')];
fill(x, y, [0.7, 0.7, 1], 'EdgeColor', 'none','facealpha',0.5);
hold on
plot(fake_age_neg,(yfit_neg(:,3)),'b','linewidth',2.5); %plot 50th centile
box off
set(gca, 'FontSize', 14);

%% sub 2
for i=1:1000
    yfitF=nan(1000,5);
    load([inpath3,...
        'GAMLSS_','t2asighc_',num2str(i),'_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
    faf_age=age;
    yfit=(yfitF);
    fake_age_boot(:,i)=(faf_age);
    fit_boot(:,i)=yfit(:,3);
    % sub 1
    yfitF=nan(1000,5);
    load([inpath3,...
        'GAMLSS_','t2asigsub2_',num2str(i),'_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
    faf_age=age;
    yfit=(yfitF);
    fake_age_sub1_boot(:,i)=(faf_age);
    fit_sub1_boot(:,i)=yfit(:,3);
end
% 95%ci
for i=1:1000
    mean_normal = mean(fit_boot(i,:));
    std_normal = std(fit_boot(i,:));
    
    confidence_level = 0.95; % 95%
    
    df = (length(fit_boot(i,:))-1);
    critical_value = tinv((1 + confidence_level) / 2, df);
    
    ci_normal_hc(i,:) = [mean_normal+ - critical_value * std_normal , ...
        mean_normal + critical_value * std_normal];
    
    
    mean_normal = mean(fit_sub1_boot(i,:));
    std_normal = std(fit_sub1_boot(i,:));
    
    confidence_level = 0.95; % 95%
    
    df = (length(fit_sub1_boot(i,:))-1);
    critical_value = tinv((1 + confidence_level) / 2, df);
    
    ci_normal_sub1(i,:) = [mean_normal+ - critical_value * std_normal , ...
        mean_normal + critical_value * std_normal];
end

yfitM=nan(1000,5);
load([inpath2,...
    'GAMLSS_','t2asighc_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;

yfit=yfitF;
fake_age=(faf_age);

% plot(fake_age,(yfit(:,3)),'k','linewidth',2.5); %plot 50th centile
% hold on
%gap<0
load([inpath2,...
    'GAMLSS_','t2asigsub2_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;

yfit_neg=(yfitF);
fake_age_neg=(faf_age);
%
x = [fake_age', fliplr(fake_age')];
y = [ci_normal_sub1(:,1)', fliplr(ci_normal_sub1(:,2)')];
fill(x, y, [1, 0.7, 0.7], 'EdgeColor', 'none','facealpha',0.5); % 灰色阴影
hold on
x = [fake_age', fliplr(fake_age')];
y1=ci_normal_hc(:,1);
y2=ci_normal_hc(:,2);
y = [(y1)', fliplr(y2')];%y的转置在里面，之前犯了先flip再转置的错误
fill(x, y, [0.5 0.5 0.5],'EdgeColor', 'none','facealpha',0.5);
hold on
plot(fake_age,(yfit(:,3)),'k','linewidth',2.5); %plot 50th centile
plot(fake_age_neg,(yfit_neg(:,3)),'r','linewidth',2.5); %plot 50th centile
hold on
% plot(fake_age,ci_normal_sub1(:,1),'b--','linewidth',1.5); %plot 50th centile
% hold on
% plot(fake_age,ci_normal_sub1(:,2),'b--','linewidth',1.5); %plot 50th centile
box off
set(gca, 'FontSize', 14);
%% permut
% sub 1
clear all
close all
inpath1=('D:\study\sub2\brain_age\variables_for_normative_modeling\all\');
inpath2=('D:\study\sub2\brain_age\variables_for_normative_modeling\output_normative_modeling\all\');
inpath3=('D:\study\sub2\brain_age\variables_for_normative_modeling\output_normative_modeling\permut_abs_new\');
for i=1:1000
    i
    yfitF=nan(1000,5);
    yfitM=nan(1000,5);
    load([inpath3,...
        'GAMLSS_','t1psighc_',num2str(i),'_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
    faf_age=age;
    load([inpath3,...
        'GAMLSS_','t1psighc_',num2str(i),'_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
    fam_age=age;
    yfit=(yfitF+yfitM)/2;
    fake_age_permut(:,i)=(faf_age+fam_age)/2;
    fit_permut(:,i)=yfit(:,3);
    % sub 1
    yfitF=nan(1000,5);
    yfitM=nan(1000,5);
    load([inpath3,...
        'GAMLSS_','t1psigsub1_',num2str(i),'_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
    faf_age=age;
    load([inpath3,...
        'GAMLSS_','t1psigsub1_',num2str(i),'_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
    fam_age=age;
    yfit=(yfitF+yfitM)/2;
    fit_sub1_permut(:,i)=yfit(:,3);
end

yfitM=nan(1000,5);
load([inpath2,...
    'GAMLSS_','t1psighc_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','t1psighc_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit=(yfitF+yfitM)/2;
real_hc=yfit(:,3);

%gap<0
load([inpath2,...
    'GAMLSS_','t1psigsub1_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','t1psigsub1_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles;
fam_age=age;
yfit_neg=(yfitF+yfitM)/2;
real_sub1_pos=yfit_neg(:,3);
real_diff_pos_sub1=(real_sub1_pos-real_hc);
permut_diff_pos=[];
for i=1:1000
    permut_diff_pos(i,:)=(fit_sub1_permut(:,i)-fit_permut(:,i)); 
end
for i=11:64
    pos=find(faf_age>=i&faf_age<i+1);
    real_sub1_pos1=mean(real_diff_pos_sub1(pos))
    permut_diff_pos1=mean(permut_diff_pos(:,pos),2); 
    sig_sub1_pos(i-10)=length(find( real_sub1_pos1< permut_diff_pos1))/1000;
end
find(sig_sub1_pos<0.05)
q=mafdr(sig_sub1_pos,'BHFDR', true);

%% sub2
clear all
close all
inpath1=('D:\study\sub2\brain_age\variables_for_normative_modeling\all\');
inpath2=('D:\study\sub2\brain_age\variables_for_normative_modeling\output_normative_modeling\all\');
inpath3=('D:\study\sub2\brain_age\variables_for_normative_modeling\output_normative_modeling\permut_abs_new\');
for i=1:1000
    i
    yfitF=nan(1000,5);
    yfitM=nan(1000,5);
    load([inpath3,...
        'GAMLSS_','t2psighc_',num2str(i),'_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
    faf_age=age;
    load([inpath3,...
        'GAMLSS_','t2psighc_',num2str(i),'_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
    fam_age=age;
    yfit=(yfitF+yfitM)/2;
    fake_age_permut(:,i)=(faf_age+fam_age)/2;
    fit_permut(:,i)=yfit(:,3);
    % sub 1
    yfitF=nan(1000,5);
    yfitM=nan(1000,5);
    load([inpath3,...
        'GAMLSS_','t2psigsub2_',num2str(i),'_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
    faf_age=age;
    load([inpath3,...
        'GAMLSS_','t2psigsub2_',num2str(i),'_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
    yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
    fam_age=age;
    yfit=(yfitF+yfitM)/2;
    fit_sub1_permut(:,i)=yfit(:,3);
end

yfitM=nan(1000,5);
load([inpath2,...
    'GAMLSS_','t2psighc_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','t2psighc_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit=(yfitF+yfitM)/2;
real_hc=yfit(:,3);

%gap<0
load([inpath2,...
    'GAMLSS_','t2psigsub2_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','t2psigsub2_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit_neg=(yfitF+yfitM)/2;
real_sub1_pos=yfit_neg(:,3);
real_diff_pos_sub1=(real_sub1_pos-real_hc);
permut_diff_pos=[];
for i=1:1000
    permut_diff_pos(i,:)=(fit_sub1_permut(:,i)-fit_permut(:,i));
end
% for i=1:1000
% sig_sub_pos(i)=length(find( real_diff_pos_sub1(i)< permut_diff_pos(i,:)))/1000;
% end
for i=11:64
    pos=find(faf_age>=i&faf_age<i+1);
    real_sub1_pos1=mean(real_diff_pos_sub1(pos));
    permut_diff_pos1=mean(permut_diff_pos(:,pos),2);
    
    sig_sub1_pos(i-10)=length(find( real_sub1_pos1>=permut_diff_pos1))/1000;
end
find(sig_sub1_pos<0.05)
q=mafdr(sig_sub1_pos,'BHFDR', true);
