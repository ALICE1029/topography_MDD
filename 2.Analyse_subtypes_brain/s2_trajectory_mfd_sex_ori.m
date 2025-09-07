
%%
addpath(genpath('/home/cxpang/matlab/code/8.brain_age'));
addpath(genpath('/home/cxpang/matlab/Collaborative_Brain_Decomposition-master/lib/NIfTI_20140122'))
load dida
voxel_num=45892;
network_num=18;
method='NGSR';
load(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/after_combat/variable_combat_remove65.mat'),'x_tmp')
lo=x_tmp;
for s=1:2
    load(strcat('/home/cxpang/matlab/code/8.brain_age/adjusted_mdd_all_gap_remove65sex',num2str(s),'.mat'),'adjusted_mdd_all_gap_sex')
    
    load ('/HeLabData2/cxpang/DIDA/using_mat/site_remove65.mat')%age,group,mfd,sex
    load ('/HeLabData2/cxpang/DIDA/using_mat/glmvariable_remove65.mat')%age,group,mfd,sex
  
    pos_hc=intersect(find(sex==s),find(group==1));
        pos_mdd=intersect(find(sex==s),find(group==2));

  sub1_index=pos_mdd(find(adjusted_mdd_all_gap_sex<0));
    sub2_index=pos_mdd(find(adjusted_mdd_all_gap_sex>0));
    
    sex=sex-1;
    hc_sex=sex(pos_hc);
    hc_mfd=mfd(pos_hc);
    hc_age=age(pos_hc);
    
    mdd1_sex=(sex(sub1_index));
    mdd1_mfd=mfd(sub1_index);
    mdd1_age=age(sub1_index);
    mdd2_sex=(sex(sub2_index));
    mdd2_mfd=mfd(sub2_index);
    mdd2_age=age(sub2_index);
    
    
    
    
    
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
    %posmax_sub1=find(t1(:)==max(t1(:)));
    neg_sub1=find(t1(:)<0);
    lo_pos=lo(pos_sub1,:);
    lo_neg=lo(neg_sub1,:);
    
    
    
    weight_p_pos=pca(lo_pos');
    
    weight_pca_pos = lo_pos'* weight_p_pos(:,:);
    weight_pca_pos=weight_pca_pos(:,1);
    
    weight_p_neg=pca(lo_neg');
    
    weight_pca_neg = lo_neg'* weight_p_neg(:,:);
    weight_pca_neg=weight_pca_neg(:,1);
    
    
   
    if(min(weight_pca_pos)<0)
         weight_new=weight_pca_pos+abs(min(weight_pca_pos))+0.1;
    else
          weight_new=weight_pca_pos;    
    end
   
    
    hc_weight=weight_new(pos_hc);
    
    sub1_weight=weight_new(sub1_index);
    
    
    vars_loso('all','osext1psighc', hc_weight,hc_sex, hc_age, hc_mfd,s)
    vars_loso('all','osext1psigsub1', sub1_weight,mdd1_sex, mdd1_age, mdd1_mfd,s)
    
    if(min(weight_pca_neg)<0)
        weight_new=weight_pca_neg+abs(min(weight_pca_neg))+0.1;
    else
        weight_new=weight_pca_neg;
    end
        hc_weight=weight_new(pos_hc);
    sub1_weight=weight_new(sub1_index);
    vars_loso('all','osext1npsighc', hc_weight,hc_sex, hc_age, hc_mfd,s)
    vars_loso('all','osext1npsigsub1', sub1_weight,mdd1_sex, mdd1_age, mdd1_mfd,s)
    
    %%  sub 2
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
    all_sub2=find(t2(:)~=0);
    pos_sub2=find(t2(:)>0);
    neg_sub2=find(t2(:)<0);
  
    lo_pos=lo(pos_sub2,:);
    lo_neg=lo(neg_sub2,:);
    
    weight_p_pos2=pca(lo_pos');
    
    weight_pca_pos2 = lo_pos'* weight_p_pos2(:,:);
    weight_pca_pos2=weight_pca_pos2(:,1);
    
    weight_p_neg2=pca(lo_neg');
    
    weight_pca_neg2 = lo_neg'* weight_p_neg2(:,:);
    weight_pca_neg2=weight_pca_neg2(:,1);
    
    if(min(weight_pca_pos2)<0)
         weight_new=weight_pca_pos2+abs(min(weight_pca_pos2))+0.1;
    else
          weight_new=weight_pca_pos2;    
    end
   
    
    hc_weight=weight_new(pos_hc);
    sub2_weight=weight_new(sub2_index);
    
    
    vars_loso('all','osext2psighc', hc_weight,hc_sex, hc_age, hc_mfd,s)
    vars_loso('all','osext2psigsub2', sub2_weight,mdd2_sex, mdd2_age, mdd2_mfd,s)
    
    if(min(weight_pca_neg2)<0)
        weight_new=weight_pca_neg2+abs(min(weight_pca_neg2))+0.1;
    else
        weight_new=weight_pca_neg2;
    end
        hc_weight=weight_new(pos_hc);
    sub2_weight=weight_new(sub2_index);
    vars_loso('all','osext2npsighc', hc_weight,hc_sex, hc_age, hc_mfd,s)
    vars_loso('all','osext2npsigsub2', sub2_weight,mdd2_sex, mdd2_age, mdd2_mfd,s)
    
end
%% sex
for s=1:2
    
    
    
    
    
    
    yfit=nan(1000,5);
    load(strcat(inpath2,...
        'GAMLSS_','sext1asighc_',num2str(s),'_predicted_WHOLE_SAMPLE.mat'),'predictions_quantiles','age');
    yfit=predictions_quantiles; %estimated females at site i (assumes females = 0)
    fake_age=age;
    %plot results (average across males and females and across sites)
    hf=figure; hf.Color='w'; hf.Position=[50,50,800,400];
    
    hold on
    plot(fake_age,(yfit(:,3)),'k','linewidth',2.5); %plot 50th centile
    hold on
    yfit=nan(1000,5);
    load(strcat(inpath2,...
        'GAMLSS_','sext1asigsub1_',num2str(s),'_predicted_WHOLE_SAMPLE.mat'),'predictions_quantiles','age');
    yfit=predictions_quantiles; %estimated females at site i (assumes females = 0)
    fake_age=age;
    %plot results (average across males and females and across sites)
    hold on
    plot(fake_age,(yfit(:,3)),'b','linewidth',2.5); %plot 50th centile
    set(gca, 'FontSize', 14);
    outputFileName = strcat(outpath, 'GAMLSS_plot_sub1_abs_sex',num2str(s),'.png'); % 定义输出文件名
    saveas(hf, outputFileName, 'png'); % 保存为 PNG 格式
    
    
    
    
    
    yfit=nan(1000,5);
    load(strcat(inpath2,...
        'GAMLSS_','sext2asighc_',num2str(s),'_predicted_WHOLE_SAMPLE.mat'),'predictions_quantiles','age');
    yfit=predictions_quantiles; %estimated females at site i (assumes females = 0)
    fake_age=age;
    %plot results (average across males and females and across sites)
    hf=figure; hf.Color='w'; hf.Position=[50,50,800,400];
    hold on
    plot(fake_age,(yfit(:,3)),'k','linewidth',2.5); %plot 50th centile
    hold on
    yfit=nan(1000,5);
    load(strcat(inpath2,...
        'GAMLSS_','sext2asigsub2_',num2str(s),'_predicted_WHOLE_SAMPLE.mat'),'predictions_quantiles','age');
    yfit=predictions_quantiles; %estimated females at site i (assumes females = 0)
    fake_age=age;
    %plot results (average across males and females and across sites)
    hold on
    plot(fake_age,(yfit(:,3)),'r','linewidth',2.5); %plot 50th centile
    set(gca, 'FontSize', 14);
    outputFileName = strcat(outpath, 'GAMLSS_plot_sub2_abs_sex',num2str(s),'.png'); % 定义输出文件名
    saveas(hf, outputFileName, 'png'); % 保存为 PNG 格式
    
end

%% mfd
load(strcat('/home/cxpang/matlab/code/8.brain_age/adjusted_mdd_all_gap_remove65_mfd.mat'),'adjusted_mdd_all_gap_mfd')

load ('/HeLabData2/cxpang/DIDA/using_mat/site_remove65.mat')%age,group,mfd,sex
load ('/HeLabData2/cxpang/DIDA/using_mat/glmvariable_remove65.mat')%age,group,mfd,sex
pos_center=find(mfd<0.25);
pos_hc=intersect(pos_center,find(group==1));
pos_mdd=intersect(pos_center,find(group==2));

sub1_index=pos_mdd(find(adjusted_mdd_all_gap_mfd<0));
sub2_index=pos_mdd(find(adjusted_mdd_all_gap_mfd>0));

sex=sex-1;
hc_sex=sex(pos_hc);
hc_mfd=mfd(pos_hc);
hc_age=age(pos_hc);

mdd1_sex=(sex(sub1_index));
mdd1_mfd=mfd(sub1_index);
mdd1_age=age(sub1_index);
mdd2_sex=(sex(sub2_index));
mdd2_mfd=mfd(sub2_index);
mdd2_age=age(sub2_index);





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
all_sub1=find(t1(:)~=0);
pos_sub1=find(t1(:)>0);
%posmax_sub1=find(t1(:)==max(t1(:)));
neg_sub1=find(t1(:)<0);
lo_pos=lo(pos_sub1,:);
lo_neg=lo(neg_sub1,:);



weight_p_pos=pca(lo_pos');

weight_pca_pos = lo_pos'* weight_p_pos(:,:);
weight_pca_pos=weight_pca_pos(:,1);

weight_p_neg=pca(lo_neg');

weight_pca_neg = lo_neg'* weight_p_neg(:,:);
weight_pca_neg=weight_pca_neg(:,1);



if(min(weight_pca_pos)<0)
    weight_new=weight_pca_pos+abs(min(weight_pca_pos))+0.1;
else
    weight_new=weight_pca_pos;
end


hc_weight=weight_new(pos_hc);

sub1_weight=weight_new(sub1_index);


vars_loso('all','omfdt1psighc', hc_weight,hc_sex, hc_age, hc_mfd,1)
vars_loso('all','omfdt1psigsub1', sub1_weight,mdd1_sex, mdd1_age, mdd1_mfd,1)

if(min(weight_pca_neg)<0)
    weight_new=weight_pca_neg+abs(min(weight_pca_neg))+0.1;
else
    weight_new=weight_pca_neg;
end
hc_weight=weight_new(pos_hc);
sub1_weight=weight_new(sub2_index);
vars_loso('all','omfdt1npsighc', hc_weight,hc_sex, hc_age, hc_mfd,1)
vars_loso('all','omfdt1npsigsub1', sub1_weight,mdd1_sex, mdd1_age, mdd1_mfd,1)

%%  sub 2
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
all_sub2=find(t2(:)~=0);
pos_sub2=find(t2(:)>0);
neg_sub2=find(t2(:)<0);

lo_pos=lo(pos_sub2,:);
lo_neg=lo(neg_sub2,:);

weight_p_pos2=pca(lo_pos');

weight_pca_pos2 = lo_pos'* weight_p_pos2(:,:);
weight_pca_pos2=weight_pca_pos2(:,1);

weight_p_neg2=pca(lo_neg');

weight_pca_neg2 = lo_neg'* weight_p_neg2(:,:);
weight_pca_neg2=weight_pca_neg2(:,1);

if(min(weight_pca_pos2)<0)
    weight_new=weight_pca_pos2+abs(min(weight_pca_pos2))+0.1;
else
    weight_new=weight_pca_pos2;
end


hc_weight=weight_new(pos_hc);
sub2_weight=weight_new(sub2_index);


vars_loso('all','omfdt2psighc', hc_weight,hc_sex, hc_age, hc_mfd,1)
vars_loso('all','omfdt2psigsub2', sub2_weight,mdd2_sex, mdd2_age, mdd2_mfd,1)

if(min(weight_pca_neg2)<0)
    weight_new=weight_pca_neg2+abs(min(weight_pca_neg2))+0.1;
else
    weight_new=weight_pca_neg2;
end
hc_weight=weight_new(pos_hc);
sub2_weight=weight_new(sub2_index);
vars_loso('all','omfdt2npsighc', hc_weight,hc_sex, hc_age, hc_mfd,1)
vars_loso('all','omfdt2npsigsub2', sub2_weight,mdd2_sex, mdd2_age, mdd2_mfd,1)
%% mfd
yfitF=nan(1000,5);
yfitM=nan(1000,5);
load([inpath2,...
    'GAMLSS_','mfdt1asighc_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','mfdt1asighc_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit=(yfitF+yfitM)/2;
fake_age=(faf_age+fam_age)/2;

%plot results (average across males and females and across sites)
hf=figure; hf.Color='w'; hf.Position=[50,50,800,400];
hold on
plot(fake_age,(yfit(:,3)),'k','linewidth',2.5); %plot 50th centile
fit=yfit(:,3);

hold on
% gap>0
load([inpath2,...
    'GAMLSS_','mfdt1asigsub1_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','mfdt1asigsub1_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit_pos=(yfitF+yfitM)/2;
fake_age_pos=(faf_age+fam_age)/2;
plot(fake_age_pos,(yfit_pos(:,3)),'b','linewidth',2.5); %plot 50th centile
set(gca, 'FontSize', 14);
outputFileName = strcat(outpath, 'GAMLSS_plot_sub1_abs_mfd.png'); % 定义输出文件名
saveas(hf, outputFileName, 'png'); % 保存为 PNG 格式


yfitF=nan(1000,5);
yfitM=nan(1000,5);
load([inpath2,...
    'GAMLSS_','mfdt2asighc_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','mfdt2asighc_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit=(yfitF+yfitM)/2;
fake_age=(faf_age+fam_age)/2;

%plot results (average across males and females and across sites)
hf=figure; hf.Color='w'; hf.Position=[50,50,800,400];
hold on
plot(fake_age,(yfit(:,3)),'k','linewidth',2.5); %plot 50th centile
fit=yfit(:,3);

hold on
% gap>0
load([inpath2,...
    'GAMLSS_','mfdt2asigsub2_1_predicted_sex0_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitF=predictions_quantiles; %estimated females at site i (assumes females = 0)
faf_age=age;
load([inpath2,...
    'GAMLSS_','mfdt2asigsub2_1_predicted_sex1_WHOLE_SAMPLE.mat'],'predictions_quantiles','age');
yfitM=predictions_quantiles; %estimated males at site i (assumes males = 1)
fam_age=age;
yfit_pos=(yfitF+yfitM)/2;
fake_age_pos=(faf_age+fam_age)/2;
fit_pos=yfit_pos(:,3);
plot(fake_age_pos,(yfit_pos(:,3)),'r','linewidth',2.5); %plot 50th centile
set(gca, 'FontSize', 14);
outputFileName = strcat(outpath, 'GAMLSS_plot_sub2_abs_mfd.png'); % 定义输出文件名
saveas(hf, outputFileName, 'png'); % 保存为 PNG 格式