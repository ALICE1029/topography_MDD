load CSU

%%
method='NGSR'
voxel_num=45892;
network_num=18;
data_dir='/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/';
addpath(genpath('/home/cxpang/matlab/SVR'));
%%
load('/HeLabData2_2/cxpang/treatment_NGSR/CSU/result/regress_variable.mat','x')
load('/HeLabData2_2/cxpang/treatment_NGSR/CSU/predict/deltahamd.mat')
load('/HeLabData2_2/cxpang/treatment_NGSR/CSU/predict/covariate.mat')
load('/HeLabData2_2/cxpang/treatment_NGSR/CSU/predict/mfd.mat','mfd')
load('/HeLabData2_2/cxpang/treatment_NGSR/CSU/predict/mfd_co.mat','mfd_base','mfd_after')
%% gap of treatment group
load('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/hc_all_model_remove65sex_slm.mat','slm')
load('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/hc_all_model_remove65sex.mat')
load('/HeLabData2/cxpang/DIDA/using_mat/pos_zero.mat','pos')
x(pos,:)=[];
load('/home/cxpang/matlab/treatment/CSU_baseline.mat')
lo1=x(:,pos);
load('/home/cxpang/matlab/treatment/CSU_after.mat','pos')
lo2=x(:,pos);



test_data=lo1';
mfd_after=mfd(pos);
Covariates_test=[sex];
test_score=age;

%regress
Covariates_quantity=1;
[test_quantity, Feature_Quantity] = size(test_data);
test_data = test_data - repmat(slm.coef(1, :), test_quantity, 1);
for k = 1:Covariates_quantity
    test_data = test_data - ...
        repmat(Covariates_test(:, k), 1, Feature_Quantity) .* repmat(slm.coef(k + 1, :), test_quantity, 1);
end
%normalize
MeanValue_New = repmat(MeanValue, length(test_score), 1);
StandardDeviation_New = repmat(StandardDeviation, length(test_score), 1);
test_data_norm = (test_data - MeanValue_New) ./ StandardDeviation_New;
test_data_final = double(test_data_norm);
[Predicted_Scores, ~, ~] = svmpredict(test_score, test_data_final, model);
%yhat=predict(model,test_data);
r_predic=corr(Predicted_Scores,test_score);
save('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/csutreatment_res_remove65sex_baseline.mat','Predicted_Scores','test_score')
test_data=lo2';
Covariates_test=[sex];
test_score=age;
%regress
Covariates_quantity=1;
[test_quantity, Feature_Quantity] = size(test_data);
test_data = test_data - repmat(slm.coef(1, :), test_quantity, 1);
for k = 1:Covariates_quantity
    test_data = test_data - ...
        repmat(Covariates_test(:, k), 1, Feature_Quantity) .* repmat(slm.coef(k + 1, :), test_quantity, 1);
end
%normalize
MeanValue_New = repmat(MeanValue, length(test_score), 1);
StandardDeviation_New = repmat(StandardDeviation, length(test_score), 1);
test_data_norm = (test_data - MeanValue_New) ./ StandardDeviation_New;
test_data_final = double(test_data_norm);
[Predicted_Scores, ~, ~] = svmpredict(test_score, test_data_final, model);
%yhat=predict(model,test_data);
r_predic=corr(Predicted_Scores,test_score);
save('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/csutreatment_res_remove65sex_after.mat','Predicted_Scores','test_score')
%%
load('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/csutreatment_res_remove65sex_after.mat','Predicted_Scores','test_score')

load(strcat(data_dir,'/SVR_all_no_regress_fold20_remove65sex/Prediction.mat'))
x=[];
y=[];
id=[];
for i=1:20
x = [x;Prediction.Score{i}];
y= [y;Prediction.realScore{i}];
id=[id;Prediction.Origin_ID{i}];
end
all_gap= x-y;
beta=glmfit(y,all_gap);
age_real=test_score;
age_predic=Predicted_Scores;
mdd_all_gap=age_predic-age_real;
% regress out age dependence using beta computed from the training set
my_corr=corr(age_predic,age_real)
mae=mean(abs(mdd_all_gap))
% adjusted_mdd_age=(age_predic-d)/k;
% adjusted_mdd_all_gap=adjusted_mdd_age-age_real;
% regress 2 
yfit=glmval(beta,age_real,'identity');
adjusted_mdd_all_gap=mdd_all_gap-yfit; % corrected age gap
save('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/csuadgap_after.mat','adjusted_mdd_all_gap')
%%
load('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/csutreatment_res_remove65sex_baseline.mat','Predicted_Scores','test_score')
load(strcat(data_dir,'/SVR_all_no_regress_fold20_remove65sex/Prediction.mat'))
x=[];
y=[];
id=[];
for i=1:20
x = [x;Prediction.Score{i}];
y= [y;Prediction.realScore{i}];
id=[id;Prediction.Origin_ID{i}];
end
all_gap= x-y;
beta=glmfit(y,all_gap);
age_real=test_score;
age_predic=Predicted_Scores;
mdd_all_gap=age_predic-age_real;
% regress out age dependence using beta computed from the training set
my_corr=corr(age_predic,age_real)
mae=mean(abs(mdd_all_gap))
% adjusted_mdd_age=(age_predic-d)/k;
% adjusted_mdd_all_gap=adjusted_mdd_age-age_real;
% regress 2 
yfit=glmval(beta,age_real,'identity');
adjusted_mdd_all_gap1=mdd_all_gap-yfit; % corrected age gap
save('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/csuadgap_base.mat','adjusted_mdd_all_gap1')


%% topography changes
% fisrt: find voxel that significant correlated with gap
load('/HeLabData2/cxpang/treatment_NGSR/CSU/result/regress_variable.mat','x')
load('/home/cxpang/matlab/treatment/CSU_baseline.mat')
lo1=x(:,pos)';
load('/home/cxpang/matlab/treatment/CSU_after.mat','pos')
lo2=x(:,pos)';
load('/HeLabData2/cxpang/treatment_NGSR/CSU/predict/group.mat','group')
for i=1:826056
    lo1(:,i)=regress_out(lo1(:,i), mfd_base);
       lo2(:,i)=regress_out(lo2(:,i), mfd_after);
end
for i=1:826056
     [h,p1(i),ci,stats]= ttest(lo1(find(group==1),i),lo2(find(group==1),i)) ;
    r1(i)=stats.tstat;
    [h,p2(i),ci,stats]= ttest(lo1(find(group==2),i),lo2(find(group==2),i)) ;
    r2(i)=stats.tstat;
end
save('/home/cxpang/matlab/code/new_brain_age/topo_change_csu.mat','r1','p1','r2','p2')
p1(isnan(p1)) = 1;
r1(isnan(r1)) = 0;
p2(isnan(p2)) =1;
r2(isnan(r2)) = 0;
p1_new=reshape(p1,[45892,18]);
p2_new=reshape(p2,[45892,18]);
r1_new=reshape(r1,[45892,18]);
r2_new=reshape(r2,[45892,18]);
for i=1：18
    tmp=p1_new(:,i);
    pos=find(tmp>0.05);
    p1_new(pos,i)=1;
    outDir = strcat ('/HeLabData2/cxpang/DIDA/',method,'/result_xy/change_sub1');
    saveFig = 0;
    maskName ='/HeLabData2/cxpang//fmri/mask/GMMask_3mm.nii';
    maskNii = load_untouch_nii(maskName);
    if ~exist(outDir,'dir')
        mkdir(outDir);
    end
    kNii = maskNii;
    kNii.img(maskNii.img~=0) = p1_new(:,i);
    outName = [outDir,filesep,strcat('p',num2str(i),'.nii.gz')];
    save_untouch_nii(kNii,outName);
        tmp=p2_new(:,i);
    pos=find(tmp>0.05);
    p2_new(pos,i)=1;
      outDir = strcat ('/HeLabData2/cxpang/DIDA/',method,'/result_xy/change_sub2');
    saveFig = 0;
    maskName ='/HeLabData2/cxpang//fmri/mask/GMMask_3mm.nii';
    maskNii = load_untouch_nii(maskName);
    if ~exist(outDir,'dir')
        mkdir(outDir);
    end
    kNii = maskNii;
    kNii.img(maskNii.img~=0) = p2_new(:,i);
    outName = [outDir,filesep,strcat('p',num2str(i),'.nii.gz')];
    save_untouch_nii(kNii,outName);
end
