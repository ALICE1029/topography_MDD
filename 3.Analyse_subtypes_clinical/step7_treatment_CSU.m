load CSU
%%
addpath(genpath('/home/cxpang/matlab/code/8.brain_age'));
method='NGSR'
voxel_num=45892;
network_num=18;
data_dir='/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/';
addpath(genpath('/home/cxpang/matlab/SVR'));
%%
load('/HeLabData2/cxpang/treatment_NGSR/CSU/result/regress_variable.mat','x')
load('/HeLabData2/cxpang/treatment_NGSR/CSU/predict/deltahamd.mat')
load('/HeLabData2/cxpang/treatment_NGSR/CSU/predict/covariate.mat')
load('/HeLabData2/cxpang/treatment_NGSR/CSU/predict/mfd.mat','mfd')
load('/HeLabData2/cxpang/treatment_NGSR/CSU/predict/mfd_co.mat','mfd_base','mfd_after')
%% gap of treatment group
load('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/hc_all_model_remove65sex_slm.mat','slm')
load('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/hc_all_model_remove65sex.mat')
load('/HeLabData2/cxpang/DIDA/using_mat/pos_zero.mat','pos')
x(pos,:)=[];
load('/home/cxpang/matlab/treatment/CSU_baseline.mat','pos')
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
% regress out age dependence using beta computed from the training sets
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
%%
for i=1:43
    if(adjusted_mdd_all_gap1(i)<0)
        group(i)=1;
    else
        group(i)=2;
    end
end
group = group';
des=[mfd_base];
[adjusted_mdd_all_gap1, b, stats] = regress_out(adjusted_mdd_all_gap1, des);
des=[mfd_after];
[adjusted_mdd_all_gap, b, stats] = regress_out(adjusted_mdd_all_gap, des);
dataall =[(adjusted_mdd_all_gap1),(adjusted_mdd_all_gap)];


tbl = array2table(dataall, 'VariableNames', {'Day1', 'Day2'});
tbl.Group = group;
tbl.Group = categorical(tbl.Group);
withinDesign = table([1 2]', 'VariableNames', {'Time'});
rm = fitrm(tbl, 'Day1, Day2 ~ Group', 'WithinDesign', withinDesign);
ranovatbl2 = ranova(rm);
disp(ranovatbl2);
[t,p]=ttest(adjusted_mdd_all_gap1(find(group==2)),adjusted_mdd_all_gap(find(group==2)));
[t,p]=ttest(adjusted_mdd_all_gap1(find(group==1)),adjusted_mdd_all_gap(find(group==1)));
mcomp = multcompare(rm, 'Time', 'By', 'Group');
disp(mcomp);
%% 
load('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/hamd_base.mat')
delta1=dataall(find(group==2),2)-dataall(find(group==2),1);
deltahamd2=deltahamd(find(group==2));
[r,p]=corr(delta1,deltahamd2./hamd(find(group==2)))

delta1=dataall(find(group==1),2)-dataall(find(group==1),1);
deltahamd1=deltahamd(find(group==1));
[r,p]=corr(delta1,deltahamd1./hamd(find(group==1)))
%% topography changes
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
% generate change_sub1/2 folders, 18folder with sub_num subfolders, for
% paired ttest with GRF corrected using dpabi