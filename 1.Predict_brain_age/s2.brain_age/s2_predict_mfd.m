addpath(genpath('/home/cxpang/matlab/code/8.brain_age'));
method='NGSR'
load dida
voxel_num=45892;
network_num=18;
data_dir='/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/';
addpath(genpath('/home/cxpang/matlab/SVR'));
load(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/brain_age/loading_mdd_remove65.mat'),'loading')
mdd_loading=loading;
load(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/brain_age/loading_hc_remove65.mat'))
%%
load(strcat('/home/cxpang/matlab/code//8.brain_age/using_mat/covariate_hc_remove65.mat'))
delete_pos=find(mfd>=0.25);
sex(delete_pos)=[];
mfd(delete_pos)=[];
age(delete_pos)=[];
loading(delete_pos,:)=[];
sbj_num=length(sex);

Subjects_Data=loading;
Subjects_Scores =age;
covariate=[sex];
ResultantFolder = [data_dir,'/SVR_all_no_regress_fold20_remove65sex_deletemfd'];
FoldQuantity =20;
Pre_Method = 'Normalize';
C_Range =1%power(10, -5:5);
Weight_Flag = 0;
Permutation_Flag = 0;
randNet = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18];
[Prediction] = pca_svr(Subjects_Data, Subjects_Scores,covariate,FoldQuantity, Pre_Method, C_Range, network_num,voxel_num,voxel_num,network_num,randNet,Weight_Flag, Permutation_Flag, ResultantFolder);
%% all hc to mdd model
ResultantFolder = [data_dir,'/SVR_all_no_regress_fold20_remove65sex_deletemfd'];
load([ResultantFolder filesep 'Prediction.mat'], 'Prediction');
Pre_Method = 'Normalize';
[slm,model,MeanValue,StandardDeviation] = SVR_all(covariate,Subjects_Data, Subjects_Scores,Pre_Method,0, ResultantFolder);
clear Subjects_Data
clear Subjects_Scores
load(strcat('/home/cxpang/matlab/code/8.brain_age/using_mat/covariate_mdd_remove65.mat'))
delete_pos=find(mfd>=0.25);
sex(delete_pos)=[];
mfd(delete_pos)=[];
age(delete_pos)=[];
mdd_loading(delete_pos,:)=[];
sbj_num=length(sex);

test_data=mdd_loading;
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
%% mdd gap
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
% compute the age gap
age_real=test_score;
age_predic=Predicted_Scores;
mdd_all_gap=age_predic-age_real;
yfit=glmval(beta,age_real,'identity');
adjusted_mdd_all_gap_mfd=mdd_all_gap-yfit; % corrected age gap
load('/home/cxpang/matlab/code/8.brain_age/adjusted_mdd_all_gap_remove65sex.mat','adjusted_mdd_all_gap')
load(strcat('/home/cxpang/matlab/code/8.brain_age/using_mat/covariate_mdd_remove65.mat'))
pos=find(mfd<0.25);
tmp=adjusted_mdd_all_gap(pos);
number_center=0;
for j=1:length(pos)
    if(sign(tmp(j))==sign(adjusted_mdd_all_gap_mfd(j)))
        number_center=number_center+1;
    end
end
percent=number_center/length(pos);
save(strcat('/home/cxpang/matlab/code/8.brain_age/adjusted_mdd_all_gap_remove65_mfd.mat'),'adjusted_mdd_all_gap_mfd')