%%
addpath(genpath('/home/cxpang/matlab/surfstat'));
addpath(genpath('/home/cxpang/matlab/code/8.brain_age'));
method='NGSR'
load dida
voxel_num=45892;
network_num=18;
data_dir='/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/';
addpath(genpath('/home/cxpang/matlab/SVR'));
%%
load(strcat(data_dir,'loading_hc_remove65.mat'))
Subjects_Data=loading;
load(strcat('/home/cxpang/matlab/code//8.brain_age/using_mat/covariate_hc_remove65.mat'))
Subjects_Scores =age;
covariate=[sex];
ResultantFolder = [data_dir,'/SVR_all_no_regress_fold20_remove65sex'];
FoldQuantity =20;
Pre_Method = 'Normalize';
C_Range =1%power(10, -5:5);
Weight_Flag = 0;    `
Permutation_Flag = 0;
randNet = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18];
[Prediction] = pca_svr(Subjects_Data, Subjects_Scores,covariate,FoldQuantity, Pre_Method, C_Range, network_num,voxel_num,voxel_num,network_num,randNet,Weight_Flag, Permutation_Flag, ResultantFolder);
x=[];
y=[];
for i=1:20
    x = [x;Prediction.Score{i}];
    y= [y;Prediction.realScore{i}];
end
my_corr=corr(x,y)
mae=mean(abs(x-y));
all_gap= x-y;
mean(all_gap)
std(all_gap)
%% hc gap analyse
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
yfit=glmval(beta,y,'identity');
adjusted_all_gap=all_gap-yfit; % corrected age gap
load ('/HeLabData2/cxpang/DIDA/using_mat/site_remove65.mat')
load ('/HeLabData2/cxpang/DIDA/using_mat/glmvariable_remove65.mat')
hc=find(group==1);
site=site(hc);
for i=1:length(id)
    site_hc(i)=site(id(i));
end
for i=1:length(id)
    age_hc(i)=age(id(i));
end
adjusted_all_gap_age=cell(8,1);
ad_gap=[];
corr_age=[];
mae_age=[];
thre=10;
for i=1:6
    pos=find(y>=thre&y<thre+10);
    thre=thre+10;
    x_tmp=x(pos);
    y_tmp=y(pos);
    adjusted_all_gap_age{i}=adjusted_all_gap(pos);
    ad_gap(i)=mean(abs(adjusted_all_gap(pos)));
    mae_age(i)=mean(abs(x_tmp-y_tmp));
end
site_hc=site_hc';
adjusted_all_gap_site=cell(10,1);
ad_gap=[];
mae_site=[];
for i=1:10
    pos=find(site_hc==i);
    x_tmp=x(pos);
    y_tmp=y(pos);
    save(strcat('/home/cxpang/matlab/code/8.brain_age/plot/',num2str(i),'xysex'),'x_tmp','y_tmp')
    adjusted_all_gap_site{i}=adjusted_all_gap(pos);
    ad_gap(i)=mean(abs(adjusted_all_gap(pos)));
    mae_site(i)=mean(abs(x_tmp-y_tmp));
    mean_age(i)=mean(y(pos));
    std_age(i)=std(y(pos));
end
%% regress age 
load(strcat(data_dir,'/SVR_all_no_regress_fold20_remove65sex/Prediction.mat'))
x_all=[];
y_all=[];
id=[];
for i=1:20
    x_all = [x_all;Prediction.Score{i}];
    y_all= [y_all;Prediction.realScore{i}];
    id=[id;Prediction.Origin_ID{i}];
end
hc_gap=x_all-y_all;
data = [id, x_all];
sorted_data = sortrows(data, 1);
x_all=sorted_data(:,2);
data = [id, y_all];
sorted_data = sortrows(data, 1);
y_all=sorted_data(:,2);

ad_gap_new=[];
for i=1:20
    x = Prediction.Score{i};
    y= Prediction.realScore{i};
    id=Prediction.Origin_ID{i};
    for j=1:1065
        tmp(j)=j;
    end
    
    id_diff=setdiff(tmp,id);
    
    x_train=x_all(id_diff);
    y_train=y_all(id_diff);
    gap=x_train-y_train;
    beta=glmfit(y_train,gap);
    ad_gap=x-y;
    yfit=glmval(beta,y,'identity');
    ad_gap_new=[ad_gap_new;ad_gap-yfit]; % corrected age gap
end
id=[];
for i=1:20
    id=[id;Prediction.Origin_ID{i}];
end
data = [id, ad_gap_new];
sorted_data = sortrows(data, 1);
ad_gap_new_sort=sorted_data(:,2);

mean(abs(ad_gap_new_sort))%adj mae
load ('/HeLabData2/cxpang/DIDA/using_mat/site_remove65.mat')
load ('/HeLabData2/cxpang/DIDA/using_mat/glmvariable_remove65.mat')
hc=find(group==1);
site=site(hc);
age=age(hc);
thre=10;
for i=1:6
    pos=find(age>=thre&age<thre+10);
    thre=thre+10;
    ad_age_mae(i)=mean(abs(ad_gap_new_sort(pos)));
    ad_age_gap(i)=mean((ad_gap_new_sort(pos)));
    ad_age_gap_s(i)=std((ad_gap_new_sort(pos)));
end
site_hc=site';
for i=1:10
    pos=find(site_hc==i);
    ad_site_mae(i)=mean(abs(ad_gap_new_sort(pos)));
    ad_site_gap(i)=mean(ad_gap_new_sort(pos));
    ad_site_gap_s(i)=std(ad_gap_new_sort(pos));
end
save('/home/cxpang/matlab/code/8.brain_age/adjusted_all_gap_remove65sex_new.mat','ad_gap_new')
Idx=zeros(1065,1);
Idx(find(ad_gap_new_sort<0))=1;
Idx(find(ad_gap_new_sort>0))=2;
save('/home/cxpang/matlab/code/8.brain_age/index_hc_remove65sex.mat','Idx') %for validation
%%
x_all=[];
y_all=[];
id=[];
for i=1:20
x_all = [x_all;Prediction.Score{i}];
y_all= [y_all;Prediction.realScore{i}];
id=[id;Prediction.Origin_ID{i}];
end
p=polyfit(y_all,ad_gap_new,1);
yfit=polyval(p,y_all);%
plot(y_all,ad_gap_new,'bo',y_all,yfit,'k-','Markersize',3.5);
[r,p]=corr(y_all,ad_gap_new)

p=polyfit(y_all,hc_gap,1);
yfit=polyval(p,y_all);%
plot(y_all,hc_gap,'bo',y_all,yfit,'k-','Markersize',3.5);
[r,p]=corr(y_all,hc_gap)
%% correlation OF MAE between sites cross folds
load(strcat(data_dir,'/SVR_all_no_regress_fold20_remove65sex/Prediction.mat'))
load ('/HeLabData2/cxpang/DIDA/using_mat/site_remove65.mat')
load ('/HeLabData2/cxpang/DIDA/using_mat/glmvariable_remove65.mat')
hc=find(group==1);
site=site(hc);
mae_site_all=zeros(20,10);
x_all=[];
y_all=[];
id=[];
for i=1:20
    x_all = [x_all;Prediction.Score{i}];
    y_all= [y_all;Prediction.realScore{i}];
    id=[id;Prediction.Origin_ID{i}];
end
data = [id, x_all];
sorted_data = sortrows(data, 1);
x_all=sorted_data(:,2);
data = [id, y_all];
sorted_data = sortrows(data, 1);
y_all=sorted_data(:,2);

for i=1:20
    x = Prediction.Score{i};
    y= Prediction.realScore{i};
    id=Prediction.Origin_ID{i};
    site_hc=site(id);
    for j=1:1065
        tmp(j)=j;
    end
    id_diff=setdiff(tmp,id);
    x_train=x_all(id_diff);
    y_train=y_all(id_diff);
    gap=x_train-y_train;
    beta=glmfit(y_train,gap);
    ad_gap=x-y;
    yfit=glmval(beta,y,'identity');
    ad_gap_new=ad_gap-yfit; % corrected age gap
    
    for k=1:10
        pos=find(site_hc==k);
        ad_gap_new_site=ad_gap_new(pos);
        mae_site_all(i,k)=mean(abs(ad_gap_new_site));
    end
end
nan_rows = any(isnan(mae_site_all), 2);
mae_site_all = mae_site_all(~nan_rows, :);
[p,tbl,stats] = anova1(mae_site_all);%单因素一元方差分析
[c,~,~,gnames] = multcompare(stats);
%% all hc to mdd model
ResultantFolder = [data_dir,'/SVR_all_no_regress_fold20_remove65sex'];
load(strcat('/home/cxpang/matlab/code/8.brain_age/using_mat/covariate_hc_remove65.mat'))
covariate=[sex];
Subjects_Scores =age;
load(strcat(data_dir,'loading_hc_remove65.mat'))
Subjects_Data=loading;
Pre_Method = 'Normalize';
[slm,model,MeanValue,StandardDeviation] = SVR_all(covariate,Subjects_Data, Subjects_Scores,Pre_Method,0, ResultantFolder);
clear Subjects_Data
clear Subjects_Scores
save('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/hc_all_model_remove65sex.mat','model','MeanValue','StandardDeviation')
load(strcat('/home/cxpang/matlab/code/8.brain_age/using_mat/covariate_mdd_remove65.mat'))
load(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/brain_age/loading_mdd_remove65.mat'),'loading')
test_data=loading;
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
r_predic=corr(Predicted_Scores,test_score);
save('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/mdd_all_res_remove65sex.mat','Predicted_Scores','test_score')
%% mdd gap
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
load('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/mdd_all_res_remove65sex.mat','Predicted_Scores','test_score')
% compute the age gap
age_real=test_score;
age_predic=Predicted_Scores;
mdd_all_gap=age_predic-age_real;
% regress out age dependence using beta computed from the training set
my_corr=corr(age_predic,age_real)
mae=mean(abs(mdd_all_gap))
% adjusted_mdd_age=(age_predic-d)/k;
% adjusted_mdd_all_gap=adjusted_mdd_age-age_real;
% regress 
yfit=glmval(beta,age_real,'identity');
adjusted_mdd_all_gap=mdd_all_gap-yfit; % corrected age gap
mean(adjusted_mdd_all_gap)
std(adjusted_mdd_all_gap)
mean(mdd_all_gap)
std(mdd_all_gap)
save('/home/cxpang/matlab/code/8.brain_age/adjusted_mdd_all_gap_remove65sex.mat','adjusted_mdd_all_gap')
[h,p,~,stats]=ttest2(ad_gap_new,adjusted_mdd_all_gap);
%% mdd analyse
load('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/mdd_all_res_remove65sex.mat','Predicted_Scores','test_score')
load ('/HeLabData2/cxpang/DIDA/using_mat/site_remove65.mat')
load ('/HeLabData2/cxpang/DIDA/using_mat/glmvariable_remove65.mat')
mdd=find(group==2);
site_mdd=site(mdd);
age_mdd=age(mdd);
adjusted_all_gap_age=cell(6,1);
thre=10;
for i=1:6
    pos=find(age_mdd>=thre&age_mdd<thre+10);
    thre=thre+10;
    x_tmp=Predicted_Scores(pos);
    y_tmp=test_score(pos);
    adjusted_all_gap_age{i}=adjusted_mdd_all_gap(pos);
    age_ad_mae(i)=mean(abs(adjusted_mdd_all_gap(pos)));
    age_ad_gap(i)=mean((adjusted_mdd_all_gap(pos)));
    age_ad_gap_s(i)=std((adjusted_mdd_all_gap(pos)));
    
    %corr_age(i)=corr(x_tmp,y_tmp);
    mae_age(i)=mean(abs(x_tmp-y_tmp));
end


for i=1:10
    pos=find(site_mdd==i);
    x_tmp=Predicted_Scores(pos);
    y_tmp=test_score(pos);
    save(strcat('/home/cxpang/matlab/code/8.brain_age/plot/mdd_',num2str(i),'xy'),'x_tmp','y_tmp')
    site_ad_mae(i)=mean(abs(adjusted_mdd_all_gap(pos)));
    mae_site(i)=mean(abs(x_tmp-y_tmp));
    site_ad_gap(i)=mean((adjusted_mdd_all_gap(pos)));
    site_ad_gap_s(i)=std((adjusted_mdd_all_gap(pos)));
    mean_age(i)=mean(y_tmp);
    std_age(i)=std(y_tmp);
end