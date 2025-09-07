
%%
addpath(genpath('/home/cxpang/matlab/code/8.brain_age'));
method='NGSR'
voxel_num=45892;
network_num=18;
data_dir='/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/';
addpath(genpath('/home/cxpang/matlab/SVR'));
%%
load('/HeLabData2/cxpang/treatment_NGSR/CSU/result/regress_variable.mat','x')
sbj_num=size(x,2);
y=x;
clear x
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
load('/HeLabData2/cxpang/treatment_NGSR/CSU/predict/deltahamd.mat')
load('/HeLabData2/cxpang/treatment_NGSR/CSU/predict/covariate.mat')
load('/HeLabData2/cxpang/treatment_NGSR/CSU/predict/mfd.mat','mfd')
load('/HeLabData2/cxpang/treatment_NGSR/CSU/predict/group.mat','group')
%% gap of treatment group
for s=1:18
    s
    load(strcat('/home/cxpang/matlab/code/8.brain_age/hc_all_model_',num2str(s),'_.mat'),'slm','model')    
    load(strcat('/home/cxpang/matlab/code/8.brain_age/hc_all_model_',num2str(s),'_norm.mat'),'MeanValue','StandardDeviation')    
    load('/home/cxpang/matlab/treatment/CSU_baseline.mat')
    lo1=y_all{s}(pos,:);
    load('/home/cxpang/matlab/treatment/CSU_after.mat')
    lo2=y_all{s}(pos,:);
    test_data=lo1;
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
    save(strcat('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/csutreatment_res_remove65sex_baseline_',num2str(s),'_.mat'),'Predicted_Scores','test_score')
    test_data=lo2;
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
    save(strcat('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/csutreatment_res_remove65sex_after_',num2str(s),'_.mat'),'Predicted_Scores','test_score')
end
%%
load('/HeLabData2/cxpang/treatment_NGSR/CSU/predict/mfd_co.mat','mfd_base','mfd_after')
for s= 1:18
    s
    load(strcat('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/csutreatment_res_remove65sex_after_',num2str(s),'_.mat'),'Predicted_Scores','test_score')
    
    load(strcat(data_dir,'/SVR_all_regress_sex_fold20',num2str(s),'/Prediction.mat'))
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
    yfit=glmval(beta,age_real,'identity');
    adjusted_mdd_all_gap=mdd_all_gap-yfit; % corrected age gap
  save(strcat('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/csuadgap_after_',num2str(s),'.mat'),'adjusted_mdd_all_gap')
load(strcat('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/csutreatment_res_remove65sex_baseline_',num2str(s),'_.mat'),'Predicted_Scores','test_score')
    load(strcat(data_dir,'/SVR_all_regress_sex_fold20',num2str(s),'/Prediction.mat'))
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
    yfit=glmval(beta,age_real,'identity');
    adjusted_mdd_all_gap1=mdd_all_gap-yfit; % corrected age gap
    save(strcat('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/csuadgap_base_',num2str(s),'.mat'),'adjusted_mdd_all_gap1')
    des=[age,sex,mfd_base];
    [adjusted_mdd_all_gap1, b, stats] = regress_out(adjusted_mdd_all_gap1, des);
    des=[age,sex,mfd_after];
    [adjusted_mdd_all_gap, b, stats] = regress_out(adjusted_mdd_all_gap, des);
    data{s} =[(adjusted_mdd_all_gap1),(adjusted_mdd_all_gap)];
% filename =strcat(num2str(s), '_mfd_treatmentcsu.txt');
% writematrix( [(adjusted_mdd_all_gap1),(adjusted_mdd_all_gap)], filename);
load('/HeLabData2/cxpang/treatment_NGSR/CSU/predict/group.mat','group')
    tbl = array2table(data{s}, 'VariableNames', {'Day1', 'Day2'});
    tbl.Group = group';
    tbl.Group = categorical(tbl.Group);
    withinDesign = table([1 2]', 'VariableNames', {'Time'});
    rm = fitrm(tbl, 'Day1, Day2 ~ Group', 'WithinDesign', withinDesign);
    ranovatbl2{s} = ranova(rm);
   % disp(ranovatbl);
   mcomp_time{s} = multcompare(rm, 'Time', 'By', 'Group');
    %disp(mcomp_time);
end
for i=1:18
       p1(i)=table2array(ranovatbl2{i}('(Intercept):Time','pValue'));
       p2(i)=table2array(ranovatbl2{i}('Group:Time','pValue')); 
end
q1=mafdr(p1,'BHFDR', true);
q2=mafdr(p2,'BHFDR', true);
q3=mafdr(p3,'BHFDR', true);
load('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/hamd_base.mat')
load('/HeLabData2/cxpang/treatment_NGSR/CSU/predict/deltahamd.mat','deltahamd')
des = [age,sex];
rate=zeros(43,1);
rate=deltahamd./hamd;
rate=deltahamd;
[rate, b, stats] = regress_out(rate, des);
[deltahamd, b, stats] = regress_out(deltahamd, des);
delta1=zeros(18,length(find(group==2)));
delta2=zeros(18,length(find(group==1)));
for s=1:18
    delta1(s,:)=data{s}(find(group==2),2)-data{s}(find(group==2),1);
    deltahamd1=rate(find(group==2));
    [r1(s),p1(s)]=corr((delta1(s,:))',deltahamd1);
    
    delta2(s,:)=data{s}(find(group==1),2)-data{s}(find(group==1),1);
    deltahamd2=rate(find(group==1));
    [r2(s),p2(s)]=corr(delta2(s,:)',deltahamd2);
end
