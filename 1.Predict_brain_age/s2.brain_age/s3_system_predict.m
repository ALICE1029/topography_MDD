%% 18
voxel_num=45892;
network_num=18;
method='NGSR';
load('/HeLabData2/cxpang/DIDA/using_mat/glmvariable_remove65.mat','age','sex','mfd','group')
load(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/y_all_remove65.mat'),'y_all')
hc=find(group==1);
for i=1:18
    loading_hc{i}=y_all{i}(hc,:);
end
save(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/brain_age/loading_hc_18.mat'),'loading_hc')
mdd=find(group==2);
for i=1:18
    loading_mdd{i}=y_all{i}(mdd,:);
end
save(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/brain_age/loading_mdd_18.mat'),'loading_mdd')
data_dir='/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/';
%%
load(strcat(data_dir,'loading_hc_18.mat'))
load(strcat('/home/cxpang/matlab/code/8.brain_age/using_mat/covariate_hc_remove65.mat'))
Subjects_Scores =age;
for i=1:18
    Subjects_Data=loading_hc{i};
    ResultantFolder = strcat(data_dir,'/SVR_all_regress_sex_fold20/',num2str(i));
    covariate=[sex];
    FoldQuantity =20;
    Pre_Method = 'Normalize';
    C_Range =1
    Weight_Flag = 1;
    Permutation_Flag = 0;
    randNet = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18];
    [Prediction] = pca_svr(Subjects_Data, Subjects_Scores,covariate,FoldQuantity, Pre_Method, C_Range, network_num,voxel_num,voxel_num,network_num,randNet,Weight_Flag, Permutation_Flag, ResultantFolder);
    x=[];
    y=[];
    for j=1:20
        x = [x;Prediction.Score{j}];
        y= [y;Prediction.realScore{j}];
        
    end
    my_corr=corr(x,y)
    save(strcat(ResultantFolder,'/res.mat'),'x','y')
end
%% 18 HC GAP
ad_gap_network=cell(18,1);
for k=1:18
    ResultantFolder = strcat(data_dir,'/SVR_all_regress_sex_fold20',num2str(k));
    load(strcat(ResultantFolder,'/res.mat'),'x','y')
    
    load(strcat(ResultantFolder,'/Prediction.mat'))
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
    ad_gap_new_sort(:,k)=sorted_data(:,2);
    adj_mae(k)=mean(abs(ad_gap_new_sort(:,k)))%adj mae
    ad_gap_network{k}=ad_gap_new;
end
save(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/brain_age/ad_gap_18.mat'),'ad_gap_network')
%% draw scatter plot
load ('/HeLabData2/cxpang/DIDA/using_mat/site_remove65.mat')
load ('/HeLabData2/cxpang/DIDA/using_mat/glmvariable_remove65.mat')
hc=find(group==1);
site=site(hc);
for k=1:18
    k
    
    ResultantFolder = strcat(data_dir,'/SVR_all_regress_sex_fold20',num2str(k));
    load(strcat(ResultantFolder,'/Prediction.mat'))
    x=[];
    y=[];
    id=[];
    for i=1:20
        x = [x;Prediction.Score{i}];
        y= [y;Prediction.realScore{i}];
        id=[id;Prediction.Origin_ID{i}];
    end
    my_corr(k)=corr(x,y);
    mae(k)=mean(abs(x-y));
    for i=1:length(id)
        site_hc(i)=site(id(i));
    end
    for i=1:10
        pos=find(site_hc==i);
        x_tmp=x(pos);
        y_tmp=y(pos);
        save(strcat('/home/cxpang/matlab/code/8.brain_age/plot/',num2str(k),'_',num2str(i),'xysex'),'x_tmp','y_tmp')
    end
    
end
%% all sub 18
load(strcat(data_dir,'loading_hc_18.mat'))
load(strcat('/HeLabData2/cxpang/DIDA/',method,'/result_xy/brain_age/loading_mdd_18.mat'),'loading_mdd')
for i=1:18
    i
    ResultantFolder = [data_dir,'/SVR_all_no_regress_fold20_remove65sex',num2str(i)];
    load(strcat('/home/cxpang/matlab/code/8.brain_age/using_mat/covariate_hc_remove65.mat'))
    covariate=[sex];
    Subjects_Scores =age;
    
    Subjects_Data=loading_hc{i};
    Pre_Method = 'Normalize';
    [slm,model,MeanValue,StandardDeviation] = SVR_all(covariate,Subjects_Data, Subjects_Scores,Pre_Method,0, ResultantFolder);
    clear Subjects_Data
    save(strcat('/home/cxpang/matlab/code/8.brain_age/hc_all_model_',num2str(i),'.mat'),'slm','model')
    save(strcat('/home/cxpang/matlab/code/8.brain_age/hc_all_model_',num2str(i),'_norm.mat'),'MeanValue','StandardDeviation')
    clear Subjects_Scores
    
    load(strcat('/home/cxpang/matlab/code/8.brain_age/using_mat/covariate_mdd_remove65.mat'))
    test_data=loading_mdd{i};
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
    MeanValue_New = repmat(MeanValue, length(test_score), 1);
    StandardDeviation_New = repmat(StandardDeviation, length(test_score), 1);
    test_data_norm = (test_data - MeanValue_New) ./ StandardDeviation_New;
    test_data_final = double(test_data_norm);
    [Predicted_Scores, ~, ~] = svmpredict(test_score, test_data_final, model);  
    save(strcat('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/mdd_all_res_remove65sex_',num2str(i),'_.mat'),'Predicted_Scores','test_score')
end
%% 18 mdd gap
for i=1:18
    load(strcat(data_dir,'/SVR_all_regress_sex_fold20',num2str(i),'/Prediction.mat'))
    x=[];
    y=[];
    id=[];
    for k=1:20
        x = [x;Prediction.Score{k}];
        y= [y;Prediction.realScore{k}];
    end
    all_gap= x-y;
    beta=glmfit(y,all_gap);
    load(strcat('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/mdd_all_res_remove65sex_',num2str(i),'_.mat'),'Predicted_Scores','test_score')
    mdd_all_gap(:,i)=Predicted_Scores-test_score;
    yfit=glmval(beta,test_score,'identity');
    adjusted_mdd_network(:,i)=mdd_all_gap(:,i)-yfit; % corrected age gap
end

%%  compare
load('/home/cxpang/matlab/code/8.brain_age/using_mat/index_all_remove65sex.mat','Idx')
sub1_network_gap=adjusted_mdd_network(find(Idx==1),:);
sub2_network_gap=adjusted_mdd_network(find(Idx==2),:);
c_net=cell(18,1);

for i=1:18
    all_data = [ad_gap_new_sort(:,i);sub1_network_gap(:,i); sub2_network_gap(:,i)];
    group_labels = [repmat({'Group1'}, length(ad_gap_network{i}), 1); ...
        repmat({'Group2'}, length(sub1_network_gap(:,i)), 1); ...
        repmat({'Group3'}, length(sub2_network_gap(:,i)), 1)];
    %p = vartestn(all_data, group_labels, 'TestType', 'LeveneAbsolute', 'Display', 'off');
    [p_net(i), tbl, stats] = anova1(all_data, (group_labels), 'off');
    tbl{2, 5}
    [c, m, h, gnames_net] = multcompare(stats);
    c_net{i}=c;
    
end
q=mafdr(p_net,'BHFDR', true);

save(strcat('/HeLabData2/cxpang/DIDA/NGSR/result_xy/brain_age/compare_network.mat'),'p_net','c_net')
for i=1:18
    if(c_net{i}(1,6)>0.05)
        i
    end
        
end
for i=1:18
    if(c_net{i}(2,6)>0.05)
        i
    end
        
end
q=mafdr(p_net,'BHFDR', true);

save('/home/cxpang/matlab/network_gap.mat','sub1_network_gap','sub2_network_gap','hc_network_gap')