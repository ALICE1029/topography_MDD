function [Prediction] = pca_svr(Subjects_Data, Subjects_Scores,covariate,FoldQuantity, Pre_Method, C_Range, ResultantFolder)

if ~exist(ResultantFolder, 'dir')
    mkdir(ResultantFolder);
end

[Subjects_Quantity, Feature_Quantity] = size(Subjects_Data);
%Split into N folds according to the behavioral scores
EachPart_Quantity = fix(Subjects_Quantity / FoldQuantity);
[~, SortedID] = sort(Subjects_Scores);
for j = 1:FoldQuantity
    Origin_ID{j} = SortedID([j : FoldQuantity : Subjects_Quantity]);
end

w_Brain = cell(FoldQuantity, 1);
Mask = cell(FoldQuantity, 1);

for j = 1:FoldQuantity
    
    disp(['The ' num2str(j) ' fold!']);
    
    Training_data = Subjects_Data;
    Training_scores = Subjects_Scores;
    
    % Select training data and testing data
    test_data = Training_data(Origin_ID{j}, :);
    test_score = Training_scores(Origin_ID{j});
    co_test=covariate(Origin_ID{j}, :);
    Training_data(Origin_ID{j}, :) = [];
    Training_scores(Origin_ID{j}) = [];
    co_train=covariate;
    co_train(Origin_ID{j},:)=[];
    %regress out covariate
    Covariates_training = covariate;
    Covariates_training(Origin_ID{j}, :) = [];
    
    b=zeros(size(Training_data,2),size(Covariates_training,2)+1);
    for i=1:size(Training_data,2)
        [Training_data(:,i), b(i,:), stats] = regress_out(Training_data(:,i), co_train);
    end
    y_regress=zeros(length(test_score),size(Training_data,2));
    b0=zeros(length(test_score),size(Training_data,2));
    for i=1:length(test_score)
        for s=1:size(Training_data,2)
            b0(i,s)=b(s,1);
        end
    end
    for i=1:size(Training_data,2)
        y_regress(:,i)=b0(:,i)+(b(i,2)*co_test(:,1)')';
    end
    test_data=test_data-y_regress;
    
    C_Optimal=1;  
    if strcmp(Pre_Method, 'Normalize')
        % Normalizing
        MeanValue = mean(Training_data);
        StandardDeviation = std(Training_data);
        [~, columns_quantity] = size(Training_data);
        for k = 1:columns_quantity
            Training_data(:, k) = (Training_data(:, k) - MeanValue(k)) / StandardDeviation(k);
        end
    elseif strcmp(Pre_Method, 'Scale')
        % Scaling to [0 1]
        MinValue = min(Training_data);
        MaxValue = max(Training_data);
        [~, columns_quantity] = size(Training_data);
        for k = 1:columns_quantity
            Training_data(:, k) = (Training_data(:, k) - MinValue(k)) / (MaxValue(k) - MinValue(k));
        end
    end
    Training_data(isnan(Training_data))=0;
    
    % SVR training
    Training_data_final = double(Training_data);
    model = svmtrain(Training_scores, Training_data_final, ['-s 3 -t 0 -c ' num2str(C_Optimal)]);
    % Normalize test data
    if strcmp(Pre_Method, 'Normalize')
        % Normalizing
        MeanValue_New = repmat(MeanValue, length(test_score), 1);
        StandardDeviation_New = repmat(StandardDeviation, length(test_score), 1);
        test_data = (test_data - MeanValue_New) ./ StandardDeviation_New;
    elseif strcmp(Pre_Method, 'Scale')
        % Scale
        MaxValue_New = repmat(MaxValue, length(test_score), 1);
        MinValue_New = repmat(MinValue, length(test_score), 1);
        test_data = (test_data - MinValue_New) ./ (MaxValue_New - MinValue_New);
    end
    test_data(isnan(test_data))=0;   
    test_data_final = double(test_data);
    % Predict test data
    [Predicted_Scores, ~, ~] = svmpredict(test_score, test_data_final, model);
    Prediction.Origin_ID{j} = Origin_ID{j};
    Prediction.Score{j} = Predicted_Scores;
    Prediction.realScore{j} = test_score;
    Prediction.Corr(j) = corr(Predicted_Scores, test_score)
    Prediction.MAE(j) = mean(abs(Predicted_Scores - test_score));
    Prediction.C_Optimal(j) = C_Optimal;
end

Prediction.Mean_Corr = mean(Prediction.Corr);
Prediction.Mean_MAE = mean(Prediction.MAE);
Prediction.weight = w_Brain;
Prediction.mask =Mask;
%
save([ResultantFolder filesep 'Prediction.mat'], 'Prediction');

end
