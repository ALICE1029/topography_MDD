function [Prediction] = pca_svr_noco(Subjects_Data, Subjects_Scores,FoldQuantity, Pre_Method, C_Range, Netnum, Netfea,Voxelfe,Network_fs_num,randNet,Weight_Flag, Permutation_Flag, ResultantFolder)

if nargin >= 8
    if ~exist(ResultantFolder, 'dir')
        mkdir(ResultantFolder);
    end
end
[Subjects_Quantity, Feature_Quantity] = size(Subjects_Data);
%[Subjects_Quantity, ~] = size(Subjects_Data);
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
    Training_data(Origin_ID{j}, :) = [];
    Training_scores(Origin_ID{j}) = [];
    
    
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
%     w_Brain{j}=zeros(1, size(Training_data,2));
%     for t = 1 : model.totalSV
%         w_Brain{j} = w_Brain{j} + model.sv_coef(t) * model.SVs(t, :);
%     end
%     w_Brain{j} = w_Brain{j} / norm(w_Brain{j});
    %        w(:,j) = model.SVs' * model.sv_coef;
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
    %Predicted_Scores=test_data_final*model(1)+model(2);
    [Predicted_Scores, ~, ~] = svmpredict(test_score, test_data_final, model);
    
    %rescale
  %  Predicted_Scores=rescale(Predicted_Scores,13,30);
    
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

if nargin >= 8
    save([ResultantFolder filesep 'Prediction.mat'], 'Prediction');
    disp(['The correlation is ' num2str(Prediction.Mean_Corr)]);
    disp(['The MAE is ' num2str(Prediction.Mean_MAE)]);
    % Calculating w
    if Weight_Flag
        %  W_Calculate_SVR_CSelect(Subjects_Data, Subjects_Scores, Pre_Method, C_Range, ResultantFolder,FoldQuantity);
    end
end
