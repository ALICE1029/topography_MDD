function [model,MeanValue,StandardDeviation] = SVR_all(Subjects_Data, Subjects_Scores,Pre_Method,Weight_Flag, ResultantFolder)

[Subjects_Quantity, Feature_Quantity] = size(Subjects_Data);
Training_data = Subjects_Data;
Training_scores = Subjects_Scores;
%Training_scores=Training_scores(randperm(length(Training_scores)))

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
Training_data_final = double(Training_data);
C_Optimal=1;
model = svmtrain(Training_scores, Training_data_final, ['-s 3 -t 0 -c ' num2str(C_Optimal)]);
% [Predicted_Scores, ~, ~] = svmpredict(Training_scores, Training_data_final, model);
% % r_predic=corr(Predicted_Scores,Training_scores)