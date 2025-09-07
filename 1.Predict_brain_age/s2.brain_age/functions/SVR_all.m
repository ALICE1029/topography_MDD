function [slm,model,MeanValue,StandardDeviation] = SVR_all(covariate,Subjects_Data, Subjects_Scores,Pre_Method,Weight_Flag, ResultantFolder)

[Subjects_Quantity, Feature_Quantity] = size(Subjects_Data);
Training_data = Subjects_Data;
Training_scores = Subjects_Scores;
%Training_scores=Training_scores(randperm(length(Training_scores)))
Covariates_training=covariate;
% 
% %regress out covariate
M = 1;
[Training_quantity, Covariates_quantity] = size(Covariates_training);
    for k = 1:Covariates_quantity
        M = M + term(Covariates_training(:, k));
    end
    slm = SurfStatLinMod(Training_data, M);
    
    Training_data = Training_data - repmat(slm.coef(1, :), Training_quantity, 1);
    for k = 1:Covariates_quantity
        Training_data = Training_data - ...
            repmat(Covariates_training(:, k), 1, Feature_Quantity) .* repmat(slm.coef(k + 1, :), Training_quantity, 1);
    end


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