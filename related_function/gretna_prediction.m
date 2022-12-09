function Results = gretna_prediction(Predictors, Responses, K_fold, C_type, C_thr, R_type, Ncomp, Lambda)

%==========================================================================
% This function is used to predict different response variables using
% multiple predictive variables based on different models.
%
% Syntax: function Results = gretna_prediction(Predictors, Responses, K_fold, C_type, C_thr, R_type, Ncomp, Lambda)
%
% Inputs:
% Predictors:
%          N*P1 data array with N being observations and P1 being variables.
% Responses:
%          N*P2 data array with N being observations and P2 being variables.
%   K_fold:
%          A positive integer indicating the number of folds for
%          cross-validation.
%   C_type:
%          The methods used for dimensionality reduction.
%            'No':       No dimensionality reduction.
%            'PCA':      Principal component analysis.
%            'Pearson':  Pearson correlation.
%            'Spearman': Spearman correlation.
%    C_thr:
%          [] if C_type = 'No'.
%          The minimum variance explained if C_type = 'PCA'
%          (>0 and <=100).
%          The p-value for significant correlations if C_type = 'Pearson'
%          (>0 and <1).
%          The p-value for significant correlations if C_type = 'Spearman'
%          (>0 and <1).
%   R_type:
%          The regression methods used for the prediction.
%            'glm':   General linear model.
%            'svr':   Support vector regression.
%            'pls':   Partial least-squares regression (with parameter ncomp).
%            'ridge': Ridge regression (with parameter Lambda).
%            'lasso': Lasso regression (with parameter Lambda).
%    Ncomp:
%          [] if R_type = 'glm';
%          [] if R_type = 'svr';
%          A vector of nonnegative values if R_type = 'pls' (default = 10);
%          [] if R_type = 'ridge';
%          [] if R_type = 'lasso'.
%   Lambda:
%          [] if R_type = 'glm';
%          [] if R_type = 'svr';
%          [] if R_type = 'pls';
%          A vector of nonnegative values if R_type = 'ridge'
%          (default = 1:10^4) (Joshua Sarfaty Siegel et al., 2016 PNAS);
%          A vector of nonnegative values if R_type = 'lasso'
%          (default = 1:10^4) (Joshua Sarfaty Siegel et al., 2016 PNAS).
%
% Outputs:
% Results.Train_size:
%          The number of subjects in the train data of each fold.
% Results.Explained_variance:
%          The percentage of total variance explained by the selected
%          principal components if C_type = 'PCA'.
% Results.R:
%          The Pearson correlation coefficients between actual responses
%          and predicted responses.
% Results.P:
%          The significance levels of Results.R.
% Results.Response_predicted
%          The predicted response scores.
% Results.Beta_fold
%          The beta coeifficent of each predictive variable of each fold.
% Results.Consensus_numbers:
%          The number of each predictive variable that is selected as a
%          valid feature across all folds.
% Results.Consensus_weights:
%          The mean weigths of each predictive variable across all folds,
%          indicating the extent to which the predictive variables are
%          positively or negatively related to the response variables.
%
% Zhen Li,     IBRR, SCNU, GuangZhou, 2021/1/15, zhen.li.dada@gmail.com
% Jinhui Wang, IBRR, SCNU, GuangZhou, 2021/1/15, jinhui.wang.1982@gmail.com
%==========================================================================

if nargin < 6
    error('At least 6 arguments are required!');
end

if nargin > 8
    error('At most 7 arguments are allowed!');
end

if nargin == 6
    if strcmpi(R_type,'pls')
        Ncomp = 10;
    end
    
    if strcmpi(R_type,'ridge') || strcmpi(R_type,'lasso')
        Lambda = 1:10^4;
    end
end

[Nsub_pre, Nvar_pre] = size(Predictors);
[Nsub_res, Nvar_res] = size(Responses);

if Nsub_pre ~= Nsub_res
    error('The number of observations are not equal between Predictors and Responses!');
end

if K_fold <= 1 || mod(K_fold,1) ~= 0
    error('K_fold must be a >1 positive integer!');
end

if ~(strcmpi(C_type,'No') || strcmpi(C_type,'pca') || strcmpi(C_type,'pearson') || strcmpi(C_type,'spearman'))
    error('Unrecognized input for C_type!');
elseif strcmpi(C_type,'pca')
    if C_thr <= 0 || C_thr > 100
        error('The range of C_thr should be (0 100] when C_type is PCA!');
    end
end

if ~(strcmpi(R_type,'glm') || strcmpi(R_type,'svr') || strcmpi(R_type,'pls') || strcmpi(R_type,'ridge') || strcmpi(R_type,'lasso'))
    error('Unrecognized input for R_type!');
end

%% prediction
Explained_variance       = zeros(K_fold,1);
Responses_predicted_test = zeros(size(Responses));
Beta_fold                = cell(Nvar_res,1);
Train_size               = zeros(Nvar_res,K_fold);

for ivar = 1:Nvar_res
    
%     fprintf('Predicting the %uth variable\n', ivar);
    
    Beta_fold{ivar,1}  = zeros(Nvar_pre,K_fold);
    Cvpar              = cvpartition(zeros(Nsub_pre,1),'KFold',K_fold);
    Train_size(ivar,:) = Cvpar.TrainSize;
    
    for ifold = 1:K_fold
        Idx_train        = Cvpar.training(ifold);
        Idx_test         = ~Idx_train;
        Predictors_train = Predictors(Idx_train,:);
        Responses_train  = Responses(Idx_train,ivar);
        Predictors_test  = Predictors(Idx_test,:);
        
        % feature standarization
        mean_Predictors_train = mean(Predictors_train);
        std_Predictors_train  = std(Predictors_train);
        Predictors_train = (Predictors_train - repmat(mean_Predictors_train,size(Predictors_train,1),1))./repmat(std_Predictors_train,size(Predictors_train,1),1);
        Predictors_test  = (Predictors_test - repmat(mean_Predictors_train,size(Predictors_test,1),1))./repmat(std_Predictors_train,size(Predictors_test,1),1);
        
        switch lower(C_type)
            case 'pca'
                [Coeff_train, Score_train, ~, ~, Explained_train, mu] = pca(Predictors_train);
                Score_test                = (Predictors_test - repmat(mu,size(Predictors_test,1),1)) * Coeff_train;
                Sum_explained_train       = cumsum(Explained_train);
                Idx_com                   = find(Sum_explained_train >= C_thr,1);
                Explained_variance(ifold) = sum(Explained_train(1:Idx_com));
                Predictors_train          = Score_train(:,1:Idx_com);
                Predictors_test           = Score_test(:,1:Idx_com);
            case 'pearson'
                [~,p]                     = corr(Predictors_train,Responses_train);
                Index_feature             = p < C_thr;
                Predictors_train          = Predictors_train(:,Index_feature);
                Predictors_test           = Predictors_test(:,Index_feature);
            case 'spearman'
                [~,p]                     = corr(Predictors_train,Responses_train,'type','Spearman');
                Index_feature             = p < C_thr;
                Predictors_train          = Predictors_train(:,Index_feature);
                Predictors_test           = Predictors_test(:,Index_feature);
        end
        
        if isempty(Predictors_train)
            Responses_predicted_test(Idx_test,ivar) = nan;
            Beta_fold{ivar,1}(:,ifold)              = 0;
        else
            switch lower(R_type)
                case 'glm'
                    Stats                                   = regstats(Responses_train,Predictors_train,'linear','beta');
                    Beta_2_end                              = Stats.beta(2:end);
                    Responses_predicted_test(Idx_test,ivar) = Stats.beta(1) + Predictors_test * Beta_2_end;
                case 'svr'
                    Model                                   = fitrsvm(Predictors_train,Responses_train);
                    Beta_2_end                              = Model.Beta;
                    Responses_predicted_test(Idx_test,ivar) = predict(Model,Predictors_test);
                case 'pls'
                    [~,~,~,~,Beta_2_end]                    = plsregress(Predictors_train,Responses_train,Ncomp);
                    Responses_predicted_test(Idx_test,ivar) = [ones(size(Predictors_test,1),1) Predictors_test]*Beta_2_end;
                    Beta_2_end                              = Beta_2_end(2:end);
                case 'ridge'
                    Beta                                    = ridge(Responses_train,Predictors_train,Lambda,0);
                    Responses_predicted_train               = repmat(Beta(1,:),Cvpar.TrainSize(ifold),1) + Predictors_train * Beta(2:end,:);
                    MSE                                     = mean((repmat(Responses_train,1,length(Lambda))- Responses_predicted_train).^2,1)';
                    Beta_final                              = Beta(:,find(MSE == min(MSE),1));
                    Beta_2_end                              = Beta_final(2:end);
                    Responses_predicted_test(Idx_test,ivar) = Beta_final(1) + Predictors_test * Beta_2_end;
                case 'lasso'
                    [Beta,Fitinfo]                          = lasso(Predictors_train,Responses_train,'lambda',Lambda);
                    Responses_predicted_train               = repmat(Fitinfo.Intercept(1,:),Cvpar.TrainSize(ifold),1) + Predictors_train * Beta;
                    MSE                                     = mean((repmat(Responses_train,1,length(Lambda))- Responses_predicted_train).^2,1)';
                    Idx_beta                                = find(MSE == min(MSE),1);
                    Beta_2_end                              = Beta(:,Idx_beta);
                    Responses_predicted_test(Idx_test,ivar) = Fitinfo.Intercept(Idx_beta) + Predictors_test * Beta_2_end;
            end
            
            switch lower(C_type)
                case 'no'
                    Beta_fold{ivar,1}(:,ifold)             = Beta_2_end;
                case 'pca'
                    Beta_fold{ivar,1}(:,ifold)             = Coeff_train(:,1:Idx_com) * Beta_2_end;
                case {'pearson', 'spearman'}
                    Beta_fold{ivar,1}(Index_feature,ifold) = Beta_2_end;
            end
        end
    end
    
%     fprintf('Predicting the %uth variable ...... is done\n', ivar);
    
end

Consensus_numbers = zeros(Nvar_res,Nvar_pre);
Consensus_weights = zeros(Nvar_res,Nvar_pre);
for ivar = 1:Nvar_res
    Consensus_numbers(ivar,:) = sum(Beta_fold{ivar,1} ~= 0,2);
    Consensus_weights(ivar,:) = mean(Beta_fold{ivar,1},2);
end

Results            = struct;
Results.Train_size = Train_size;

if strcmpi(C_type,'pca')
    Results.Explained_variance = Explained_variance;
end

R = zeros(Nvar_res,1);
P = zeros(Nvar_res,1);
for ivar = 1:Nvar_res
    if sum(~isnan(Responses_predicted_test(:,ivar))) >= 2
        [r,p]     = corr(Responses(:,ivar),Responses_predicted_test(:,ivar),'rows','pairwise');
        R(ivar,1) = r;
        P(ivar,1) = p;
    else
        warning('There are no valid predicted values for %uth response variable due to no features left in all folds!\n',ivar);
    end
end

Results.R                  = R;
Results.P                  = P;
Results.Response_predicted = Responses_predicted_test;
Results.Beta_fold          = Beta_fold;
Results.Consensus_numbers  = Consensus_numbers;
Results.Consensus_weights  = Consensus_weights;

return