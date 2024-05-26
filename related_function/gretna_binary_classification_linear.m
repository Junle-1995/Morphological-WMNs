function Results = gretna_binary_classification_linear(Features, Group_label, K_fold, Method, C_type, C_thr)

%==========================================================================
% This function is used to classify two groups based on linear models.
%
% Syntax: function Results = gretna_binary_classification_linear(Features, Group_label, K_fold, Method, C_type, C_thr)
%
% Inputs:
%       Features:
%                N*P data array with N being observations and P being
%                variables.
%    Group_label:
%                N*1 data array indicating group labels (0 for controls and
%                1 for patients).
%         K_fold:
%                A positive integer indicating the number of folds for
%                cross-validation.
%         Method:
%                A string indicating the classification method.
%                  'SVM':      Support vector machine (default).
%                  'Logistic': Logistic regression.
%         C_type:
%                A string indicating the method for feature reduction.
%                  'No':    No dimensionality reduction (default).
%                  'PCA':   Principal component analysis.
%                  'Ttest': Independent t test.
%          C_thr:
%                The parameters used for feature reduction.
%                  [] if C_type = 'No'.
%                  The minimum variance explained if C_type = 'PCA'
%                  (>0 and <=100, default = 80).
%                  The p-value for significant effects if C_type = 'Ttest'
%                  (>0 and <1, default = 0.05).
%
% Outputs:
% Results.Label_predicted:
%         The predicted labels.
% Results.Matrix:
%         The 2*2 results matrix of the classification, indicating:
%         [TN  FP
%          FN  TP].
% Results.Accuracy:
%         The accuracy of the classification.
% Results.Sensitivity:
%         The sensitivity of the classification.
% Results.Specificity:
%         The specificity of the classification.
% Results.Train_size:
%         The number of subjects in the train data of each fold.
% Results.Method:
%         The classification method.
% Results.C_type:
%         The method for feature reduction.
% Results.C_thr:
%         The parameters of the method for feature reduction.
% Results.Consensus_numbers:
%         The number of each feature that is selected as a valid feature
%         across all folds if C_type = 'Ttest'.
% Results.Beta_fold:
%         The beta values in each fold.
% Results.Consensus_weights:
%         The mean weigths of each feature across all folds.
% Results.Explained_variance:
%         The percentage of total variance explained by the selected
%         principal components in each fold if C_type = 'PCA'.
%
% Junle Li,    IBRR, SCNU, GuangZhou, 2021/1/19, lijunle.1995@gmail.com
% Jinhui Wang, IBRR, SCNU, GuangZhou, 2021/1/19, jinhui.wang.1982@gmail.com
%==========================================================================

%% Default setting
if nargin < 3
    error('At least 3 arguments are required!');
end

if nargin == 3
    Method = 'SVM';
    C_type = 'No';
    C_thr  = [];
end

if nargin == 4
    C_type = 'No';
    C_thr  = [];
end

if nargin == 5
    if strcmpi(C_type,'No')
        C_thr = [];
    elseif strcmpi(C_type,'PCA')
        C_thr = 80;
    elseif strcmpi(C_type,'Ttest')
        C_thr = 0.05;
    end
end

if nargin > 6
    error('At most 6 arguments are required!');
end

%% Examining inputs
[Nsub, Nvar] = size(Features);

if length(Group_label) ~= Nsub
    error('The number of observations are not equal between Features and Group_label!');
end

if length(unique(Group_label)) ~= 2
    error('The number of Group_label type must be 2!');
end

if K_fold <= 1 || mod(K_fold,1) ~= 0
    error('K_fold must be a >1 positive integer!');
end

if ~(strcmpi(Method,'svm') || strcmpi(Method,'logistic'))
    error('Unrecognized input for Method!');
end

if ~(strcmpi(C_type,'No') || strcmpi(C_type,'PCA') || strcmpi(C_type,'Ttest'))
    error('Unrecognized input for C_type!');
end

if strcmpi(C_type,'No') && ~isempty(C_thr)
    warning('The C_thr will be set to empty!');
end

if strcmpi(C_type,'PCA')
    if C_thr <= 0 || C_thr > 100
        error('The range of C_thr for PCA should be (0 100]!');
    end
end

if strcmpi(C_type,'Ttest')
    if C_thr <= 0 || C_thr >= 1
        error('The range of C_thr for T_test should be (0 1)!');
    end
end

%% Classification
Label_predict = zeros(Nsub,1);

Feature_fold = zeros(Nvar,K_fold);
Beta_fold    = zeros(Nvar,K_fold);
if strcmpi(C_type,'PCA')
    Explained_variance = zeros(K_fold,1);
end

Cvpar = cvpartition(Group_label,'KFold',K_fold);
for ifold = 1:K_fold
    
%     disp(['Classifying the ' num2str(ifold) 'th fold    |' datestr(clock)])
    
    Idx_train      = Cvpar.training(ifold);
    Idx_test       = ~Idx_train;
    Features_train = Features(Idx_train,:);
    Label_train    = Group_label(Idx_train);
    Features_test  = Features(Idx_test,:);
    
    % selection to feature
    if strcmpi(C_type,'No')
        Feature_fold(:,ifold) = 1;
    end
    
    if strcmpi(C_type,'PCA')
        [Coeff_train, Score_train, ~, ~, Explained_train, mu] = pca(Features_train);
        Score_test                                            = (Features_test - repmat(mu,size(Features_test,1),1)) * Coeff_train;
        sum_Explained_train       = cumsum(Explained_train);
        Idx                       = find(sum_Explained_train >= C_thr,1);
        Explained_variance(ifold) = sum_Explained_train(Idx);
        Features_train            = Score_train(:,1:Idx);
        Features_test             = Score_test(:,1:Idx);
        Feature_fold(:,ifold)     = 1;
    end
    
    if strcmpi(C_type,'Ttest')
        [~,p]                          = ttest2(Features_train(Label_train == 0,:),Features_train(Label_train == 1,:));
        Features_train                 = Features_train(:,p <= C_thr);
        Features_test                  = Features_test(:,p <= C_thr);
        Feature_fold(p <= C_thr,ifold) = 1;
        if sum(Feature_fold(:,ifold)) == 0
            error('C_thr is too strict for feature selection!');
        end
    end
    
    % Model construction
    Mdl = fitclinear(Features_train',Label_train,'ObservationsIn','columns','Learner',Method);
    
    if strcmpi(C_type,'PCA')
        Beta_fold(:,ifold) = Coeff_train(:,1:Idx) * Mdl.Beta;
    else
        Beta_fold(logical(Feature_fold(:,ifold)),ifold) = Mdl.Beta;
    end
    
    Label_predict(Idx_test) = predict(Mdl,Features_test);
    
%     disp(['Classifying the ' num2str(ifold) 'th fold ...... is done    |' datestr(clock)])
    
end

%% Results
Results_matrix      = zeros(2);
Results_matrix(1,1) = sum(Group_label == 0 & Label_predict == 0);
Results_matrix(1,2) = sum(Group_label == 0 & Label_predict == 1);
Results_matrix(2,1) = sum(Group_label == 1 & Label_predict == 0);
Results_matrix(2,2) = sum(Group_label == 1 & Label_predict == 1);

Results                        = struct;
Results.Label_predicted        = Label_predict;
Results.Matrix                 = Results_matrix;
Results.Accuracy               = (Results_matrix(1,1) + Results_matrix(2,2)) / Nsub;
Results.Sensitivity            = Results_matrix(2,2) / (Results_matrix(2,1) + Results_matrix(2,2));
Results.Specificity            = Results_matrix(1,1) / (Results_matrix(1,1) + Results_matrix(1,2));
Results.Train_size             = Cvpar.TrainSize;
Results.Method                 = Method;
Results.C_type                 = C_type;
Results.C_thr                  = C_thr;
if strcmpi(C_type,'Ttest')
    Results.Consensus_numbers  = sum(Feature_fold,2);
end
Results.Beta_fold              = Beta_fold;
Results.Consensus_weights      = sum(Beta_fold,2);
if strcmpi(C_type,'PCA')
    Results.Explained_variance = Explained_variance;
end

return