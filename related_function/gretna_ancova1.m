function [F_ancova, P_ancova, DF_ancova, T_post, P_post, DF_post] = gretna_ancova1(Data, Cova)

%==========================================================================
% This function is used to perform one-way ancova on independent samples.
%
%
% Syntax: function [F_ancova, P_ancova, DF_ancova, T_post, P_post, DF_post] = gretna_ancova1(Data, Cova)
%
% Inputs:
%           Data:
%                N*1 cell with each cell containing a m*n array (N, number
%                of groups; m, number of observations; n, number of
%                variables).
%      Cova (opt):
%                N*1 cell with each cell containing a m*l array (N, number
%                of groups; m, number of observations; l, number of
%                covariates).
%
% Outputs:
%       F_ancova:
%                F statistics.
%       P_ancova:
%                Significance levels for F tests.
%      DF_ancova:
%                Degree of freedom for F tests.
%         T_post:
%                Post-hoc T statistics.
%         P_post:
%                Significance levels for T tests.
%        DF_post:
%                Degree of freedom for T tests.
%
%
% Jinhui WANG, IBRR, SCNU, Guangzhou, 2019/11/11, jinhui.wang.1982@gmail.com
%==========================================================================

Ngroup = size(Data,1);

Ind_group = [];
Nsub = zeros(Ngroup,1);

for igroup = 1:Ngroup
    Nsub(igroup,1) = size(Data{igroup},1);
    Ind_group = cat(1, Ind_group, ones(size(Data{igroup},1),1)*igroup);
end

Ind_group = dummyvar(Ind_group);

Covariate = [];
if nargin == 2
    Covariate = cell2mat(Cova);
end

Predictors = [Ind_group Covariate];

Data = cell2mat(Data);

Nvar = size(Data,2);
F_ancova = zeros(1,Nvar);
P_ancova = zeros(1,Nvar);
DF_ancova.group = zeros(1,Nvar);
DF_ancova.error = zeros(1,Nvar);
T_post = zeros(Ngroup*(Ngroup-1)/2,Nvar);
P_post = zeros(Ngroup*(Ngroup-1)/2,Nvar);
DF_post = zeros(Ngroup*(Ngroup-1)/2,Nvar);

for ivar = 1:Nvar
    
    Ind_nonan_ivar = ~isnan(Data(:,ivar)) & ~isinf(Data(:,ivar));
    
    if prod(sum(Ind_group(Ind_nonan_ivar,:))) == 0
        error('Error. \nAll data are NANs for %dth variable in at least one group!',ivar);
    end
    
    % unrestricted model
    DF_ancova.group(ivar) = Ngroup - 1;
    DF_ancova.error(ivar) = sum(Ind_nonan_ivar) - size(Predictors, 2);
    
    [Q,R] = qr(Predictors(Ind_nonan_ivar,:),0);
    beta = R \ (Q'*Data(Ind_nonan_ivar,ivar));
    
    r = Data(Ind_nonan_ivar,ivar) - Predictors(Ind_nonan_ivar,:)*beta;
    SSE = sum(r.^2);
    
    % restricted model
    Predictors_res = ones(sum(Ind_nonan_ivar),1);
    if nargin == 2
        Predictors_res = [ones(sum(Ind_nonan_ivar),1) Covariate(Ind_nonan_ivar,:)];
    end
    
    beta_res = (Predictors_res'*Predictors_res)^(-1) * Predictors_res'*Data(Ind_nonan_ivar,ivar);
    r_res = Data(Ind_nonan_ivar,ivar) - Predictors_res*beta_res;
    SSE_res = sum(r_res.^2);
    
    % F/P value
    F_ancova(ivar) = ((SSE_res-SSE)/DF_ancova.group(ivar))/(SSE/DF_ancova.error(ivar));
    P_ancova(ivar) = (1 - fcdf(F_ancova(ivar), DF_ancova.group(ivar), DF_ancova.error(ivar)));
    
    % post-hoc
    Num = Ngroup-1:-1:1;
    Index_start = 1;
    Contrast_post = zeros(Ngroup*(Ngroup-1)/2,size(Predictors, 2));
    
    for icon = 1:Ngroup-1
        
        Contrast_post(Index_start:sum(Num(1:icon)),icon) = 1;
        Contrast_post(Index_start:sum(Num(1:icon)),icon+1:Ngroup) = -eye(Ngroup-icon);
        
        Index_start = Index_start + sum(Num(icon));
        
    end
    
    std_e = sqrt(SSE/(sum(Ind_nonan_ivar)-size(Predictors,2)));
    d = sqrt(Contrast_post*(Predictors(Ind_nonan_ivar,:)'*Predictors(Ind_nonan_ivar,:))^(-1)*Contrast_post');
    
    T_post(:,ivar) = (Contrast_post*beta)./(std_e*diag(d));
    P_post(:,ivar) = 2*(tcdf(-abs(T_post(:,ivar)), DF_ancova.error(ivar)));
    DF_post(:,ivar) = DF_ancova.error(ivar);
    
end

return