function [AUC] = gretna_auc(Netpara,Delta)

%==========================================================================
% This function is used to calculate area under curve.
%
%
% Syntax: function [AUC] = gretna_auc(Netpara,Delta)
%
% Inputs:
%           Netpara:
%                   Network parameters (m*n, m = observations; n = variables).
%             Delta:
%                   The interval of the threshold.
%
% Output:
%               AUC:
%                   Area Under Curve.
%
% Jinhui WANG, NKLCNL, BNU, BeiJing, 2011/10/23, Jinhui.Wang.1982@gmail.com
%==========================================================================

Nvar = size(Netpara,2);
AUC = zeros(Nvar,1);

for ivar = 1:Nvar
    Para_ivar = Netpara(:,ivar);
    Para_ivar(isnan(Para_ivar)) = [];
    
    if length(Para_ivar) <= 1
        AUC(ivar,1) = NaN ;
    else
        AUC(ivar,1) = (sum(abs(Para_ivar)) - sum(abs(Para_ivar([1 end])))/2) * Delta;
    end
end

return