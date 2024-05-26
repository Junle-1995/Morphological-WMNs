function [A,W] = gretna_R2b_MST(Original,Stype,Ttype,Thr)

%==========================================================================
% This function is used to threshold a weighted matrix into a binary matrix
% and weighted matrix by identifying the minimum spanning tree subgraph
% first, followed by adding new connections in descending order of
% connectivity strength among the left elements.
%
%
% Syntax: function [A,W] = gretna_R2b_MST(Original,Stype,Ttype,Thr)
%
% Inputs:
%        Original:
%                 The symmetric weighted matrix.
%          Stype:
%                'pos' for only positive elements;
%                'neg' for only negative elements using absolute values;
%                'abs' for all elements using absolute values.
%          Ttype:
%                'r' for correlation threshold;
%                's' for sparsity threshold;
%                'k' for edge threshold.
%            Thr:
%                 Correlation value if 'r';
%                 Sparsity value    if 's';
%                 Number of edges   if 'k'.
%
% Outputs:
%             A:
%               The resulting binary matrix.
%             W:
%               The resulting weighted matrix.
%
% Junle Li,    IBRR, SCNU, GuangZhou, 2019/10/20, lijunle.1995@gmail.com
% Jinhui WANG, IBRR, SCNU, Guangzhou, 2019/10/20, jinhui.wang.1982@gmail.com
%==========================================================================

Original = Original - diag(diag(Original));

if strcmpi(Stype,'pos')
    Original(Original < 0) = 0;
elseif strcmpi(Stype,'neg')
    Original(Original > 0) = 0;
    Original = abs(Original);
elseif strcmpi(Stype,'abs')
    Original = abs(Original);
else
    error('Stype must be ''pos'',''neg'' or ''abs''!');
end

N = length(Original);
A = zeros(N);
W = Original;

Weight = 1./Original;
Weight(Original == 0) = 0;
G = graph(Weight);
T = minspantree(G);
MST = T.Edges.EndNodes;
MST_n = size(MST,1);

for iMST = 1:MST_n
    A(MST(iMST,1),MST(iMST,2)) = 1;
    A(MST(iMST,2),MST(iMST,1)) = 1;
    Original(MST(iMST,1),MST(iMST,2)) = 0;
    Original(MST(iMST,2),MST(iMST,1)) = 0;
end

Vec = sort(Original(:),'descend');

if lower(Ttype) == 'r'
    if Thr < 0
        error('The following limits should be met: r >= 0');
    end
    rThr = Thr;
    
elseif lower(Ttype) == 's'
    if Thr < 2*MST_n/N/(N-1) || Thr > (2*MST_n+sum(Vec ~= 0))/N/(N-1)
        error(['The following limits should be met: ' num2str(2*MST_n/N/(N-1)) ' <= s <= ' num2str((2*MST_n+sum(Vec ~= 0))/N/(N-1))]);
    end
    K = ceil(Thr*N*(N-1));
    if mod(K,2) ~= 0
        K = K + 1;
    end
    K = K - 2*MST_n;
    rThr = Vec(K+1);
    
elseif lower(Ttype) == 'k'
    if Thr < MST_n || Thr > MST_n+sum(Vec ~= 0)/2
        error(['The following limits should be met: ' num2str(MST_n) ' <= k <= ' num2str(MST_n+sum(Vec ~= 0)/2 )])
    end
    K = 2*(Thr-MST_n);
    rThr = Vec(K+1);
else
    error('Ttype must be ''r'', ''s'' or ''k''');
end

A(Original > rThr) = 1;
W = W.*A;

if lower(Ttype) ~= 'r' && K ~= 0
    if rThr == Vec(K)
        warning('The number of edges in the resultant network is not exactly equal to the specificed because multiple edges share the same weight of critical threshold!');
    end
end

return