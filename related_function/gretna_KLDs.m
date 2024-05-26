function KLDs = gretna_KLDs(Data, varargin)

%==========================================================================
% This function is used to calculate KL divergence based similarity among
% all pairs of probability distribution functions in the Data.
%
% Syntax: function KLDs = gretna_KLDs(Data, varargin)
%
% Inputs:
%       Data:
%            M*N1 data array of probability distribution function with rows
%            denoting observations, and columns denoting variables. The
%            probability distribution functions can be otained with the
%            gretna_PDF.m function.
%
%   varargin:
%            M*N2 data array of probability distribution function with rows
%            denoting observations, and columns denoting variables. The
%            probability distribution functions can be otained with the
%            gretna_PDF.m function. NOTE, the number of rows must be the
%            same with that of Data.
% Output:
%       KLDs:
%            Pairwise KL divergence based similarity.
%            (N1*N1 for one argument; N1*N2 for two arguments)
%
% References:
% 1) Kong X, Liu Z, Huang L, et al. Mapping Individual Brain Networks Using
%    Statistical Similarity in Regional Morphology from MRI[J]. PloS one,
%    2015, 10(11):e0141840.
% 2) Wang Hao, Jin Xiaoqing,Zhang Ye and Wang Jinhui. Single-subject
%    morphological brain networks: connectivity mapping, topological
%    characterization and test¨Cretest reliability.  Brain and Behavior,
%    2016, 6(4):e00448.
%
% Hao WANG, HZNU, Hangzhou, China, 2017/01/13, hall.wong@outlook.com
% Jinhui WANG, HZNU, Hangzhou, China, 2017/01/13, jinhui.Wang.1982@gmail.com
%==========================================================================

if nargin < 1
    error('At least one argument is needed!');
end

Px = Data;
logPx = log(Px);
PxlogPx = dot(Px,logPx)';

if nargin == 1
    KLD = bsxfun(@plus,-Px'*logPx,PxlogPx);
    KLD = KLD + KLD';
    
    KLDs = exp(-KLD);
    KLDs(1:size(Data,2)+1:end) = 0;
    KLDs = (KLDs + KLDs')/2;
else
    Py = varargin{1};
    
    if size(Px,1) ~= size(Py,1)
        error('The number of rows must be equal between the two inputs!')
    end
    
    logPy = log(Py);
    PylogPy = dot(Py,logPy);
    
    KLD = bsxfun(@plus,-Px'*logPy,PxlogPx)+ bsxfun(@plus,-(Py'*logPx)',PylogPy);
    KLDs = exp(-KLD);
    KLDs = (KLDs + KLDs')/2;
end

return