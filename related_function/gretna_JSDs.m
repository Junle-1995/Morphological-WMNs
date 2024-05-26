function JSDs = gretna_JSDs(Data, varargin)
%
%==========================================================================
% This function is used to calculate the Jensen每Shannon (JS) distance based
% similarity among all pairs of probability distribution functions in the
% Data. The JS distance is the square root of the JS divergence, which is
% based on the Kullback每Leibler divergence; it is symmetric and always a
% finite value.
%
% Syntax: function JSDs = gretna_JSDs(Data, varargin)
%
% Inputs:
%       Data:
%            M*N1 data array of probability distribution function with rows
%            denoting observations, and N1 denoting variables/regions. The
%            probability distribution functions can be otained with the
%            gretna_PDF.m function.
%
%   varargin:
%            M*N2 data array of probability distribution function with rows
%            denoting observations, and N2 denoting variables/regions. The
%            probability distribution functions can be otained with the
%            gretna_PDF.m function. NOTE, the number of rows must be the
%            same with that of Data.
% Output:
%       JSDs:
%            Pairwise Jensen每Shannon distance based similarity.
%            (N1*N1 for one argument; N1*N2 for two arguments)
%
% References:
% 1) Endres, D. M.; J. E. Schindelin (2003). "A new metric for probability
% distributions". IEEE Trans. Inf. Theory. 49 (7): 1858每1860.
% doi:10.1109/TIT.2003.813506.
% 2) Osterreicher, F.; I. Vajda (2003). "A new class of metric divergences
% on probability spaces and its statistical applications". Ann. Inst.
% Statist. Math. 55 (3): 639每653. doi:10.1007/BF02517812.
%
% Hao WANG, HZNU, Hangzhou, China, 2017/04/06, hall.wong@outlook.com
% Jinhui WANG, HZNU, Hangzhou, China, 2017/04/06, jinhui.wang.1982@gmail.com
%==========================================================================

if nargin < 1
    error('At least one argument is needed!');
end

Px = Data;

if nargin < 2
    Py = Data;
else
    Py = varargin{1};
    if size(Px,1) ~= size(Py,1)
        error('The number of rows must be equal between the two inputs!')
    end
end

% Specify the accuracy of the floating-point precision(+eps), avoiding
% imaginary number in the following calculation process.
logPx = log2(Px + eps);
logPy = log2(Py + eps);

PxlogPx = dot(Px, logPx)';
PylogPy = dot(Py, logPy);

JSD = zeros(size(Px, 2),size(Py, 2));

for i = 1:size(Px, 2)
    M = bsxfun(@plus, Px(:,i),Py)/2;
    logM = log2(M + eps);
    JSD(i,:) = bsxfun(@plus, -Px(:,i)'*logM, PxlogPx(i)) - (sum(Py.*logM,1)) + PylogPy;
end

JSD = JSD./2;
JSD = sqrt(JSD);
JSD = real(JSD);
JSDs = 1 - JSD;

if nargin < 2
    JSDs = (JSDs + JSDs')/2;
    JSDs(1:size(Data,2)+1:end) = 0;
end

return