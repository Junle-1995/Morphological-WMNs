function KSDs = gretna_KSDs(Data)

%==========================================================================
% This function is used to calculate the Kolmogorov-Smirnov distance based
% similarity among all pairs of probability distribution functions in the
% Data.
%
% Syntax: function  KSDs = gretna_KSDs(Data)
%
% Input:
%        Data:
%            M*N data array (M, the number of samples; N, the number of
%            variables.
%
% Outputs:
%       KSDs:
%            Pairwise Kolmogorov-Smirnov distance based similarity (N*N).
%
% Example:
%        KSDs = gretna_KSD(rand(120,10))
%
% Hao WANG, HZNU, Hangzhou, China, 2016/05/11, hall.wong@outlook.com
% Jinhui WANG, HZNU, Hangzhou, China, 2016/05/11, jinhui.wang.1982@gmail.com
%==========================================================================

[~, N] = size(Data);
D = zeros(N);

for ireg = 1:N-1
    for jreg = ireg+1:N
        [~,~,D(ireg,jreg)] = kstest2(Data(:,ireg),Data(:,jreg),'Alpha',0.05);
    end
end

KSD = D + D';

KSDs = 1 - KSD;
KSDs(1:N+1:end) = 0;
KSDs = (KSDs + KSDs')/2;

return