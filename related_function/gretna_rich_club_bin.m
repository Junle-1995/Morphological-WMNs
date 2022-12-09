function [rc] = gretna_rich_club_bin(A)

%==========================================================================
% This function is used to calculate rich club metric for a binary graph or
% network G.
%
%
% Syntax:  function [rc] = gretna_rich_club(A)
%
% Input:
%        A:
%                The adjencent matrix of G.
%
% Outputs:
%        rc.deg:
%                Degree threshold used to determine rich club coefficient of G.
%        rc.fi:
%                Rich club coefficient of G.
%
% References:
% 1. Martijn P. van den Heuvel and Olaf Sporns (2011) Rich-Club Organization
%    of the Human Connectome. The Journal of Neuroscience, 31:15775¨C15786.
%
% Jinhui WANG, CCBD, HZNU, HangZhou, 2013/04/25, Jinhui.Wang.1982@gmail.com
%==========================================================================

A = A - diag(diag(A));
A = abs(A);

N = length(A);
K  = sum(A);

Kthr = 1:N-1;
rc.deg = Kthr;

for i = 1:length(Kthr)
    
    if Kthr(i) >= max(K)
        rc.fi(1,i) = 0;
    else
        
        ind = find(K > Kthr(i));
        
        if length(ind) == 1
            rc.fi(1,i) = 0;
        else
            net = A(ind, ind);
            rc.fi(1,i) = sum(net(:))/length(ind)/(length(ind)-1);
        end
    end
end

return