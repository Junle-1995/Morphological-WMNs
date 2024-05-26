function [PDF, Tsamp] = gretna_PDF(Data, Npoints, Type)

%==========================================================================
% This function is used to calculate the probability distribution function
% of 1-D variable based on an optimized kernel density estimate.
%
% Syntax: function [PDF, Tsamp] = gretna_PDF(Data, Npoints, Type)
%
% Inputs:
%       Data:
%            A m-by-n data array with rows denoting observations and columns
%            denoting variables.
%    Npoints:
%            A scalar for the number of sampling points (e.g., 2^8).
%       Type:
%            'ksdensity': ksdensity works best for continuously
%                         distributed samples (matlab function, default).
%            'tde':       Topological density estimation (Huntsman, 2017).
%            'kde':       Reliable and extremely fast kernel density
%                         estimation (Botev et al., 2010).
%
% Outputs:
%        PDF:
%            The resulant probability distribution function.
%      Tsamp:
%            Optimized points at which the PDF are estimated.
%
% Hao WANG, CCBD, HZNU, Hangzhou, 2015/04/29, hall.wong@outlook.com
% Jinhui WANG, SCNU, Guangzhou,2018/11/27, jinhui.wang.1982@gmail.com
%==========================================================================

if nargin < 1
    error('At least one argument is needed!'); end

if nargin == 1
    Npoints = 2^8; Type = 'ksdensity'; end

if nargin == 2
    Type = 'ksdensity'; end

if nargin > 3
    error('At most three arguments are permitted!'); end

[~, N] = size(Data);
PDF    = zeros(Npoints, N);
Tsamp  = zeros(Npoints, N);

for i = 1:N
    
    Data_tem = Data(:,i);
    Data_tem(isnan(Data_tem)) = [];
    
    if isempty(Data_tem)
        error('No valid data for the inputs, all are NAN!');
    else
        % exclude valus 2.58 standard deviation larger or smaller than the
        % mean (i.e., 99% confidence interval)
        mean_tem = mean(Data_tem); std_tem = std(Data_tem);
        Data_tem(Data_tem > mean_tem + 2.58*std_tem | Data_tem < mean_tem - 2.58*std_tem) = [];
        
        switch lower(Type)
            
            case 'ksdensity'
                [KDE, t] = ksdensity(Data_tem, 'npoints', Npoints);
                KDE = KDE'; KDE(KDE <= 0) = eps; t = t';
                
            case 'tde'
                tde = tde1d(Data_tem, Npoints, 1);
                t   = tde.x;
                KDE = tde.y;
                
            case 'kde'
                % Default window parameter is optimal for normal distribution
                sig = median(abs(Data_tem-median(Data_tem)))/0.6745;
                
                if sig <= 0
                    sig = max(Data_tem)-min(Data_tem); end
                
                if sig > 0
                    Data_u = sig * (4/(3*length(Data_tem)))^0.2;
                else
                    Data_u = 1;
                end
                
                % 'Unbounded' if the density can extend over the whole real
                % line (it is the same with matlab funcion ksdensity)
                t = linspace(min(Data_tem)-3*Data_u, max(Data_tem)+3*Data_u, Npoints)';
                
                % Performing kernel density estimation based on a normal kernel function
                % with bandwidths adapted to data
                [~,KDE,~,~] = kde(Data_tem, Npoints, min(t), max(t));
                KDE(KDE <= 0) = eps;
                
            otherwise
                error('The inputted Type is not recognized, please check it!')
        end
        
        % Calculate probability distribution functions
        PDF(:,i) = KDE./sum(KDE);
        Tsamp(:,i) = t;
        
    end
end

return