function [Net_para, Node_para] = gretna_batch_networkanalysis_bin(Mat_path, File_filter, Thr1, Thr2, Delta, N_rand, Stype, Ttype)

%==========================================================================
% This function is used to calculate several global and nodal network
% metrics for both binary and weighted networks over a continuous threshold
% range for multiple subjects.
%
%
% Syntax: function [Net_para, Node_para] = gretna_batch_networkanalysis(Mat_path, File_filter, Thr1, Thr2, Delta, N_rand, Stype, Ttype)
%
% Inputs:
%          Mat_path:
%                   The directory where individual matrices are sorted.
%       File_filter:
%                   The prefix of those matrices.
%              Thr1:
%                   The lower bound of the threshold range.
%              Thr2:
%                   The upper bound of the threhsold range.
%             Delta:
%                   The interval of the threshold range.
%            N_rand:
%                   The number of random networks.
%             Stype:
%                   'pos' for only positive elements;
%                   'neg' for only negative elements using absolute values;
%                   'abs' for all elements using absolute values.
%             Ttype:
%                   'r' for correlation threshold;
%                   's' for sparsity threshold;
%                   'k' for edge threshold.
%
%
% Jinhui WANG, IBRR, SCNU, Guangzhou, 2019/10/31, jinhui.wang.1982@gmail.com
%==========================================================================

Thres = Thr1:Delta:Thr2;

cd(Mat_path)
Mats = ls([File_filter '*.mat']);

for isub = 1:size(Mats,1) % subjects
    
    fprintf('------------------------------------------------------------------------------------\n')
    fprintf('Calculating network parameters for %s\n', Mats(isub,:));
    fprintf('------------------------------------------------------------------------------------\n')
    
    Matrix = load(deblank(Mats(isub,:)));
    Field_name = fieldnames(Matrix);
    Matrix = Matrix.(Field_name{1});
    N_region = length(Matrix);
    
    for ithres = 1:length(Thres) % thresholds
        
        [Bin, ~] = gretna_R2b_MST(Matrix,Stype,Ttype,Thres(ithres));
        
        [Net_para.bin.cp(ithres,isub), ~] = gretna_node_clustcoeff(Bin);
        [Net_para.bin.lp(ithres,isub), ~] = gretna_node_shortestpathlength(Bin);

        [M_bin] = gretna_rich_club_bin(Bin);
        Net_para.bin.fi(ithres,:,isub) = M_bin.fi;

        Net_para.bin.c(:,:,ithres,isub) = expm(Bin);
        Net_para.bin.dis(:,:,ithres,isub) = gretna_distance(Bin);
        
        [~, Node_para.bin.deg(ithres,isub,:)] = gretna_node_degree(Bin);
        
        Qrand_bin = zeros(1,length(Bin)-1,N_rand);
        for irand = 1:N_rand % random networks
            
            Bin_rand = gretna_gen_random_network1(Bin);
            
            [Net_para.bin.cprand(ithres,isub,irand), ~] = gretna_node_clustcoeff(Bin_rand);
            [Net_para.bin.lprand(ithres,isub,irand), ~] = gretna_node_shortestpathlength(Bin_rand);

            [M_bin] = gretna_rich_club_bin(Bin_rand);
            Qrand_bin(1,:,irand) = M_bin.fi;
            
        end % random networks
        Net_para.bin.zfi(ithres,:,isub) = Net_para.bin.fi(ithres,:,isub)./mean(Qrand_bin,3);
        
        fprintf ('Threshold = %.3f ...... done\n', Thres(ithres));
        
    end % thresholds
    
end % subjects

% normalization
Net_para.bin.cpratio = Net_para.bin.cp ./ nanmean(Net_para.bin.cprand,3);
Net_para.bin.lpratio = Net_para.bin.lp ./ nanmean(Net_para.bin.lprand,3);

% area under curve
for isub = 1:size(Mats,1)
    
    Net_para.bin.acp(isub,1)      = gretna_auc(Net_para.bin.cp(:,isub),Delta);
    Net_para.bin.alp(isub,1)      = gretna_auc(Net_para.bin.lp(:,isub),Delta);
    Net_para.bin.acpratio(isub,1) = gretna_auc(Net_para.bin.cpratio(:,isub),Delta);
    Net_para.bin.alpratio(isub,1) = gretna_auc(Net_para.bin.lpratio(:,isub),Delta);
    
    Node_para.bin.adeg(isub,:)   = gretna_auc(squeeze(Node_para.bin.deg(:,isub,:)),Delta);

    ac = zeros(N_region);
    adis = zeros(N_region);
    for iregion = 1:N_region
        for jregion = iregion+1:N_region
            c_data = squeeze(Net_para.bin.c(iregion,jregion,:,isub));
            ac(iregion,jregion) = gretna_auc(c_data,Delta);     
            
            c_data = squeeze(Net_para.bin.dis(iregion,jregion,:,isub));
            adis(iregion,jregion) = gretna_auc(c_data,Delta);   
        end
    end
    ac = ac + ac';
    Net_para.bin.ac(:,:,isub) = ac;
    
    adis = adis + adis';
    Net_para.bin.adis(:,:,isub) = adis;
    
end

return