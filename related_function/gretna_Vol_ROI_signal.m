function gretna_Vol_ROI_signal(Data_path, File_filter, Template_path, Template_filter, RoiIndex, Output_path, Level, Num_thres)

%==========================================================================
% This function is used to extract signals for each voxel within each ROI.
%
%
% Syntax: function gretna_Vol_ROI_signal(Data_path, File_filter, Template_path, Template_filter, RoiIndex, Output_path, Level, Num_thres)
%
% Inputs:
%         Data_path:
%                   The directory & filename of a .txt file that contains
%                   the directory of those files to be processed (can be
%                   obtained by gretna_gen_data_path.m).
%       File_filter:
%                   The prefix of those files to be processed.
%     Template_path:
%                   The directory & filename of a .txt file that contains
%                   the directory of individual templates (can be obtained 
%                   by gretna_gen_data_path.m).
%   Template_filter:
%                   The prefix of those templates.
%          RoiIndex:
%                   The index of ROIs in the brain template or atlas (e.g.,
%                   [1:40 60:90]). Note that the order of extracted time
%                   course is the same as the order entered in RoiIndex.
%       Output_path:
%                   The directory where the resulting MTC are sorted.
%             Level:
%                   'all':  extract all signals in each ROI;
%                   'mean': extract mean signal in each ROI;
%                   'both': extract both all and mean signals in each ROI.
%         Num_thres:
%                   Positive integer threshold for regional number of
%                   voxels. Regions with voxels less than the threshold are
%                   excluded (Default = 128).
%
% Output:
%    Signal_xxx.mat:
%                   N (# of images) * M (# of ROIs) cell array for each
%                   subject with each cell containing a data array. NB. the
%                   outputted files are sorted in the same order as that in
%                   Data_path.
% mean_Signal_xxx.mat:
%                   N (# of images) * M (# of ROIs) data array for each
%                   subject. NB. the outputted files are sorted in the same
%                   order as that in Data_path.
%
% Jinhui WANG,  NKLCNL, BNU,  BeiJing,   2011/10/23, Jinhui.Wang.1982@gmail.com
% Hao WANG,     CCBD,   HZNU, Hangzhou,  2014/10/28, hall.wong@outlook.com
% Ningkai WANG, IBRR,   SCNU, Guangzhou, 2020/04/02, Ningkai.Wang.1993@gmail.com
% Yuping YANG,  IBRR,   SCNU, Guangzhou, 2019/10/15, Yupingyanghvd@gmail.com
%==========================================================================

if  nargin == 7
    Num_thres = 128;
end

fid = fopen(Data_path);
Dir_data = textscan(fid, '%s');
fclose(fid);
Num_subs = size(Dir_data{1}, 1);

fid = fopen(Template_path);
Dir_template = textscan(fid, '%s');
fclose(fid);
Num_templates = size(Dir_template{1}, 1);

if Num_templates > 1
    if Num_templates ~= Num_subs
        error('The number of templates is not equal to the number of subjects!');
    end
end

Num_regs = length(RoiIndex);
Ind_regs = cell(Num_templates, Num_regs);
Num_voxs = zeros(Num_templates, Num_regs);

% Record regional index for each template
for i_tem = 1:Num_templates
    
    cd([Dir_template{1}{i_tem}])
    Template_name = spm_select('ExtList',pwd, ['^' Template_filter '.*\.img$'],inf);
    if isempty(Template_name)
        Template_name = spm_select('ExtList',pwd, ['^' Template_filter '.*\.nii$'],inf);
    end
    
    [Vtem, Ytem, ~] = gretna_read_image(Template_name);
    Ytem(isnan(Ytem)) = 0;
    
    for i_reg = 1:Num_regs
        Region = RoiIndex(i_reg);
        Ind_regs{i_tem, i_reg} = find(Region == Ytem(:));
        Num_voxs(i_tem, i_reg) = length(Ind_regs{i_tem, i_reg});
        if isempty(Ind_regs)
            error (['There are no voxels in ROI' blanks(1) num2str(RoiIndex(i_reg)) ', please specify ROIs again']);
        end
    end
end
clear Ytem

if Num_templates == 1
    Ind_regs = repmat(Ind_regs, Num_subs, 1);
end

% Extracting regional signal for each subject
Num_vox_reg       = zeros(Num_subs, Num_regs);
removed_subs_reg  = zeros(Num_subs, Num_regs);
Name_sub          = cell(Num_subs, 1);
Cell_removed_regs = cell(Num_subs, 1);

for i_sub = 1:Num_subs
    
    fprintf('Extracting regional signal for %s\n', Dir_data{1}{i_sub});
    
    cd ([Dir_data{1}{i_sub}])
    File_name = spm_select('ExtList', pwd, ['^' File_filter '.*\.img$'],inf);
    if isempty(File_name)
        File_name = spm_select('ExtList', pwd, ['^' File_filter '.*\.nii$'],inf);
    end
    
    Vin = spm_vol(File_name);
    if  ~isequal(Vtem.dim, Vin(1).dim)
        error('The dimensions must be the same between the brain template (atlas) and the data!')
    end
    
    Num_imgs      = size(Vin, 1);
    Sig_roi       = cell(Num_imgs, Num_regs);
    Sig_roi_check = cell(1, Num_regs);
    mean_Sig_roi  = zeros(Num_imgs, Num_regs);
    
    for i_img = 1:Num_imgs
        [Ydata, ~] = spm_read_vols(Vin(i_img));
        for i_reg = 1:Num_regs          
            reg_Ydata = Ydata(Ind_regs{i_sub,i_reg});
            Sig_roi{i_img, i_reg}       = reg_Ydata(~isnan(reg_Ydata));
            Sig_roi_check{1, i_reg}     = Sig_roi{1, i_reg};
            mean_Sig_roi(i_img, i_reg) = nanmean(Ydata(Ind_regs{i_sub, i_reg}));
        end
    end
    
    Num_vox_reg(i_sub,:)         = cellfun('size', Sig_roi_check, 1);
    removed_reg                  = Num_vox_reg(i_sub,:) < Num_thres;
    removed_subs_reg(i_sub, :)   = removed_reg;
    Sig_roi(:, removed_reg)      = [];
    mean_Sig_roi(:, removed_reg) = [];
    
    Name_sub{i_sub, 1} = ['sub_' num2str(i_sub,'%04d') '.mat'];
    
    switch lower(Level)
        case 'all'
            save([Output_path filesep 'Signal_'      Name_sub{i_sub, 1}],'Sig_roi');
        case 'mean'
            save([Output_path filesep 'mean_Signal_' Name_sub{i_sub, 1}],'mean_Sig_roi');
        case 'both'
            save([Output_path filesep 'Signal_'      Name_sub{i_sub, 1}],'Sig_roi');
            save([Output_path filesep 'mean_Signal_' Name_sub{i_sub, 1}],'mean_Sig_roi');
    end
    
    fprintf('Extracting regional signal for %s ...... is done\n', [Dir_data{1}{i_sub}]);
    
end

% Generate an excel file containing important atlas informations
if  Num_templates == 1
    Altas_name = '';
    Supported_atlas = gretna_label('vol');
    Supported_atlas = Supported_atlas.abbr;
    
    for i_vol = 1:length(Supported_atlas)
        if ~isempty(strfind(Template_name, Supported_atlas{i_vol}));
            Altas_name  = Supported_atlas{i_vol};
            if strcmpi(Altas_name, 'craddock_200')
                Region_list    = num2cell(RoiIndex(:));
            elseif strcmpi(Altas_name, 'power_264')
                Region_list    = gretna_label(Altas_name);
                Region_list    = Region_list.functional_classification(RoiIndex,1);
            else
                Region_list    = gretna_label(Altas_name);
                Region_list    = Region_list.abbr(RoiIndex,1);
            end
        end
    end
    
    if ~exist('Region_list','var')
        Region_list = num2cell(RoiIndex(:));
    end
    
    Output_switch        = mean(removed_subs_reg,1);
    Num_mean_vox_regmean = mean(Num_vox_reg,1);
    Cell_full_regs = cat(2, num2cell(RoiIndex'), Region_list, num2cell(Num_mean_vox_regmean'));
    
    if mean(bsxfun(@eq, Output_switch, abs(fix(Output_switch)))) == 1
        Num_rem_regs          = sum(~removed_reg);
        Cell_ind_rem          = num2cell(1:Num_rem_regs);
        Cell_ori_ind_rem      = Cell_full_regs(~removed_reg,1);
        Cell_reserved         = Cell_full_regs(~removed_reg,2:3);
        Cell_reserved         = cat(2, Cell_ind_rem', Cell_ori_ind_rem, Cell_reserved);
        IND_ori_ind_removed   = 1:sum(removed_reg);
        Cell_ori_ind_removed  = num2cell(IND_ori_ind_removed');
        Cell_removed          = Cell_full_regs(removed_reg,:);
        
        Output_label_csv                              = cell(Num_regs+3, 13);
        Output_label_csv(4:end, 1:2)                  = Cell_full_regs(:,[1, 2]);
        Output_label_csv(4:end, 3)                    = num2cell(Num_voxs);
        Output_label_csv(4:3+Num_rem_regs, 5:8)       = Cell_reserved;
        Output_label_csv(4:3+sum(removed_reg), 10)    = Cell_ori_ind_removed;
        Output_label_csv(4:3+sum(removed_reg), 11:13) = Cell_removed;
        
        Output_label_csv(3, [1, 5, 10]) = {'Index'};
        Output_label_csv(3, [2, 7, 12]) = {'Label'};
        Output_label_csv(3, [3, 8, 13]) = {'Voxels'};
        Output_label_csv(3, [6, 11])    = {'Original_Index'};
        Output_label_csv(2, [1, 5, 10]) = {'Original', 'Reserved', 'Removed'};
    else
        Output_label_csv = cell(Num_regs+3, 3);
        Output_label_csv(3, 1) = {'Index'};
        Output_label_csv(3, 2) = {'Label'};
        Output_label_csv(3, 3) = {'Voxels'};
        Output_label_csv(2, 1) = {'Original'};
    end
    
    Output_label_csv(1,1) = {['The threshold you set is: ' num2str(Num_thres)]};
    
    for i_sub = 1:Num_subs
        Cell_removed_regs{i_sub,1} = find(removed_subs_reg(i_sub,:) == 1);
    end
    
    Output_subs_csv = table(Dir_data{1} ,Name_sub, sum(~removed_subs_reg,2), sum(removed_subs_reg,2),...
        Cell_removed_regs, 'VariableNames', {'Location_file', 'ID_sub', 'Num_reserved_regs', 'Num_removed_regs', 'List_removed_regs'});
    
    writetable(Output_subs_csv, [Output_path filesep sprintf('Table_Subjects_Region_Modified_%s.csv', Altas_name)], 'WriteVariableNames', 1);
    fprintf('\nTable_Subjects_Region_Modified_%s.csv has been output to %s successfully, please check it!\n', Altas_name, Output_path);
    
    writetable(cell2table(Output_label_csv), [Output_path filesep sprintf('Table_Labels_%s.csv',Altas_name)], 'WriteVariableNames', 0);
    fprintf('\nTable_Labels_%s.csv has been output to %s successfully, please check it!\n', Altas_name, Output_path);
    
    cd(Output_path);
end

return