%% construct the JS-based individual morphological WM networks based on the preprocessed structural images
%% required toolbox: SPM12
path = pwd;

% extract the voxel-vise signal within the roi 
gretna_Vol_ROI_signal('Z:\Data\HCP1200_T1\path_anat_mri.txt', 'mwp2',...
    [path '\Template.txt'], 'JHU-WhiteMatter-labels_1.5mm', 1:48, [path '\VBM_nosmooth'],'both',100)
gretna_Vol_ROI_signal('Z:\Data\HCP1200_T1\path_anat_mri.txt', 'wj',...
    [path '\Template.txt'], 'JHU-WhiteMatter-labels_1.5mm', 1:48, [path '\DBM_nosmooth'],'both',100)
gretna_Vol_ROI_signal('Z:\Data\HCP1200_T1\path_anat_mri.txt', 'smo888_mwp2',...
    [path '\Template.txt'], 'JHU-WhiteMatter-labels_1.5mm', 1:48, [path '\VBM_smooth'],'both',100)
gretna_Vol_ROI_signal('Z:\Data\HCP1200_T1\path_anat_mri.txt', 'smo888_wj',...
    [path '\Template.txt'], 'JHU-WhiteMatter-labels_1.5mm', 1:48, [path '\DBM_smooth'],'both',100)

% calculate the PDF of the voxel-vise signal within the roi
gretna_Sim_PDF_matrix([path '\VBM_nosmooth'], 'Signal')
gretna_Sim_PDF_matrix([path '\DBM_nosmooth'], 'Signal')
gretna_Sim_PDF_matrix([path '\VBM_smooth'], 'Signal')
gretna_Sim_PDF_matrix([path '\DBM_smooth'], 'Signal')

