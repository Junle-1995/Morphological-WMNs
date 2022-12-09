function [Vout, Ydata, XYZ] = gretna_read_image(Path_filename)

%==========================================================================
% This function is used to read volume data of single nifti image.
%
%
% Syntax: function [Vout Ydata XYZ] = gretna_read_vol(Path_filename)
%
% Input:
%          Path_filename:
%                The directory & filename of the nifti image.
%
% Outputs:
%          Vout:
%                The header information of the nifti image.
%          Ydata:
%                The 3D data sorted in the nifti image.
%          XYZ:
%                The 3D coordination.
%
% Jinhui WANG, NKLCNL, BNU, BeiJing, 2011/10/23, Jinhui.Wang.1982@gmail.com
%==========================================================================

Vout = spm_vol(Path_filename);
[Ydata, XYZ] = spm_read_vols(Vout);

return