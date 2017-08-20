function [padded_Img, img_Info] = padNifti(unpadded_Nifti)
% This function is used in order to reslice the image, getting rid of 
% potential variable voxel dimensions issues and also pads the edges of
% the image in order to ensure that there's space for the generation and 
% placement of electrodes.

% Also collects file information for ScanIP script generation

org_File = unpadded_Nifti;

% This uses the reslice tool from the Tools for Nifti and ANALYZE Library 
pad_spec = struct('pad_from_L', 20, 'pad_from_R', 20, 'pad_from_P', 20, ...
    'pad_from_A', 20, 'pad_from_I', 0, 'pad_from_S', 20, 'bg', 0);
evalc('reslice_nii(unpadded_Nifti, unpadded_Nifti)');
unpadded_Nifti = load_nii(unpadded_Nifti);
padded_Nifti = pad_nii(unpadded_Nifti, pad_spec);

%Collect HDR information for Python script
x_spc = padded_Nifti.hdr.dime.pixdim(2); 
y_spc = padded_Nifti.hdr.dime.pixdim(3);
z_spc = padded_Nifti.hdr.dime.pixdim(4);
x_voxs = padded_Nifti.hdr.dime.dim(2);
y_voxs = padded_Nifti.hdr.dime.dim(3);
z_voxs = padded_Nifti.hdr.dime.dim(4);
img_Info = [x_spc, y_spc, z_spc, x_voxs, y_voxs, z_voxs];

%Saves the padded image and returns the location of the image to segment
save_nii(padded_Nifti, 'nii_for_seg.img');
[nii_Path, ~, ~] = fileparts(org_File);
padded_Img = [nii_Path filesep 'nii_for_seg.img'];

end