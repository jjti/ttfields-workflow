function img_Info = spmCoregistration(nii_Img)

ref_Nii = fullfile(spm('Dir'), 'templates', 'EPI.nii');

spm_jobman('initcfg');

matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {ref_Nii};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {nii_Img};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 1;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

spm_jobman('run',matlabbatch);

% This uses the reslice tool from the Tools for Nifti and ANALYZE Library 
evalc('reslice_nii(nii_Img, nii_Img)');
rNii = load_nii(nii_Img);

%Collect HDR information for Python script
x_spc = rNii.hdr.dime.pixdim(2); 
y_spc = rNii.hdr.dime.pixdim(3);
z_spc = rNii.hdr.dime.pixdim(4);
x_voxs = rNii.hdr.dime.dim(2);
y_voxs = rNii.hdr.dime.dim(3);
z_voxs = rNii.hdr.dime.dim(4);
img_Info = [x_spc, y_spc, z_spc, x_voxs, y_voxs, z_voxs];

%Saves the padded image and returns the location of the image to segment
save_nii(rNii, 'nii_for_seg.img');
end