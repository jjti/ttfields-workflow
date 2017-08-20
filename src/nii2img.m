function img = nii2img(nii)
%Used to convert all the nifti (nii) segmentations from spmSegment to img
%format. ScanIP is unable to properly import nii files.

[~,fnm,~]=fileparts(nii);

V=spm_vol(nii);
outgoing=spm_read_vols(V);
V.fname=[fnm,'.img'];
spm_write_vol(V,outgoing);

img = V.fname;
delete(nii);

end