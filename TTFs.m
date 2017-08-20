
function TTFs

addpath(genpath('E:\spm8&TTFs'));

img_Path = spm_select;
[img_Loc, ~, ~] = fileparts(img_Path(1,:));
cd(img_Loc);

if length(img_Path) > 1
    %Converts the DICOM images to img/hdr format
    img_Path = DICOMConvert(img_Path);
end

%Coregister the padded image to a template
img_Info = spmCoregistration(img_Path);

%Segments the image
spmSegment('nii_for_seg.img');

%Applies a modified binary mask generation workflow 
postSegment(img_Loc);

%Generates and places electrodes on the surface of the scalp
TTF = fullfile(spm('Dir'),'tpm','TTF.nii');
electrodeGeneration('mask_scalp.img', TTF, 'iy_nii_for_seg.nii');

%Generates a python script for ScanIP import and smoothing
pyScanIPScript(img_Info, img_Loc);
 
%Optional step for bat script generation/automated ScanIP import
ScanIP_Loc = 'E:\Program Files\Simpleware\ScanIP';
batScript(img_Loc, ScanIP_Loc);

end