function TTFs

img_Path = spm_select;
[img_Loc, ~, ~] = fileparts(img_Path(1,:));
cd(img_Loc);

if length(img_Path) > 1
    % Converts the DICOM images to img/hdr format
    img_Path = DICOMConvert(img_Path);
end

% Coregister the padded image to a template
img_Info = spmCoregistration(img_Path);

% Segments the image using nii_for_seg.img in last step
tpmsPath = [pwd filesep 'tpms'];
spmSegment(tpmsPath, 'nii_for_seg.img');

% Applies a modified binary mask generation workflow 
postSegment(img_Loc);

% Generates and places electrodes on the surface of the scalp
elecSeeds = [pwd filesep 'tpms' filesep 'TTF.nii']; % TTField transducer seed coordinates
electrodeGeneration('mask_scalp.img', elecSeeds, 'iy_nii_for_seg.nii');

% Generates a python script for ScanIP import and smoothing
% if uncommented, this would generate a python script that will
% automated mask import and postprocessing in ScanIP
%
% pyScanIPScript(img_Info, img_Loc);


% Optional step for bat script generation/automated ScanIP import
% useful for Windows if using ScanIP. Generates a batch script for
% and calls it to initiate ScanIP and the import script in the last
% step
%
% ScanIP_Loc = 'C:\Program Files\Simpleware\ScanIP';
% batScript(img_Loc, ScanIP_Loc);

end