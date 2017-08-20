function outgoing_img = DICOMConvert(img_Path)
% Uses SPM's DICOMImport functions to convert a DICOM folder of images
% into a single img (img format rather than nifti because ScanIP won't 
% import nifti format segmentations.
% 
% We have had issues with there being two niftis per patient file
% (one being a guide for the radiologist)
% so the smaller of the two img files is removed post-conversion

header_Array = spm_dicom_headers(img_Path); 
nifti_Images = spm_dicom_convert(header_Array, 'all', 'flat', 'nii');

% If two files, deletes the smaller of the two
if (length(nifti_Images.files) == 2)
    File1 = dir(char(nifti_Images.files{1}));
    File2 = dir(char(nifti_Images.files{2}));
    if File1.bytes > File2.bytes
        delete(nifti_Images.files{2});
        delete([nifti_Images.files{2}, '.hdr']);
        outgoing_img = char(nifti_Images.files{1});
    elseif File2.bytes > File1.bytes
        delete(nifti_Images.files{1});
        outgoing_img = char(nifti_Images.files{2});
    end
else
    outgoing_img = char(nifti_Images.files{1});
end

end