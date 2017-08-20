function postSegment(img_Loc)
%
% This function performs automated clean up on the output of Segmentation.
% It does the following:
%
% 1. Smooths the mask if it's not clean enough, to avoid convergence problem
% in the subsequent meshing using ScanIP.
% 2. Creates binary masks for each tissue class.
% 3. Fixes csf continuity.
% 4. Identify the disconnected voxels, and make them unassigned voxels
% (empty voxels).
% 5. Relabels each empty voxel to its nearest tissue type.
%
% Adapted from Yu (Andy) Huang, 2013.05

cd(img_Loc);


gray = nii2img('c1nii_for_seg.nii');
white = nii2img('c2nii_for_seg.nii');
csf = nii2img('c3nii_for_seg.nii');
bone = nii2img('c4nii_for_seg.nii');
skin = nii2img('c5nii_for_seg.nii');
cere = nii2img('c6nii_for_seg.nii');
orbits = nii2img('c7nii_for_seg.nii');
air = nii2img('c8nii_for_seg.nii');

white = load_untouch_nii(white);
gray = load_untouch_nii(gray);
csf = load_untouch_nii(csf);
bone = load_untouch_nii(bone);
skin = load_untouch_nii(skin);
cere = load_untouch_nii(cere);
orbits = load_untouch_nii(orbits);
air = load_untouch_nii(air);

gray_temp = gray.img; white_temp = white.img; csf_temp = csf.img;
cere_temp = cere.img; bone_temp = bone.img; skin_temp = skin.img; 
orbits_temp = orbits.img; 
air_temp = air.img;

% Create binary masks for each tissue class
% at each position, what has the greatest probability
[~,gray_temp,white_temp,csf_temp,bone_temp,...
    skin_temp,cere_temp,air_temp,orbits_temp] = binaryMaskGenerate(...
    gray_temp,white_temp,csf_temp,bone_temp,...
    skin_temp,cere_temp,air_temp,orbits_temp);

% save the results with the same header info as the input
gray.img = uint8(gray_temp)*255; save_untouch_nii(gray,'mask_gray.img');
white.img = uint8(white_temp)*255; save_untouch_nii(white,'mask_white.img');
csf.img = uint8(csf_temp)*255; save_untouch_nii(csf,'mask_csf.img');
bone.img = uint8(bone_temp)*255; save_untouch_nii(bone,'mask_bone.img');
skin.img = uint8(skin_temp)*255; save_untouch_nii(skin,'mask_scalp.img');
cere.img = uint8(cere_temp)*255; save_untouch_nii(cere,'mask_cere.img');
orbits.img = uint8(orbits_temp)*255; save_untouch_nii(orbits,'mask_orbits.img');
air.img = uint8(air_temp)*255; save_untouch_nii(air,'mask_air.img');

end

function [empt,mask1,mask2,mask3,mask4,mask5,mask6,mask7,mask8,mask9]...
    = binaryMaskGenerate(data1,data2,data3,data4,data5,data6,data7,data8,data9)

% Use max operation to generate binary mask for each tissue type from the
% results generated from SPM8. It can also be used to update each mask
% after one mask is processed.
% Output: empty mask represents those empty voxels which do not
% belong to any tissue type
%
% (c) Yu Huang (Andy), May 2011
% The Neural Engineering Lab, Dept. of Biomedical Engineering, 
% City College of New York

if nargin < 2
    data2 = [];
end
if nargin < 3
    data3 = [];
end
if nargin < 4
    data4 = [];
end
if nargin < 5
    data5 = [];
end
if nargin < 6
    data6 = [];
end
if nargin < 7
    data7 = [];
end
if nargin < 8
    data8 = [];
end
if nargin < 9
    data9 = [];
end

data = [];
data(:,:,:,1) = zeros(size(data1));
data(:,:,:,2) = data1;
if ~isempty(data2)
    data(:,:,:,3) = data2;
end
if ~isempty(data3)
    data(:,:,:,4) = data3;
end
if ~isempty(data4)
    data(:,:,:,5) = data4;
end
if ~isempty(data5)
    data(:,:,:,6) = data5;
end
if ~isempty(data6)
    data(:,:,:,7) = data6;
end
if ~isempty(data7)
    data(:,:,:,8) = data7;
end
if ~isempty(data8)
    data(:,:,:,9) = data8;
end
if ~isempty(data9)
    data(:,:,:,10) = data9;
end


[~,maxind] = max(data,[],4);

empt = (maxind==1);
mask1 = (maxind==2);
if ~isempty(data2)
    mask2 = (maxind==3);
end
if ~isempty(data3)
    mask3 = (maxind==4);
end
if ~isempty(data4)
    mask4 = (maxind==5);
end
if ~isempty(data5)
    mask5 = (maxind==6);
end
if ~isempty(data6)
    mask6 = (maxind==7);
end
if ~isempty(data7)
    mask7 = (maxind==8);
end
if ~isempty(data8)
    mask8 = (maxind==9);
end
if ~isempty(data9)
    mask9 = (maxind==10);
end


end