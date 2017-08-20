function mySegment(fnm)
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

[dir, bnm, ~] = spm_fileparts(fnm);

disp('loading data...')

white = load_untouch_nii(['c1' bnm '.nii']);
gray = load_untouch_nii(['c2' bnm '.nii']);
csf = load_untouch_nii(['c3' bnm '.nii']);
bone = load_untouch_nii(['c4' bnm '.nii']);
skin = load_untouch_nii(['c5' bnm '.nii']);
cere = load_untouch_nii(['c6' bnm '.nii']);
orbits = load_untouch_nii(['c7' bnm '.nii']);
air = load_untouch_nii(['c8' bnm '.nii']);
% load the masks

gray_temp = gray.img; white_temp = white.img; csf_temp = csf.img;
cere_temp = cere.img; bone_temp = bone.img; skin_temp = skin.img; 
orbits_temp = orbits.img; air_temp = air.img;


% Smooth the mask if it's not clean enough, to avoid convergence problem
% in the subsequent meshing using ScanIP
disp('smoothing images...')
smt_fil = fspecial('gaussian', 3, 0.2);
for i = 1:size(gray_temp,3)
    gray_temp(:,:,i) = imfilter(gray_temp(:,:,i),smt_fil);
end

smt_fil = fspecial('gaussian', 3, 0.1);
for i = 1:size(white_temp,3)
    white_temp(:,:,i) = imfilter(white_temp(:,:,i),smt_fil);
end

smt_fil = fspecial('gaussian', 5, 0.3);
for i = 1:size(csf_temp,3)
    csf_temp(:,:,i) = imfilter(csf_temp(:,:,i),smt_fil);
end

smt_fil = fspecial('gaussian', 3, 0.4);
for i = 1:size(bone_temp,3)
    bone_temp(:,:,i) = imfilter(bone_temp(:,:,i),smt_fil);
end

smt_fil = fspecial('gaussian', 5, 1);
for i = 1:size(skin_temp,3)
    skin_temp(:,:,i) = imfilter(skin_temp(:,:,i),smt_fil);
end

smt_fil = fspecial('gaussian', 5, 1);
for i = 1:size(air_temp,3)
    air_temp(:,:,i) = imfilter(air_temp(:,:,i),smt_fil);
end

smt_fil = fspecial('gaussian', 5, 0.2);
for i = 1:size(cere_temp,3)
    cere_temp(:,:,i) = imfilter(cere_temp(:,:,i),smt_fil);
end

smt_fil = fspecial('gaussian', 5, 0.4);
for i = 1:size(orbits_temp,3)
    orbits_temp(:,:,i) = imfilter(orbits_temp(:,:,i),smt_fil);
end

% Create binary masks for each tissue class
[empt_temp,cere_temp,gray_temp,white_temp,orbits_temp,...
    csf_temp,bone_temp,skin_temp,air_temp]...
    = binaryMaskGenerate(cere_temp,gray_temp,...
    white_temp,orbits_temp,csf_temp,bone_temp,skin_temp,air_temp);

% Fix CSF continuity
se=ones(4,4,4);
dcsf=imdilate(csf_temp, se);
dbone=imdilate(bone_temp, se);
contin=(empt_temp&dcsf)|(dbone&gray_temp);
csf_temp=csf_temp|contin;
[~,csf_temp,bone_temp,gray_temp]...
    = binaryMaskGenerate(csf_temp,bone_temp,gray_temp);

disp('removing disconnected voxels...')

siz = sizeOfObject(cere_temp);
thres = siz(1)-1;
cere_temp=bwareaopen(cere_temp,thres);

siz = sizeOfObject(gray_temp);
thres = siz(1)-1;
gray_temp = bwareaopen(gray_temp,thres);

siz = sizeOfObject(white_temp);
thres = siz(1)-1;
white_temp = bwareaopen(white_temp,thres);

siz = sizeOfObject(csf_temp);
thres = siz(4)+1;
csf_temp = bwareaopen(csf_temp,thres);

siz = sizeOfObject(bone_temp);
thres = siz(1)-1;
bone_temp = bwareaopen(bone_temp,thres);

siz = sizeOfObject(skin_temp);
thres = siz(1)-1;
skin_temp = bwareaopen(skin_temp,thres);

siz = sizeOfObject(air_temp);
thres = siz(5)+1;
air_temp=bwareaopen(air_temp,thres);

% Identify the disconnected voxels
disp('generating and labeling empty voxels...')
empt_temp = binaryMaskGenerate(cere_temp,gray_temp,white_temp,...
    orbits_temp,bone_temp,csf_temp,skin_temp,air_temp);

% Generate unassigned voxels (empty voxels)
for j = 1:5 % usually all empty voxels will be labelled ina  few loops
    gray_fil = uint8(gray_temp)*255; white_fil = uint8(white_temp)*255;
    csf_fil = uint8(csf_temp)*255; bone_fil = uint8(bone_temp)*255;
    skin_fil = uint8(skin_temp)*255; air_fil = uint8(air_temp)*255;
    cere_fil = uint8(cere_temp)*255; orbits_fil = uint8(orbits_temp)*255;
    smt_fil = fspecial('gaussian', 5, 1.3);
    for i = 1:size(csf_fil,3)
        img = csf_fil(:,:,i);
        csf_fil(:,:,i) = imfilter(img,smt_fil);
    end
    smt_fil = fspecial('gaussian', 5, 1);
    for i = 1:size(orbits_fil,3)
        img = orbits_fil(:,:,i);
        orbits_fil(:,:,i) = imfilter(img,smt_fil);
    end
    smt_fil = fspecial('gaussian', 5, 1);
    for i = 1:size(cere_fil,3)
        img = cere_fil(:,:,i);
        cere_fil(:,:,i) = imfilter(img,smt_fil);
    end
    smt_fil = fspecial('gaussian', 5, 1);
    for i = 1:size(gray_fil,3)
        img = gray_fil(:,:,i);
        gray_fil(:,:,i) = imfilter(img,smt_fil);
    end
    smt_fil = fspecial('gaussian', 5, 1);
    for i = 1:size(white_fil,3)
        img = white_fil(:,:,i);
        white_fil(:,:,i) = imfilter(img,smt_fil);
    end
    smt_fil = fspecial('gaussian', 5, 1);
    for i = 1:size(bone_fil,3)
        img = bone_fil(:,:,i);
        bone_fil(:,:,i) = imfilter(img,smt_fil);
    end
    smt_fil = fspecial('gaussian', 5, 1);
    for i = 1:size(skin_fil,3)
        img = skin_fil(:,:,i);
        skin_fil(:,:,i) = imfilter(img,smt_fil);
    end
    smt_fil = fspecial('gaussian', 5, 1);
    for i = 1:size(air_fil,3)
        img = air_fil(:,:,i);
        air_fil(:,:,i) = imfilter(img,smt_fil);
    end
    
    [~,cere_fil,gray_fil,white_fil,orbits_fil,csf_fil,bone_fil,skin_fil,air_fil]...
        = binaryMaskGenerate(cere_fil,gray_fil,white_fil,orbits_fil,csf_fil,bone_fil,skin_fil,air_fil);

    gray_temp = (empt_temp&gray_fil)|gray_temp;
    white_temp = (empt_temp&white_fil)|white_temp;
    orbits_temp = (empt_temp&orbits_fil)|orbits_temp;
    csf_temp = (empt_temp&csf_fil)|csf_temp;
    skin_temp = (empt_temp&skin_fil)|skin_temp;
    cere_temp = (empt_temp&cere_fil)|cere_temp;
    bone_temp = (empt_temp&bone_fil)|bone_temp;
    air_temp = (empt_temp&air_fil)|air_temp;
    
    empt_temp = xor(empt_temp,((empt_temp&gray_fil) | (empt_temp&skin_fil)...
        | (empt_temp&white_fil) | (empt_temp&csf_fil) | (empt_temp&white_fil)...
        | (empt_temp&orbits_fil) | (empt_temp&bone_fil) |...
        (empt_temp&cere_fil)) | (empt_temp&air_fil)); % update empty voxels


% Relabel each empty voxel to its nearest tissue type
% The Gaussian filter is used to calculate distances, and max operation
% relabels each empty voxel based on the distances.
end
cd(dir);
disp('saving masks...')
gray.img = uint8(gray_temp)*255; save_untouch_nii(gray,'mask_gray.nii');
white.img = uint8(white_temp)*255; save_untouch_nii(white,'mask_white.nii');
orbits.img = uint8(orbits_temp)*255; save_untouch_nii(orbits,'mask_orbits.nii');
csf.img = uint8(csf_temp)*255; save_untouch_nii(csf,'mask_csf.nii');
bone.img = uint8(bone_temp)*255; save_untouch_nii(bone,'mask_bone.nii');
skin.img = uint8(skin_temp)*255; save_untouch_nii(skin,'mask_skin.nii');
air.img = uint8(air_temp)*255; save_untouch_nii(air,'mask_air.nii');
cere.img = uint8(cere_temp)*255; save_untouch_nii(cere,'mask_cere.nii');
% save the results with the same header info as the input

dcm2nii = fullfile(spm('Dir'),'Other/mricron');
cd(dcm2nii);

disp('converting to imgs...')
nii2img([dir '\mask_gray.nii']);
nii2img([dir '\mask_white.nii']);
nii2img([dir '\mask_orbits.nii']);
nii2img([dir '\mask_csf.nii']);
nii2img([dir '\mask_bone.nii']);
nii2img([dir '\mask_skin.nii']);
nii2img([dir '\mask_air.nii']);
nii2img([dir '\mask_cere.nii']);

disp('Done');
end

function nii2img(fnm)
first = ['dcm2nii -n n -m n ' fnm];
system(first);
end

function [empt,mask1,mask2,mask3,mask4,mask5,mask6,mask7,mask8,mask9]...
    = binaryMaskGenerate(data1,data2,data3,data4,data5,data6,data7,data8,data9)

% [empt,mask1,mask2,mask3,mask4,mask5,mask6] = binaryMaskGenerate(data1,data2,data3,data4,data5,data6)
%
% Use max operation to generate binary mask for each tissue type from the
% results generated from SPM8. It can also be used to update each mask
% after one mask is processed.
% Input: data1~data6: 6 tissue types, at least need 1 tissue type;
% Output: empt: empty mask represents those empty voxels which do not
% belong to any tissue type; mask1~mask6: 6 masks corresponding to
% data1~data6.
%
% (c) Yu Huang (Andy), May 2011
% The Neural Engineering Lab, Dept. of Biomedical Engineering, City College of New York
% Send bugs to yhuang16@citymail.cuny.edu

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