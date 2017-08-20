function elec_structs = electrodeGeneration(scalp_mask, arrayTemplate, inverseField)
% scalpRead(scalp_mask, arrayTemplate, inverseField)
% Joshua Timmons, 3/30/2016

nii = load_untouch_nii(scalp_mask);
loc = pwd;

%Applies inverse deformation field to elecorode template
spmDeformation(arrayTemplate, inverseField); 

%Gathers starting points for seed voxels based on the deformation field
array_info = seedFinder(loc);

%Makes filled sphere of points
sphere_coords = sphereMaker(nii);

%Determines the electrodes starting points against surface of scalp
seed_voxels = seedAdjuster(array_info, nii);

%Normal vector determination for each seed point
elec_structs = electrodeDesigner(nii, seed_voxels, sphere_coords, array_info);

%Creates the electrodes and gel masks
electrodeMaker(nii, elec_structs);

end

function spmDeformation(arrayTemplate, deformField)
%Calls SPM's built-in deformation utility which warps TTF template to brain
%using the inverse deformation field

spm_jobman('initcfg');

matlabbatch{1}.spm.util.defs.comp{1}.def = {deformField};
matlabbatch{1}.spm.util.defs.ofname = '';
matlabbatch{1}.spm.util.defs.fnames = {arrayTemplate};
matlabbatch{1}.spm.util.defs.savedir.savepwd = 1;
matlabbatch{1}.spm.util.defs.interp = 1;
    
spm_jobman('run',matlabbatch);

end

function array_info = seedFinder(location)

cd(location);
init_template = load_untouch_nii('wTTF.nii');
init_template = init_template.img;

array_list = zeros(36,3); 
center_list = zeros(6,3);

%Requires Image Processing Toolbox. Assigning array location based on
%maximum gray scale intensity of each connected region
Connected_regions = bwconncomp(init_template); 
Region_centers = regionprops(Connected_regions, 'centroid'); 
Cen_ind = 0; Array_ind = 0; 

for i=2:size(Region_centers,1)
    %Find each centroid in the 3D image
    centroid = round([Region_centers(i).Centroid(2), ...
        Region_centers(i).Centroid(1), Region_centers(i).Centroid(3)]);
    %Convert to a linear index
    linear_indexes = Connected_regions.PixelIdxList{i};
    [x_index,y_index,z_index]= ...
        ind2sub(size(init_template), linear_indexes(:));
    img_coord_list = [x_index(:), y_index(:), z_index(:)];
    intensity_list = zeros(size(img_coord_list));
    %Create a list of pixels in this centroid
    for j=1:size(img_coord_list,1)
        intensity_list(j) = init_template(img_coord_list(j,1), ...
            img_coord_list(j,2), img_coord_list(j,3));
    end
    %Find the max intensity and assign this centroid that value
    intensity = max(intensity_list(:,1));
    if (intensity > 0.95)
        Array_ind = Array_ind + 1;
        array_list(Array_ind,:) = centroid;
    end
    if (intensity > 0.70) && (intensity < 0.81)
        Cen_ind = Cen_ind + 1;
        center_list(Cen_ind,:) = centroid; 
    end 
end

%Center of the brain determined by the center of the proximal points
center_point = round(mean(center_list));
array_info = struct('total_list', array_list, 'center', center_point);

end

function sphere_coords = sphereMaker(nii)

% Gets the voxel spacing 
spacing = nii.hdr.dime.pixdim(2);

%An approximate radius of the electrode. 
%A larger radius was found to work best on masks with uneven surfaces
sph_radius = 10/spacing; 
[sph_x, sph_y, sph_z] = sphere(49); 

%Fills the sphere
adj_sph_x = zeros(500,50); 
adj_sph_y = zeros(500,50); 
adj_sph_z = zeros(500,50);
for i = 1:10
    adj_sph_x((i*50-49):(i*50),:) = round(i*sph_x*0.1*sph_radius);
    adj_sph_y((i*50-49):(i*50),:) = round(i*sph_y*0.1*sph_radius);
    adj_sph_z((i*50-49):(i*50),:) = round(i*sph_z*0.1*sph_radius);
end
sphere_coords = struct('sph_x', adj_sph_x, 'sph_y', adj_sph_y, ...
    'sph_z', adj_sph_z);
end

function seed_voxels = seedAdjuster(array_info, nii)

%Gets the voxel spacing in 3 dimensions
spacing = nii.hdr.dime.pixdim(2); 

%Moves seeds to the edge of the scalp as starting points
elec_list = zeros(36,3);
for i = 1:36
    elec_list(i,:) = ...
        seedToEdge(nii, array_info.total_list(i,:), array_info);
end

%Reorders based on proximity and sets their distance as 4th value
adj_elec = zeros(36,4); row_ind = 1;
warning('off','all')
for i = 1:36
[index, dist] = knnsearch(elec_list(i,:), elec_list, 3);
dist = dist*spacing;
    if dist(2) < 30 && dist(3) < 30
        adj_elec(row_ind,:) = [elec_list(index(2),:), dist(2)-22];
        adj_elec(row_ind+1,:) = [elec_list(index(1),:), dist(1)];
        adj_elec(row_ind+2,:) = [elec_list(index(3),:), dist(3)-22];
        row_ind = row_ind + 3;
    end
end
warning('on','all')

%Checks whether pairwise electrode seed points are the proper distance
%apart (to within half a mm). Adjusts if not.
for i=2:3:35
   if adj_elec(i-1,4) < -0.5
       dist = abs(adj_elec(i-1,4));
       vect = (adj_elec(i-1,1:3)-adj_elec(i,1:3)); vect = vect/norm(vect);
       adj_elec(i-1,1:3) = adj_elec(i-1,1:3) + dist*vect;
   end
   if adj_elec(i-1,4) > 0.5
       dist = abs(adj_elec(i-1,4));
       vect = (adj_elec(i,1:3)-adj_elec(i-1,1:3)); vect = vect/norm(vect);
       adj_elec(i-1,1:3) = adj_elec(i-1,1:3) + dist*vect;
   end
   if adj_elec(i+1,4) < -0.5
       dist = abs(adj_elec(i+1,4));
       vect = (adj_elec(i+1,1:3)-adj_elec(i,1:3)); vect = vect/norm(vect);
       adj_elec(i+1,1:3) = adj_elec(i+1,1:3) + dist*vect;
   end
   if adj_elec(i+1,4) > 0.5
       dist = abs(adj_elec(i+1,4));
       vect = (adj_elec(i,1:3)-adj_elec(i+1,1:3)); vect = vect/norm(vect);
       adj_elec(i+1,1:3) = adj_elec(i+1,1:3) + dist*vect;
   end
end

adj_elec = round(adj_elec);
outgoing_seeds = zeros(36,3);
for i = 1:36
outgoing_seeds(i,1:3) = seedToEdge(nii, adj_elec(i,1:3), array_info);
end

%Recollects all the points
seed_voxels = outgoing_seeds;

end

function seed_voxel = seedToEdge(nii, seed_voxel, array_info)
% Creates a vector between the seed_voxel and the center of the image.
% Then checks whether the seed_voxel is on the inside or outside of the
% scalp and adjusts accordingly. 

img_center = array_info.center;

%Creates a vector to center. Corrects direction if needed
inward_vector = [img_center(1)-seed_voxel(1), ...
    img_center(2)-seed_voxel(2), img_center(3)-seed_voxel(3)];
inward_vector = inward_vector/(norm(inward_vector));
test_voxel = seed_voxel + 10*inward_vector;
distance_to_center = sqrt((img_center(1)-seed_voxel(1))^2 + ...
    (img_center(2)-seed_voxel(2))^2 + (img_center(3)-seed_voxel(3))^2);
test_distance = sqrt((img_center(1)-test_voxel(1))^2 + ...
    (img_center(2)-test_voxel(2))^2 + (img_center(3)-test_voxel(3))^2);

if test_distance > distance_to_center
    outward_vector = inward_vector;
    inward_vector = -inward_vector;
else
    outward_vector = -inward_vector;
end

%Makes line of points between seed and center to determine seed's position
dir_check = 0;
distance_to_center = round(distance_to_center);
for i=0:distance_to_center
    voxel_step = seed_voxel + i*inward_vector;
    voxel_step = round(voxel_step);
    output_val = nii.img(voxel_step(1), voxel_step(2), voxel_step(3));
    dir_check = dir_check + output_val;
end

%If the voxel is inside scalp, moves it outwards until it reaches the scalp
if dir_check == 0
init_voxel = seed_voxel; inc = 0;
while nii.img(seed_voxel(1), seed_voxel(2), seed_voxel(3)) == 0 
    inc_vector = outward_vector;
    inc = inc + 1;
    inc_vector = inc*inc_vector;
    inc_vector = round(inc_vector);
    seed_voxel = init_voxel + inc_vector;
end

%Moves the seed_voxel to one loop outside the mask
init_voxel = seed_voxel; inc = 0;
while nii.img(seed_voxel(1), seed_voxel(2), seed_voxel(3)) > 0 
    inc_vector = outward_vector;
    inc = inc + 1;
    inc_vector = inc*inc_vector;
    inc_vector = round(inc_vector);
    seed_voxel = init_voxel + inc_vector;
end

% Moves the seed_voxel in to start point
init_voxel = seed_voxel; inc = 0;
while nii.img(seed_voxel(1), seed_voxel(2), seed_voxel(3)) == 0 
    inc_vector = inward_vector;
    inc = inc + 1;
    inc_vector = inc*inc_vector;
    inc_vector = round(inc_vector);
    seed_voxel = init_voxel + inc_vector;
end
end

if dir_check > 0  
%Moves inwards until contact with the mask
init_voxel = seed_voxel;
inc = 0;
while nii.img(seed_voxel(1), seed_voxel(2), seed_voxel(3)) == 0 
    inc_vector = inward_vector;
    inc = inc + 1;
    inc_vector = inc*inc_vector;
    inc_vector = round(inc_vector);
    seed_voxel = init_voxel + inc_vector;
end
end
end

function elec_struct = electrodeDesigner(nii, seed_voxels, ...
    sphere_coords, array_info)
%Finds the guide positions for electrodes
elec_struct = cell(36,1);

for k = 1:36
seed_voxel = seed_voxels(k,:);
%Finds all the spherical locations around the initial seed point
img_sph_coords = zeros(50000,3); index = 0;
for i=1:size(sphere_coords.sph_x,1)
    for j=1:size(sphere_coords.sph_x,2)
        index = index + 1;
        img_sph_coords(index,:) = ...
            [seed_voxel(1) + sphere_coords.sph_x(i,j), ...
            seed_voxel(2) + sphere_coords.sph_y(i,j), ...
            seed_voxel(3) + sphere_coords.sph_z(i,j)];
    end
end
img_sph_coords(all(img_sph_coords==0,2),:)=[];
img_sph_coords = unique(img_sph_coords, 'rows');
img_sph_coords = img_sph_coords(...
    img_sph_coords(:,1) < size(nii.img,1)+1 & img_sph_coords(:,1) > 0 &...
    img_sph_coords(:,2) < size(nii.img,2)+1 & img_sph_coords(:,2) > 0 &...
    img_sph_coords(:,3) < size(nii.img,3)+1 & img_sph_coords(:,3) > 0, :);

%Collects the filled sphere coordinates based on the binarized scalp mask
scalp_coords = zeros(50000,3); index = 0;
for i=1:size(img_sph_coords,1)
    if nii.img(img_sph_coords(i,1), img_sph_coords(i,2), ...
            img_sph_coords(i,3)) > 0
        index = index + 1;
        scalp_coords(index,:) = img_sph_coords(i,:);
    end
end
scalp_coords(all(scalp_coords==0,2),:)=[];

%Determines the normal vector and ensures it's outward
normalv = normnd(scalp_coords);
outward_normal_vector = [normalv(1), normalv(2), normalv(3)];

%Ensures that the normal vector is facing outward
center = array_info.center;
test_voxel = seed_voxel + 10*outward_normal_vector;
dist_to_seed = ((seed_voxel(1)-center(1))^2 + ...
    (seed_voxel(2)-center(2))^2 + (seed_voxel(3)-center(3))^2)^(1/2);
dist_to_test = ((test_voxel(1)-center(1))^2 + ...
    (test_voxel(2)-center(2))^2 + (test_voxel(3)-center(3))^2)^(1/2);
if dist_to_seed > dist_to_test
    outward_normal_vector = -outward_normal_vector;
end

%Determines the gel and electrode seeds as well as adjusted normal
[elec_end_point, elec_seed_voxel, gel_seed_voxel]...
    = electrodePointFinder(nii, seed_voxel, outward_normal_vector);

%Makes a strut from the normal and both seed points
elec_struct{k} = struct('normal_vector', outward_normal_vector, 'elec_seed', ...
    elec_seed_voxel, 'gel_seed', gel_seed_voxel, 'elec_end', elec_end_point);
end
end

function [unrounded_elec_end, unrounded_elec_seed, unrounded_gel_seed]...
    = electrodePointFinder(nii, seed_voxel, outward_normal_vector)
%This function is needed to move the gel base inwards in case it
%doesn't initially lie flush again the scalp. Determines the start and end
%point for both the gel and electrode cylinders

%Gets spacing infomation
spacing = nii.hdr.dime.pixdim(2);
unrounded_gel_seed = seed_voxel; 

%Makes filled circle using initial seed_voxel location
circle_points = zeros(100000,3); theta=0:0.01:2*pi; 
v=null(outward_normal_vector); new_index = 0; radius = 10/spacing;
for i=0.01:0.11:1
point = repmat(seed_voxel',1,size(theta,2))+i*radius*(v(:,1)...
    *cos(theta)+v(:,2)*sin(theta));
point1 = reshape(point(1,:), [(size(point,2)),1]);
point2 = reshape(point(2,:), [(size(point,2)),1]);
point3 = reshape(point(3,:), [(size(point,2)),1]);
points = cat(2, point1, point2, point3);

old_index = new_index;
new_index = new_index + size(points,1);
circle_points((old_index+1):new_index,:) = points(:,:);
end
circle_points(all(circle_points==0,2),:)=[];

%Counts up the number of points that are filled in the mask
circle_points = round(circle_points);
circle_points = unique(circle_points, 'rows');
scalp_point_counter = 0;

cleaned_circle_points = zeros(100000,3); index = 0;
for i=1:(size(circle_points,1))
    if circle_points(i,1) < size(nii.img,1) && ...
            circle_points(i,2) < size(nii.img,2) && ...
            circle_points(i,3) < size(nii.img,3)
        index = index + 1;
        cleaned_circle_points(index,:) = circle_points(i,:);
    end
end 
cleaned_circle_points(all(cleaned_circle_points==0,2),:)=[];

for i=1:size(cleaned_circle_points,1)
    if nii.img(cleaned_circle_points(i,1), cleaned_circle_points(i,2), ...
            cleaned_circle_points(i,3)) > 0
        scalp_point_counter = scalp_point_counter + 1;
    end
end
total_circle_points = size(cleaned_circle_points, 1);
if scalp_point_counter > 0.96*total_circle_points
    unrounded_gel_seed = seed_voxel;
end

%Repeats a loop of the above, moving the electrode face inward at each
%increment until the entire circular base is within the skin. This is the
%start of the gel mask cylinder
vector_multiplier = 0;
while scalp_point_counter < 0.99*total_circle_points

circle_points2 = zeros(100000,3);
vector_multiplier = vector_multiplier + 1;
unrounded_gel_seed = ...
    (seed_voxel - vector_multiplier*outward_normal_vector);
for i=0.05:0.05:1
point = repmat(unrounded_gel_seed',1,size(theta,2))+i*radius*(v(:,1)...
    *cos(theta)+v(:,2)*sin(theta));
point1 = reshape(point(1,:), [(size(point,2)),1]);
point2 = reshape(point(2,:), [(size(point,2)),1]);
point3 = reshape(point(3,:), [(size(point,2)),1]);
points = cat(2, point1, point2, point3);

old_index = new_index;
new_index = new_index + size(points,1);
circle_points2((old_index+1):new_index,:) = points(:,:);
end
circle_points2(all(circle_points2==0,2),:)=[];
circle_points2 = round(circle_points2);
circle_points2 = unique(circle_points2, 'rows');

cleaned_circle_points = zeros(1000000,3); index = 0;
for i=1:(size(circle_points2,1))
    if circle_points2(i,1) < size(nii.img,1) && ...
            circle_points2(i,2) < size(nii.img,2) && ...
            circle_points2(i,3) < size(nii.img,3)
        index = index + 1;
        cleaned_circle_points((index),:) = circle_points2(i,:);
    end
end 
cleaned_circle_points(all(cleaned_circle_points==0,2),:)=[];

for i=1:size(cleaned_circle_points,1)
    if nii.img(cleaned_circle_points(i,1), cleaned_circle_points(i,2), ...
            cleaned_circle_points(i,3)) > 0
        scalp_point_counter = scalp_point_counter + 1;
    end
end
total_circle_points = size(cleaned_circle_points, 1);
end

%Now moves the plane of the circle outwards until there's no filled scalp 
%in the circle to determine the electrode seed point
vector_multiplier = 0; scalp_point_counter = 1; total_circle_points = 0;
while scalp_point_counter > 0.20*total_circle_points

circle_points3 = zeros(100000,3); 
vector_multiplier = vector_multiplier + 1;
unrounded_elec_seed = ...
    unrounded_gel_seed + vector_multiplier*outward_normal_vector;
for i=0.05:0.05:1
point=repmat(unrounded_elec_seed',1,size(theta,2))+i*radius*(v(:,1)...
    *cos(theta)+v(:,2)*sin(theta));
point1 = reshape(point(1,:), [(size(point,2)),1]);
point2 = reshape(point(2,:), [(size(point,2)),1]);
point3 = reshape(point(3,:), [(size(point,2)),1]);
points = cat(2, point1, point2, point3);
old_index = new_index;
new_index = new_index + size(points,1);
circle_points3((old_index+1):new_index,:) = points(:,:);
end
circle_points3(all(circle_points3==0,2),:)=[];
circle_points3 = round(circle_points3);
circle_points3 = unique(circle_points3, 'rows');
scalp_point_counter = 0;

cleaned_circle_points = zeros(100000,3); index = 0;
for i=1:(size(circle_points3,1))
    if circle_points3(i,1) < size(nii.img,1) && ...
            circle_points3(i,2) < size(nii.img,2) && ...
            circle_points3(i,3) < size(nii.img,3)
        index = index + 1;
        cleaned_circle_points((index),:) = circle_points3(i,:);
    end
end 
cleaned_circle_points(all(cleaned_circle_points==0,2),:)=[];

for i=1:size(cleaned_circle_points,1)
    if nii.img(cleaned_circle_points(i,1), cleaned_circle_points(i,2), ...
            cleaned_circle_points(i,3)) > 0
        scalp_point_counter = scalp_point_counter + 1;
    end
end
total_circle_points = size(cleaned_circle_points, 1);
end

%Extends gel slightly further to ensure it is one structure
unrounded_elec_seed = unrounded_elec_seed + ...
    (2/spacing)*outward_normal_vector;
%Finds the end point of the electrode
unrounded_elec_end = unrounded_elec_seed + ...
    (2/spacing)*outward_normal_vector;

end

function electrodeMaker(nii, elec_struct)

spacing = nii.hdr.dime.pixdim(2);
end_gel_index = 0; end_elec_index = 0;
gel_total = zeros(50000,3); elec_total = zeros(50000,3);

for i=1:36
gel_seed = elec_struct{i}.gel_seed;
elec_seed = elec_struct{i}.elec_seed;
elec_end = elec_struct{i}.elec_end;
r = round(10/spacing);

%Generates the gel coordinates
cyl_X = zeros(1050, 1000);
cyl_Y = zeros(1050, 1000);
cyl_Z = zeros(1050, 1000);
for j=0:r
    [cyl_X((j*50+1):(j*50+50), :), cyl_Y((j*50+1):(j*50+50), :), ...
        cyl_Z((j*50+1):(j*50+50), :)] = ...
        cylinder2P((r-j), 1000, 50, gel_seed, elec_seed);
end
cyl_X(all(cyl_X==0,2),:)=[]; 
cyl_Y(all(cyl_Y==0,2),:)=[]; 
cyl_Z(all(cyl_Z==0,2),:)=[];

%Adds gel coordinates from each electrode
gel_coors = round([cyl_X(:), cyl_Y(:), cyl_Z(:)]);
gel_coors = unique(gel_coors,'rows');
start_gel_index = end_gel_index + 1;
end_gel_index = start_gel_index + size(gel_coors,1) - 1;
gel_total(start_gel_index:end_gel_index,:) = gel_coors(:,:);

%Generates the electrode coordinates
cyl_X = zeros(1050, 1000);
cyl_Y = zeros(1050, 1000);
cyl_Z = zeros(1050, 1000);
for j=1:r
    [cyl_X((j*50+1):(j*50+50), :), cyl_Y((j*50+1):(j*50+50), :), ...
        cyl_Z((j*50+1):(j*50+50), :)] = ...
        cylinder2P((r-j), 1000, 50, elec_seed, elec_end);
end
cyl_X(all(cyl_X==0,2),:)=[]; 
cyl_Y(all(cyl_Y==0,2),:)=[]; 
cyl_Z(all(cyl_Z==0,2),:)=[];

%Adds each electrode's coordinates to a list
elec_coors = round([cyl_X(:), cyl_Y(:), cyl_Z(:)]);
elec_coors = unique(elec_coors,'rows');
start_elec_index = end_elec_index + 1;
end_elec_index = start_elec_index + size(elec_coors,1) - 1;
elec_total(start_elec_index:end_elec_index,:) = elec_coors(:,:);
end

%Removes points outside the image
gel_total = gel_total(...
    (gel_total(:,1) > 0) & (gel_total(:,1) < size(nii.img,1)+1) & ...
    (gel_total(:,2) > 0) & (gel_total(:,2) < size(nii.img,2)+1) & ...
    (gel_total(:,3) > 0) & (gel_total(:,3) < size(nii.img,3)+1), :);

%Removes points from gel that intersect the scalp
for i = size(gel_total,1):-1:1
    if nii.img(gel_total(i,1), gel_total(i,2), gel_total(i,3)) == 1
        gel_total(i,:) = [];
    end
end

%Removes points outside the image
elec_total = elec_total(...
    elec_total(:,1) > 0 & elec_total(:,1) < size(nii.img,1)+1 & ...
    elec_total(:,2) > 0 & elec_total(:,2) < size(nii.img,2)+1 & ...
    elec_total(:,3) > 0 & elec_total(:,3) < size(nii.img,3)+1,:);

%Fills the gel mask
gel_mask = zeros(size(nii.img));
for i=1:size(gel_total,1)
gel_mask(gel_total(i, 1), gel_total(i, 2), gel_total(i, 3)) = 1;
end

%Fills the electrode mask
elec_mask = zeros(size(nii.img));
for i=1:size(elec_total,1)
elec_mask(elec_total(i, 1), elec_total(i, 2), elec_total(i, 3)) = 1;
end

nii.img = im2uint8(gel_mask)*255; save_untouch_nii(nii, 'mask_gel.img');
nii.img = im2uint8(elec_mask)*255; save_untouch_nii(nii, 'mask_elec.img');

end