function meshProcess(dirname,fname,elecLoc,gelLoc)

% meshProcess(dirname,fname,elecLoc,gelLoc)
%
% Process the mesh file (.inp) generated from Simpleware.
% It reads the mesh file, labels the electrodes for solving in Abaqus, and
% labels the gel to calculate EXACT energized area for each electrode.
%
% Input:
% dirname: directory of the mesh file;
% fname: name of the mesh file;
% elecLoc: the location of electrode in the .inp file;
% gelLoc: the location of gel in the .inp file;
%
% Output:
% en.mat: nodes belonging to each electrode, used to define grounded
% electrode (cathode);
% ee.mat: elements belonging to each electrode, used to define stimulating
% electrode (anode);
% elecArea.mat: EXACT energized area for each electrode, used to PRECISELY
% scale the Abaqus solution of each electrode pair;
% 
% The 'en' and 'ee' will be used in the script 'main.m' for automatically
% solving the FEM for all possible bipolar electrode configurations.
% The 'elecArea' will be used in the optimization script for computation of
% optimal electrode montage for tDCS.
%
% (c) Yu Huang (Andy), February 2015
% (c) Yu Huang (Andy), March 2014
% (c) Yu Huang (Andy), March 2011
% (c) Jacek Dmochowski, March 2011
% The Neural Engineering Lab, Dept. of Biomedical Engineering, City College of New York
% Send bugs to yhuang16@citymail.cuny.edu

% extract mesh from .inp file
cd(dirname)
fid=fopen(fname);

% find beginning of node listing, and meanwhile count how many tissue masks
% are there in the model
numOfTissue = 0;
mystr='';
while strcmp(mystr,'*NODE')==0
    C = textscan(fid, '%5c', 1);
    mystr=C{1};
    if strcmp(mystr,'*MATE')==1
        numOfTissue = numOfTissue+1;
    end
    fgetl(fid);
end
fprintf('Found beginning of node listings.\n');
% now read in nodes
fprintf('reading in nodes...\n');
fgetl(fid);
fgetl(fid);
fgetl(fid);
N = textscan(fid, '%f %f %f %f', 'Delimiter', ',');

E = cell(numOfTissue,1);
for k = 1:numOfTissue
    mystr='';
    while strcmp(mystr,'*ELEMENT')==0
        C = textscan(fid, '%8c', 1);
        mystr=C{1};
        fgetl(fid);
    end
    fprintf('Found beginning of element listings for Tissue #%d.\n',k);
    % now read in elements
    fprintf('reading in elements for Tissue #%d...\n',k);
    E{k} = cell2mat(textscan(fid, '%d %d %d %d %d', 'Delimiter', ','));
end

fclose(fid);

nodes=cell2mat(N);
nodes(:,2:4) = nodes(:,2:4)+1; % for 1 mm model

element_electrode = E{elecLoc};

% generate coordinates for each electrode element
X = zeros(size(element_electrode,1),3);
for e = 1:size(element_electrode,1)
    X(e,:) = mean ( nodes(element_electrode(e,2:5),2:4) )';  % centroid of element
end

% load the range of coordinates for each electrode
load rnge_elec.mat;
% labelling electrodes
disp('labeling electrodes...')
label_elec = zeros(size(X,1),1);
for i = 1:size(X,1)
    offset = 0;
    while label_elec(i)==0
        for k = 1:length(rnge)
            if ~isempty(rnge{k}) && X(i,1)>rnge{k}(2,1)-offset && X(i,1)<rnge{k}(1,1)+offset && X(i,2)>rnge{k}(2,2)-offset && X(i,2)<rnge{k}(1,2)+offset && X(i,3)>rnge{k}(2,3)-offset && X(i,3)<rnge{k}(1,3)+offset
                label_elec(i) = k;
            end
        end
        offset = offset+0.1;
    end
    if mod(i,100) == 0
        fprintf('%f%% completed...\n',i*100/length(X))
    end
end
disp('labeling electrodes DONE!')

save label_elec.mat label_elec

% data preparation
en = cell(1,length(rnge));
ee = cell(1,length(rnge));
for n = 1:length(rnge)
    if ~isempty(find(label_elec==n))
        en{n}=nodes(element_electrode(label_elec==n,2:5),1); % nodes belonging to each electrode, used to define grounded electrode
        ee{n}=element_electrode(label_elec==n,1); % elements belonging to each electrode, used to define stimulating electrode
    end
end
save en.mat en
save ee.mat ee

element_gel = E{gelLoc};

% generate coordinates for each gel element
X = zeros(size(element_gel,1),3);
for g = 1:size(element_gel,1)
    X(g,:) = mean ( nodes(element_gel(g,2:5),2:4) )';  % centroid of element
end

% load the range of coordinates for each gel
load rnge_gel.mat;
% labelling gel
disp('labeling gel...')
label_gel = zeros(size(X,1),1);
for i = 1:size(X,1)
    offset = 0;
    while label_gel(i)==0
        for k = 1:length(rnge)
            if ~isempty(rnge{k}) && X(i,1)>rnge{k}(2,1)-offset && X(i,1)<rnge{k}(1,1)+offset && X(i,2)>rnge{k}(2,2)-offset && X(i,2)<rnge{k}(1,2)+offset && X(i,3)>rnge{k}(2,3)-offset && X(i,3)<rnge{k}(1,3)+offset
                label_gel(i) = k;
            end
        end
        offset = offset+0.1;
    end
    if mod(i,100) == 0
        fprintf('%f%% completed...\n',i*100/length(X))
    end
end
disp('labeling gel DONE!')

save label_gel.mat label_gel

% area calculation
disp('calculating energized area for each electrode...')
elecArea = zeros(length(rnge),1);
element_electrode = double(element_electrode);
element_gel = double(element_gel);
for n = 1:length(elecArea)
    % compute surface faces and vertices for each electrode and gel
    tri = element_electrode(label_elec==n,2:5);
    trep = TriRep(tri,nodes(:,2:4));
    [faces_elec verts_elec]=freeBoundary(trep);
    
    tri = element_gel(label_gel==n,2:5);
    trep = TriRep(tri,nodes(:,2:4));
    [faces_gel verts_gel]=freeBoundary(trep);
    
    [~,iE,iG] = intersect(verts_elec,verts_gel,'rows');
    tempTag = ismember(faces_elec,iE);
    
    faces_overlap = faces_elec(sum(tempTag,2)==3,:); % find the overlap faces

    % calculate the total surface area of each electrode
    a = verts_elec(faces_elec(:, 2), :) - verts_elec(faces_elec(:, 1), :);
    b = verts_elec(faces_elec(:, 3), :) - verts_elec(faces_elec(:, 1), :);
    c = cross(a, b, 2);
    totalElecArea = sum(0.5*sqrt(sum(c.^2, 2)));
    % calculate the area of the overlap between electrode and gel
    a = verts_elec(faces_overlap(:, 2), :) - verts_elec(faces_overlap(:, 1), :);
    b = verts_elec(faces_overlap(:, 3), :) - verts_elec(faces_overlap(:, 1), :);
    c = cross(a, b, 2);
    elecGelOverlapArea = sum(0.5*sqrt(sum(c.^2, 2)));
    
    elecArea(n) = totalElecArea-elecGelOverlapArea;
    % the EXACT energized area is the total electrode surface area minus the area of the overlap between electrode and gel
end
save elecArea.mat elecArea

disp('DONE!!!')