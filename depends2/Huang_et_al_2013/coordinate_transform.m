function coordinate_transform(excelname,outputname)

% coordinate_transform(excelname,outputname)
%
% This function transforms the spherical coordinates of electrodes from
% EasyCap (http://www.easycap.de/easycap/e/downloads/M1_ThetaPhi.htm) into
% Cartesian coordinates that will be used by electrode_placement.m
%
% You can save the coordiantes from EasyCap as an MS excel file and edit it
% as you wish. This function then reads the excel file and output as a .mat
% file.
%
% (c) Yu Huang (Andy), November 2010
% The Neural Engineering Lab, Dept. of Biomedical Engineering, City College of New York
% Send bugs to yhuang16@citymail.cuny.edu

coor_data = importdata(excelname);

if isfield(coor_data.data,'Sheet1')
    origin_coor = coor_data.data.Sheet1;
else
    origin_coor = coor_data.data;
end
origin_coor = origin_coor*pi/180; % spherical coordinates
elec_name = coor_data.textdata.Sheet1; % electrode names

X = sin(origin_coor(:,1)).*cos(origin_coor(:,2));
Y = sin(origin_coor(:,1)).*sin(origin_coor(:,2));
Z = cos(origin_coor(:,1));
% converting to Cartesian coordinates

% figure; plot3(X,Y,Z);
% hold on, plot3(X,Y,Z,'r.'); hold off

if length(elec_name) ~= length(X)
    offset = length(elec_name) - length(X);
    if offset<0
        error('Electrode labels appear shorter than number of coordinates in input file.');
    %else disp('Input file appers to have a header.')
    end
else
    offset = 0;
end

D = struct('labels',{},'X',{},'Y',{},'Z',{}); % output structure
for i = 1:(length(elec_name)-offset)
    D(i).labels = elec_name{i+offset};
    D(i).X = X(i);
    D(i).Y = Y(i);
    D(i).Z = Z(i);
end

save(outputname,'D');