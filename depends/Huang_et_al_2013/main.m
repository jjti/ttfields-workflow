% This script will call Abaqus FEM solver to automatically solve the FEM
% for all possible bipolar electrode configurations with one fixed cathode.
%
% You will need the 'en' and 'ee' Matlab files generated from meshProcess.m
%
% After running through, Abaqus will save its solution as .odb file for
% each electrode configuration.
%
% (c) Jacek Dmochowski, March 2011
% The Neural Engineering Lab, Dept. of Biomedical Engineering, City College of New York
% Send bugs to yhuang16@citymail.cuny.edu

clear all; close all; clc

% PARAMETERS
load en;  % cell array of nodes belonging to each electrode
load ee;  % cell array of elements belonging to each electrode
mkdir('./scratch') % temporary directory used by Abaqus
nsecs_to_pause=60*60; % waiting time for Matlab when Abaqus is working
pause on;

for i = 1:73 % labels of anodes (if 74 electrodes placed, and #74 is cathode)
    
    ib = 74; % label of grounded electrode (cathode)
    is = i % labels of anodes
    filename_in='enter_your_mesh_file_name_here.inp'; % input mesh file name
    filename_out=['enter_your_mesh_file_name_here_' num2str(ib) '_' num2str(is)  '.inp']; % mesh file submitted to Abaqus for solving
    
    sigma=[   ...
        2.5e-14; ...  % air
        0.465; ... % skin
        0.01; ...     % bone
        0.276; ...    % gray
        0.126; ...    % white
        1.65; ...     % csf
        5.9e7; ... % elec
        0.3 ; ... % gel
        ];
    % conductivities for all tissue types, coming in the priority order in ScanIP
    sigma=sigma(end:-1:1); cc=1;
    % because in the mesh file (inp), the 8 tissues come in an inversed order of priority order in ScanIP
    
    %  STEPS 1 & 2 -- add material properties + change all element types to thermal electric
    fin = fopen(filename_in);
    fout = fopen(filename_out,'w');
    while ~feof(fin)
        s = fgetl(fin);
        s = strrep(s, 'C3D4', 'DC3D4E'); % Solving type: thermal electric element
        fprintf(fout,'%s \n',[s]);
        if length(s)>=9
            if strcmp(s(1:9),'*MATERIAL')
                s2='*ELECTRICAL CONDUCTIVITY';
                s3=num2str(sigma(cc)); cc=cc+1;
                fprintf(fout,'%s \n %s \n',s2,s3); % write conductivity values to each tissue type in .inp file
            end
        end
    end
    fclose(fin);
    fclose(fout);
    
    % STEP 3: add boundary condition
    % define set of nodes in grounded electrode
    M='*NSET, NSET=GELEC';
    dlmwrite(filename_out, M, 'delimiter', '','-append');
    
    dlmwrite(filename_out, en{ib}, 'delimiter', ',','-append','precision','%7d');
    
    M='*BOUNDARY';
    dlmwrite(filename_out, M, 'delimiter', '','-append');
    
    M='GELEC,9,9';
    dlmwrite(filename_out, M, 'delimiter', '','-append');
    
    % STEP 4: add analysis step
    % create set of elements that will be energized
    M='*ELSET, ELSET=SELEC';
    dlmwrite(filename_out, M, 'delimiter', '','-append');
    dlmwrite(filename_out, ee{is}, 'delimiter', ',','-append','precision','%7d');
    
    % create surface
    M='*SURFACE, NAME=SSURF, TYPE=ELEMENT';
    dlmwrite(filename_out, M, 'delimiter', '','-append');
    
    M='SELEC,';
    dlmwrite(filename_out, M, 'delimiter', '','-append');
    
    % add step with load
    M='*STEP,NAME=STEP-1,AMPLITUDE=STEP,UNSYMM=NO'; % quasi-static
    dlmwrite(filename_out, M, 'delimiter', '','-append');
    M='*Coupled Thermal-Electrical, steady state, deltmx=0.';
    dlmwrite(filename_out, M, 'delimiter', '','-append')
    M='1., 1., 1e-05, 1.,  ';
    dlmwrite(filename_out, M, 'delimiter', '','-append')
    
    M='*SOLUTION TECHNIQUE,TYPE=SEPARATED';
    dlmwrite(filename_out, M, 'delimiter', '','-append')
    M='*Dsecurrent';
    dlmwrite(filename_out, M, 'delimiter', '','-append')
    M=' SSURF, CS, 1'; % injecting current density of 1 A/m^2
    dlmwrite(filename_out, M, 'delimiter', '','-append')
    M='*Restart, write, frequency=0';
    dlmwrite(filename_out, M, 'delimiter', '','-append')
    M='*Output, field, frequency=99999';
    dlmwrite(filename_out, M, 'delimiter', '','-append')
    M='*Node Output';
    dlmwrite(filename_out, M, 'delimiter', '','-append')
    M=' COORD,';
    dlmwrite(filename_out, M, 'delimiter', '','-append')
    M='*Element Output, directions=YES';
    dlmwrite(filename_out, M, 'delimiter', '','-append')
    M='EPG,'; % output EPG (electric potential gradient, i.e., electric field)
    dlmwrite(filename_out, M, 'delimiter', '','-append')
    M='*Output, history, variable=PRESELECT';
    dlmwrite(filename_out, M, 'delimiter', '','-append')
    M='*END STEP';
    dlmwrite(filename_out, M, 'delimiter', '','-append')
    
    % STEP 5: call Abaqus to solve the model
    system(['/usr/local/abaqus/Commands/abaqus job=' filename_out ' scratch=./scratch cpus=12'])
             % the directory of Abaqus                                      % the directory of 'scratch' folder
                                                                            % you can specify how many CPU cores can be used for the computation
    
    fprintf('waiting until solution completes...\n');
    pause(nsecs_to_pause);
    
end
