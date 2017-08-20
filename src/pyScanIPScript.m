function pyScanIPScript(info, dir)

template = fullfile(spm('Dir'),'Toolbox/TTFs','ScanIP_Template.txt');

cd(dir);
i1 = num2str(info(1),4);
i2 = num2str(info(2),4);
i3 = num2str(info(3),4);
i4 = num2str(info(4), '%d');
i5 = num2str(info(5), '%d');
i6 = num2str(info(6), '%d');
zero = num2str(0, '%d');

% Converting the directory 
dir = regexprep(dir, '\', '\\\');

% Creating the sizing and cropping input
% I am 90% sure this is a bug that should be fixed in ScanIP
% having to offset by approximately half the x-dimension
sizing = [i4 ',' i5 ',' i6 ',' i1 ',' i2 ',' i3 ',' zero ','];
sizing = eval(mat2str(sizing));
cropping = [zero ',' zero ',' zero ',' i4 ',' i5 ',' i6 '))'];
cropping = eval(mat2str(cropping));

% Taking in all the lines from the ScanIP template
scriptLines = fopen(template, 'r');
ScanIPScript = textscan(scriptLines, '%s', 'Delimiter', '\n');
fclose(scriptLines);

% Adjusts the template to the specifics of the segmented image
ScanIPScript{1}{4} = ['"', dir, '\\mask_elec.img"'];
ScanIPScript{1}{11} = ['"', dir, '\\mask_gel.img"'];
ScanIPScript{1}{18} = ['"', dir, '\\mask_scalp.img"'];
ScanIPScript{1}{25} = ['"', dir, '\\mask_orbits.img"'];
ScanIPScript{1}{32} = ['"', dir, '\\mask_bone.img"'];
ScanIPScript{1}{39} = ['"', dir, '\\mask_csf.img"'];
ScanIPScript{1}{46} = ['"', dir, '\\mask_cere.img"'];
ScanIPScript{1}{53} = ['"', dir, '\\mask_gray.img"'];
ScanIPScript{1}{60} = ['"', dir, '\\mask_white.img"'];
if exist('nii_for_seg.nii', 'file')
    ScanIPScript{1}{67} = ['"', dir, '\\nii_for_seg.nii"'];
else
    ScanIPScript{1}{67} = ['"', dir, '\\nii_for_seg.img"'];
end

for i = 6:7:69
   ScanIPScript{1}{i,:} = sizing; 
end

for i = 8:7:71
   ScanIPScript{1}{i,:} = cropping; 
end

ScanIPScript{1}{221} = ['"', dir, '\\Masks.sip"'];

% Prints out each line to a text file
nrows = size(ScanIPScript{1}, 1);
FormatSpec = '%s\r\n';
outgoingScript = fopen('Script_for_ScanIP.txt', 'w');
for i = 1:nrows
    fprintf(outgoingScript, FormatSpec, ScanIPScript{1}{i, :}); 
end
fclose(outgoingScript);

copyfile('Script_for_ScanIP.txt', 'Script_for_ScanIP.py');
delete('Script_for_ScanIP.txt');

end