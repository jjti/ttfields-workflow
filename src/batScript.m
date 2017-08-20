function batScript(fnm, ScanIP_Loc)

% Optional function for creating a bat script for fully automated import to ScanIP.
% Input: <Folder with segmentation and python script>, <ScanIP system
% location>

name = 'Script_for_ScanIP.py';

fileID = fopen('Upload2ScanIP.txt', 'w');

fprintf(fileID, '%s\r\n', '@ECHO OFF');
fprintf(fileID, 'PATH "%s"\r\n', ScanIP_Loc);
formatSpec = 'start ScanIP --run-script="%s\\%s" --script-lan=python \r\n';
fprintf(fileID, formatSpec, fnm, name);
fprintf(fileID, '%s\r\n', 'END');
fprintf(fileID, '%s', 'EXIT');

fclose(fileID);

copyfile('Upload2ScanIP.txt', 'Upload2SCNIP.bat');
delete('Upload2ScanIP.txt');

% Starts the program via ScanIP command
first = 'start Upload2SCNIP.bat';
system(first);

end