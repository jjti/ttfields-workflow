function spmSegment(nii_Img)

n = [1 0];
w = [0 0];

WM = fullfile(spm('Dir'),'toolbox/TTFields-Workflow/tpms','WM.nii');
GM = fullfile(spm('Dir'),'toolbox/TTFields-Workflow/tpms','GM.nii');
CSF = fullfile(spm('Dir'),'toolbox/TTFields-Workflow/tpms','CSF.nii');
Cere = fullfile(spm('Dir'),'toolbox/TTFields-Workflow/tpms','SUIT.nii');
Bone = fullfile(spm('Dir'),'toolbox/TTFields-Workflow/tpms','Bone.nii');
Scalp = fullfile(spm('Dir'),'toolbox/TTFields-Workflow/tpms','Scalp.nii');
Orbitals = fullfile(spm('Dir'),'toolbox/TTFields-Workflow/tpms','Orbitals.nii');
Air = fullfile(spm('Dir'),'toolbox/TTFields-Workflow/tpms','Air.nii');

spm_jobman('initcfg');

matlabbatch{1}.spm.tools.preproc8.channel.vols = {nii_Img};
matlabbatch{1}.spm.tools.preproc8.channel.biasreg = 0.1;
matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm = 60;
matlabbatch{1}.spm.tools.preproc8.channel.write = [1 1];
matlabbatch{1}.spm.tools.preproc8.tissue(1).tpm = {GM};
matlabbatch{1}.spm.tools.preproc8.tissue(1).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(1).native = n;
matlabbatch{1}.spm.tools.preproc8.tissue(1).warped = w;
matlabbatch{1}.spm.tools.preproc8.tissue(2).tpm = {WM};
matlabbatch{1}.spm.tools.preproc8.tissue(2).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(2).native = n;
matlabbatch{1}.spm.tools.preproc8.tissue(2).warped = w;
matlabbatch{1}.spm.tools.preproc8.tissue(3).tpm = {CSF};
matlabbatch{1}.spm.tools.preproc8.tissue(3).ngaus = 3;
matlabbatch{1}.spm.tools.preproc8.tissue(3).native = n;
matlabbatch{1}.spm.tools.preproc8.tissue(3).warped = w;
matlabbatch{1}.spm.tools.preproc8.tissue(4).tpm = {Bone};
matlabbatch{1}.spm.tools.preproc8.tissue(4).ngaus = 4;
matlabbatch{1}.spm.tools.preproc8.tissue(4).native = n;
matlabbatch{1}.spm.tools.preproc8.tissue(4).warped = w;
matlabbatch{1}.spm.tools.preproc8.tissue(5).tpm = {Scalp};
matlabbatch{1}.spm.tools.preproc8.tissue(5).ngaus = 4;
matlabbatch{1}.spm.tools.preproc8.tissue(5).native = n;
matlabbatch{1}.spm.tools.preproc8.tissue(5).warped = w;
matlabbatch{1}.spm.tools.preproc8.tissue(6).tpm = {Cere};
matlabbatch{1}.spm.tools.preproc8.tissue(6).ngaus = 4;
matlabbatch{1}.spm.tools.preproc8.tissue(6).native = n;
matlabbatch{1}.spm.tools.preproc8.tissue(6).warped = w;
matlabbatch{1}.spm.tools.preproc8.tissue(7).tpm = {Orbitals};
matlabbatch{1}.spm.tools.preproc8.tissue(7).ngaus = 4;
matlabbatch{1}.spm.tools.preproc8.tissue(7).native = n;
matlabbatch{1}.spm.tools.preproc8.tissue(7).warped = w;
matlabbatch{1}.spm.tools.preproc8.tissue(8).tpm = {Air};
matlabbatch{1}.spm.tools.preproc8.tissue(8).ngaus = 4;
matlabbatch{1}.spm.tools.preproc8.tissue(8).native = n;
matlabbatch{1}.spm.tools.preproc8.tissue(8).warped = w;
matlabbatch{1}.spm.tools.preproc8.warp.mrf = 1;
matlabbatch{1}.spm.tools.preproc8.warp.reg = 4;
matlabbatch{1}.spm.tools.preproc8.warp.affreg = 'mni'; 
matlabbatch{1}.spm.tools.preproc8.warp.samp = 3;
matlabbatch{1}.spm.tools.preproc8.warp.write = [1 1]; 

spm_jobman('run',matlabbatch);

end
