
function segment_input(vol_path,spm_path,script_path,volume, reg, fwhm)
disp(['vol path',vol_path]);
disp(['spm path',spm_path]);
disp(['script path',script_path]);
disp(['volume',volume]);
disp(['reg',reg]);
disp(['fwhm',fwhm]);
% spm path
addpath(spm_path)
%script path
addpath(script_path)
% directory of the folder
currDir = strcat(vol_path,volume,'/')
%cd(currDir)
%ls
% input strings to number 
reg = str2num(reg)
fwhm = str2num(fwhm)

%%%%%%%%%%%%%%%%%%%%%%%%
%%% SPM JOB MANAGER  %%%
%%%%%%%%%%%%%%%%%%%%%%%%
%initialize job manager
spm('defaults', 'FMRI');
spm_jobman('initcfg');

% tissue values, (standard)
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {strcat(spm_path,'tpm/TPM.nii,1')};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {strcat(spm_path,'tpm/TPM.nii,2')};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {strcat(spm_path,'tpm/TPM.nii,3')};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {strcat(spm_path,'tpm/TPM.nii,4')};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {strcat(spm_path,'tpm/TPM.nii,5')};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {strcat(spm_path,'tpm/TPM.nii,6')};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];

%warp values, standard
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];

% channel values
matlabbatch{1}.spm.spatial.preproc.channel.vols = {strcat(currDir,'highres_head.nii,1')};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = reg;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = fwhm;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];

% run spm job
disp(['Volume: ',currDir]);
disp(['spm_jobman with: biasreg=',num2str(reg),', biasfwhm=',num2str(fwhm)]);
spm_jobman('run', matlabbatch);
            
%rename files respectivly
disp('Renaming...');
files = dir(currDir); % all files in current folder
for f=1:length(files)
	if regexp(files(f).name,'^c') %filename starting with 'c'
	    oldname = files(f).name;
	    newname = strcat(oldname(2),'_',sprintf('%g_%03d.nii',reg,fwhm));
	    disp([oldname,' to ',newname]);
	    movefile(strcat(currDir,oldname), strcat(currDir,newname));
	end
end % end renaming

% exit matlab when the job is done
disp('segment finished')
exit
end % end function/ script
