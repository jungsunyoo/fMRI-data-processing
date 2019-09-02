%% 1. Defining subjects and pathway
clear; clc
warning('off', 'all')
allSubjects = {'fmri121'};
projectPath =  '/Volumes/JS/TempMem/topup_data';
sesses = {'RET-01', 'RET-02', 'RET-03', 'RET-04', 'RET-05'};
addpath(genpath('/Volumes/JS/spm12'));

%% 2. Gunzip
for isub = 1:length(allSubjects)
    currPath = fullfile(projectPath, filesep, allSubjects{isub}, filesep, 'day1', filesep, 'retrieval');
    for nrun = 1:length(sesses)%6 % for encoding
        clear matlabbatch
        cd(currPath)
        if ~exist([sesses{nrun} '.nii'])
            gunzip([sesses{nrun} '.nii.gz']);
            disp(['gunzip for ' allSubjects{isub} ' ' sesses{nrun}  ' done']);
        end
        if ~exist(sesses{nrun})
            mkdir(sesses{nrun})
        end
        movefile([sesses{nrun} '.nii'], fullfile(currPath, sesses{nrun}))
        cd(fullfile(currPath, sesses{nrun}))
        spm_file_split([sesses{nrun} '.nii']);
        disp(['Splitting for subject ' num2str(isub) ' run ' num2str(nrun) ' is finished']);
        delete([sesses{nrun} '.nii']);
    end
end

%% 3. SLICE-TIMING
tr = 0.8;
nslice = 60;
TR_temp = tr*1000;
sliceorder = [1:2:11; 13:2:23; 25:2:35; 37:2:47; 49:2:59; 2:2:12; 14:2:24; 26:2:36; 38:2:48; 50:2:60]'; % (with rows scanned simultaneously)
slicetimes = zeros(1,60);
slicetimes(sliceorder(1,:)) = 0 : TR_temp/10 : TR_temp-TR_temp/10; % number of slices (60) / MB factor (6) = 10
slicetimes(sliceorder(2,:)) = 0 : TR_temp/10 : TR_temp-TR_temp/10;
slicetimes(sliceorder(3,:)) = 0 : TR_temp/10 : TR_temp-TR_temp/10;
slicetimes(sliceorder(4,:)) = 0 : TR_temp/10 : TR_temp-TR_temp/10;
slicetimes(sliceorder(5,:)) = 0 : TR_temp/10 : TR_temp-TR_temp/10;
slicetimes(sliceorder(6,:)) = 0 : TR_temp/10 : TR_temp-TR_temp/10;
for isub = 1:length(allSubjects)
    clear matlabbatch;
    for irun=1:length(sesses)
        clear tempFiles;
    currPath = fullfile(projectPath, filesep, allSubjects{isub}, filesep, 'day1', filesep, 'retrieval', filesep, sesses{irun});%'/Users/mnd_mac/Desktop/TempMem_YJS/topup_data';
    delete(fullfile(currPath, filesep,'ra*.nii'));
        Filter = 'RET*-*';%Filter = [sesses{irun} '.nii']; % this reads in all valid text file names (which are not 4D!!!!)
        tempFiles = dir(fullfile(currPath, Filter));
        tempFiles = char({tempFiles.name}');
        tempFiles = [repmat(fullfile(currPath, filesep), size(tempFiles,1),1) tempFiles];
        
        matlabbatch{1}.spm.temporal.st.scans{irun} = cellstr(tempFiles);
        matlabbatch{1}.spm.temporal.st.nslices = nslice;
        matlabbatch{1}.spm.temporal.st.tr = tr;
        matlabbatch{1}.spm.temporal.st.ta = tr-(tr/nslice);
        matlabbatch{1}.spm.temporal.st.so = slicetimes;
        matlabbatch{1}.spm.temporal.st.refslice = tr*1000/2;
        matlabbatch{1}.spm.temporal.st.prefix = 'a';
    end
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
end

%% 4. Realign 
Filter = 'a*'; % this reads in all valid text file names
currPath = projectPath;%'/Users/mnd_mac/Desktop/TempMem_YJS/topup_data';
for isub = 1:length(allSubjects)%4:length(allSubjects)
    clear matlabbatch tempFiles1  tempFiles2  tempFiles3  tempFiles4  tempFiles5  tempFiles6 
    for irun = 1:length(sesses)
        tempdir = fullfile(currPath, filesep, allSubjects{isub},filesep, 'day1', filesep, 'retrieval', filesep, sesses{irun});
        tempFiles = dir(fullfile(tempdir,filesep, Filter));
        tempFiles = char({tempFiles.name});
        tempFiles = [repmat(fullfile(tempdir, filesep), size(tempFiles,1),1) tempFiles];
        eval(['tempFiles' num2str(irun) '=tempFiles;']);
    end
        matlabbatch{1}.spm.spatial.realign.estwrite.data = {cellstr(tempFiles1), cellstr(tempFiles2), cellstr(tempFiles3), cellstr(tempFiles4), cellstr(tempFiles5)};
        
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);

    
for irun = 1:length(sesses) %this deletes slicetimed/not realigned data in order to save disk space)
        delete(fullfile(currPath, filesep,'a*nii'));
        disp(['deleted ' allSubjects{isub} ' ' num2str(irun)])       
end
end
%%
origT1path = projectPath;%'/Volumes/MGTEC/TempMem_topup/topup_data';
for i = 1:length(allSubjects)

origpath = [origT1path filesep allSubjects{i} filesep 'day1/encoding'];
destpath = [fullfile(projectPath, filesep, allSubjects{i}, '/day1/encoding')];

cd(origpath)
filename = dir(fullfile(origpath, '*T1*.nii'));
filename = char(filename.name);
copyfile(filename, destpath)
end

%% 5. Corregistration (only available after day2 retrieval)
for isub = 1:length(allSubjects)%4:length(allSubjects)
    clear matlabbatch;
    currPath = projectPath;%'/Users/mnd_mac/Desktop/TempMem_YJS/topup_data';
    corregPath = fullfile(currPath, filesep, allSubjects{isub}, filesep, 'day2');
    refPath = fullfile(currPath, filesep, allSubjects{isub},filesep, 'day1', filesep, 'retrieval', filesep, sesses{1});

    refFilter = ['mean*'];
    corregFilter = ['*T1*'];
    t1dir = dir(fullfile(corregPath, corregFilter));
    t1File = char({t1dir.name}');
    T1source = fullfile(corregPath, filesep, t1File);
    refdir = dir(fullfile(refPath, refFilter));
    refFile = char({refdir.name}');
    refsource = fullfile(refPath, filesep, refFile);
    for irun = 1:length(sesses)
        currPath = fullfile(projectPath, filesep, allSubjects{isub}, filesep, 'day1\encoding', filesep, sesses{irun});
        matlabbatch{1}.spm.spatial.coreg.estimate.ref = cellstr(refsource);
        matlabbatch{1}.spm.spatial.coreg.estimate.source = cellstr(T1source);
        matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    end
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
    
end
%% 6. Segmentation
matlabdir = '/Volumes/JS/spm12';%'/Users/mnd_mac/Desktop/Jinyoung/spm12';
for isub = 1:length(allSubjects)
    clear matlabbatch T1source;
    corregPath = fullfile(currPath, filesep, allSubjects{isub}, filesep, 'day2');
    corregFilter = ['*T1*'];
    t1dir = dir(fullfile(corregPath, corregFilter));
    t1File = char({t1dir.name}');
    T1source = fullfile(corregPath, filesep, t1File);
    matlabbatch{1}.spm.spatial.preproc.channel.vols = cellstr(T1source);
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'/Volumes/JS/spm12/tpm/TPM.nii,1'};
    matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 1];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'/Volumes/JS/spm12/tpm/TPM.nii,2'};
    matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 1];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'/Volumes/JS/spm12/tpm/TPM.nii,3'};
    matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'/Volumes/JS/spm12/tpm/TPM.nii,4'};
    matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'/Volumes/JS/spm12/tpm/TPM.nii,5'};
    matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'/Volumes/JS/spm12/tpm/TPM.nii,6'};
    matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'eastern';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
end

%% 7. Normalization
for isub = 1:length(allSubjects)
    clear matlabbatch t1dir t1Path t1File T1source;
    normPath = fullfile(currPath, filesep, allSubjects{isub}, filesep, 'day2');%fullfile(projectPath, filesep, allSubjects{isub}, filesep, 'day2');
    normFilter = ['y_*'];
    t1dir = dir(fullfile(normPath, normFilter));
    t1File = char({t1dir.name}');
    T1source = fullfile(normPath, filesep, t1File);
    epiFilter = 'ra*.nii';
  
    allepis = [];

    for irun = 1:length(sesses)
        runepis       = dir(fullfile(projectPath,allSubjects{isub},'day1', filesep, 'retrieval', filesep, sesses{irun}, epiFilter));
        runepis       = [char({runepis.name}') ];
        runepis       = [repmat([fullfile(projectPath,allSubjects{isub},'day1', filesep, 'retrieval', filesep, sesses{irun}) filesep],size(runepis,1),1) runepis];
        runepis = cellstr(runepis);
        allepis       = [allepis;runepis];
    end
    allepis
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(T1source);
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = allepis;
    
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-90 -126 -72
        90 90 108];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
end
%% 8. Smoothing
epiFilter = 'wra*.nii';
for isub = 1:length(allSubjects) % (2019.06.12. 01:17am: run isub=19 later)
    clear matlabbatch;
    allepis = [];
    for irun = 1:length(sesses)
        runepis       = dir(fullfile(projectPath,allSubjects{isub},'day1', filesep, 'retrieval', filesep, sesses{irun},filesep, epiFilter));
        runepis       = [char({runepis.name}') ];
        runepis       = [repmat([fullfile(projectPath,allSubjects{isub},'day1', filesep, 'retrieval', filesep, sesses{irun}) filesep],size(runepis,1),1) runepis];
        runepis = cellstr(runepis);
        allepis       = [allepis;runepis];
    end
    allepis
    
    matlabbatch{1}.spm.spatial.smooth.data = allepis;
    matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
    
    for irun = 1:length(sesses) %this deletes slicetimed/not realigned data in order to save disk space)
        delete(fullfile(projectPath,allSubjects{isub},['day1\encoding\Run' sprintf('%01d',irun)],'wua*nii'));
        disp(['deleted ' allSubjects{isub} ' ' num2str(irun)])
    end
end

