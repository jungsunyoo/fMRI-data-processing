% Analyzing informational connectivity of GL
clear;clc;
SID = input('Subject number (1 to 24)');
%behavPath = ' '; path where onset for each condition is stored
%txtPath = ' '; root path
%projectPath = fullfile(txtPath, filesep, 'IC2');
%resultPath = fullfile(projectPath, filesep, 'IC_results');
%subjects = { }; subject IDs
sesses = {'run02', 'run03', 'run04', 'run05'}; % fMRI runs
perRun = 144; % Trials per run
nruns = 4; % number of runs
roilist = dlmread(fullfile(txtPath, filesep, 'rois.txt')); % list of ROIs to extract time series
normalVec = 1:1:580; % total TRs (145 TR per run * 4 runs)
excludevec = 1:145:580; % TRs to exclude
useVec = setdiff(normalVec, excludevec); % excluding
useVecRoi = 2:1:145;
%% Configuring FB / NFB index
temp_fb = [];
temp_nf = [];
for i = 1:6%*nruns
    a = 1+(i-1)*24:1:12+(i-1)*24;
    b = 13+(i-1)*24:1:24+(i-1)*24;
    temp_fb = [temp_fb a];
    temp_nf = [temp_nf b];
end
idx_fb = temp_fb; idx_nf = temp_nf; % index of two conditions (fb, nf)

%% Informational connectivity
for isub = 1:length(subjects)
    % structure of targetID = 1 * 576

    cd(projectPath)
    for r = 1:nruns
    % Configuration: Behav
    cd(behavPath)
    load([subjects{isub} '-fmri']) % load behavioral data
%     clearvars -except behavPath projectPath subjects perRun nruns roilist useVec isub iroi ifeat targetID idx_fb idx_nf resultPath sesses useVecRoi r
    targetID = targetID(13:end);


    targetID = targetID(useVec);
    targetID = targetID';

        targetID = targetID((r-1)*144+1:r*144);
        for iroi = 1:length(roilist)
            clear currROI
            % use 'classify' for LDA
            cd(projectPath)
            tempROI = load(['ICbetas.' subjects{isub} '.roi' sprintf('%03d', roilist(iroi)) '.' sesses{r}]);
            % load single-trial MVPA beta maps

            currROI = tempROI(:,useVecRoi);
            currROI = currROI';

            for ifeat = 1:size(currROI, 1)
                X_for_test = currROI(ifeat,:);
                X_for_train = currROI;
                X_for_train(ifeat,:) = [];
                Y_for_test = targetID(ifeat,:);
                Y_for_train = targetID;
                Y_for_train(ifeat,:) = [];

                % performing LDA
                discr = fitcdiscr(X_for_train, Y_for_train);
                class = discr.predict(X_for_test);
                predictedClass(1,ifeat) = class;
                ansClass(1,ifeat) = Y_for_test;
                %ifeat
            end
            %getting accuracy
            ansLDA = (predictedClass == ansClass);
            accuracy(iroi,isub) = mean(ansLDA);
            % getting answer pattern
            pattern(:,iroi) = predictedClass';

            disp(['Subject ' subjects{isub} ' Run ' sesses{r} ' ROI number ' num2str(roilist(iroi)) ' completed.'])
        end

        % now getting informational connectivity
        % First: dividing into FB / NFB
        % Also divide into four runs

        pattern_fb = pattern(idx_fb,:);
        pattern_nf = pattern(idx_nf,:);

        for i = 1:size(pattern,2)
            for j = 1:size(pattern,2)
                if j>=i
                    fb1 = pattern_fb(:,i); fb2 = pattern_fb(:,j);
                    nf1 = pattern_nf(:,i); nf2 = pattern_nf(:,j);
                    corr_fb = mean((fb1 == fb2));
                    corr_nf = mean((nf1 == nf2));
                    rval_fb = corr(fb1, fb2);
                    rval_nf = corr(nf1, nf2);
                    corrMatrix_fb(i,j,r) = corr_fb;
                    corrMatrix_nf(i,j,r) = corr_nf;
                    rMatrix_fb(i,j,r) = rval_fb;
                    rMatrix_nf(i,j,r) = rval_nf;
                    corrMatrix_fb(j,i,r) = corr_fb;%corrMatrix_fb(j,j,r);
                    corrMatrix_nf(j,i,r) = corr_nf;%corrMatrix_nf(i,j,r);
                    rMatrix_fb(j,i,r) = rval_fb;%rMatrix_fb(i,j,r);
                    rMatrix_nf(j,i,r) = rval_nf;%rMatrix_nf(i,j,r);
                end
            end
        end

    end
    cd(resultPath)
    save([subjects{isub} '_IC'], 'corrMatrix_fb', 'corrMatrix_nf', 'rMatrix_fb', 'rMatrix_nf', 'pattern', 'accuracy')

end
