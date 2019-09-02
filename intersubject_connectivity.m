clear;clc;
projectPath = '/Users/lim_clmn/Desktop/GL/timeseries';
outputdir = fullfile(projectPath, 'intersub_conn');
fbPath = '/Users/lim_clmn/Desktop/GL/GL_beta';
perRun = 144;
nruns = 4;

%getting roi which we'll use
cd /Users/lim_clmn/Desktop/GL/roi
[roilist roinames c] = xlsread('fan280_fullname_GL.xlsx', 1);%dlmread('rois_Yoo.txt');

subjects = {'GL03', 'GL04', 'GL05', 'GL06', 'GL07', 'GL08', 'GL09', 'GL10', 'GL11', 'GL12', 'GL14', 'GL15', 'GL16', ...
    'GL17', 'GL18', 'GL19', 'GL20', 'GL21', 'GL22', 'GL24', 'GL25', 'GL26', 'GL27', 'GL29'};

ct = 0;

for i = 1:length(subjects)
    for ii = 1:length(subjects)
        
        if i < ii
            ct = ct+1;
            cd(projectPath)
            
            tempfilename1 = ['timeseries.' subjects{i} '.2D'];
            tempfilename2 = ['timeseries.' subjects{ii} '.2D'];
            tempfile3_1 = load(tempfilename1);
            tempfile3_2 = load(tempfilename2);
            
            % added 2019-01-28 by YJS: to follow the z-scoring procedure of
            % original manucript (Simony et. al, 2016)
            
            tempfile3_1 = zscore(tempfile3_1);
            tempfile3_2 = zscore(tempfile3_2);
            
 
            % assign whether high or low run
            
            cd(fbPath)
            load([subjects{i} '_fb'])
            fbid2_1 = fbid2;
            clear fbid2;
            load([subjects{ii} '_fb'])
            fbid2_2 = fbid2;
            clear fbid2;
            fbid2_1 = fbid2_1(1:perRun*4);
            fbid2_2 = fbid2_2(1:perRun*4);
            % repeat for each TR
            for jj = 1:length(fbid2_1)
                fbid3_1(2*jj) = fbid2_1(jj);
                fbid3_1(2*jj-1) = fbid2_1(jj);
                fbid3_2(2*jj) = fbid2_2(jj);
                fbid3_2(2*jj-1) = fbid2_2(jj);
            end
            
            % sanity check
            if size(tempfile3_1,1) == length(fbid3_1) && size(tempfile3_2,1) == length(fbid3_2)
            else
                disp(['Length of behavioral data and beta series does not match for ' subjects{i} '.'])
                break
            end
            
            
            % rearrange (added 2019-01-22)
            
            tempfile4_1 = [roilist tempfile3_1'];
            tempfile3_1 = sortrows(tempfile4_1,1);
            tempfile3_1 = tempfile3_1(:,2:end);
            tempfile3_1 = tempfile3_1';
            
            tempfile4_2 = [roilist tempfile3_2'];
            tempfile3_2 = sortrows(tempfile4_2,1);
            tempfile3_2 = tempfile3_2(:,2:end);
            tempfile3_2 = tempfile3_2';
            
            roinames = sortrows(c, 2);
            roinames = roinames(:,1);
            
            idx_fb1 = find(fbid3_1==1); idx_nf1 = find(fbid3_1==0);
            series_fb1 = tempfile3_1(idx_fb1,:); series_nf1 = tempfile3_1(idx_nf1,:);
            idx_fb2 = find(fbid3_2==1); idx_nf2 = find(fbid3_2==0);
            series_fb2 = tempfile3_2(idx_fb2,:); series_nf2 = tempfile3_2(idx_nf2,:);
            
            
            % sanity check that high + low = 288
            if length(series_fb1) + length(series_nf1) == perRun*nruns*2 && length(series_fb2) + length(series_nf2) == perRun*nruns*2
            else
                disp(['Length of high and low trials do not add up to 288 for ' subjects{i} '.'])
                break
            end
            
            % now getting r-values of timeseries
            temp_corrMat_fb = corr(series_fb1, series_fb2);
            temp_corrMat_nf = corr(series_nf1, series_nf2);
            
            
            rt=1;
            for r = 1:length(roinames)
                for rr = 1:length(roinames)
                    
                    if r<=rr
                        temp_corrVec_fb{1,rt} = roinames{rr,1};
                        temp_corrVec_fb{2,rt} = roinames{r,1};
                        temp_corrVec_fb{3,rt} = temp_corrMat_fb(r,rr);
                        ttest_fb(i,rt) = temp_corrMat_fb(r,rr);
                        
                        
                        temp_corrVec_nf{1,rt} = roinames{rr,1};
                        temp_corrVec_nf{2,rt} = roinames{r,1};
                        temp_corrVec_nf{3,rt} = temp_corrMat_nf(r,rr);
                        ttest_nf(i,rt) = temp_corrMat_nf(r,rr);
                        rt = rt+1;
                    end
                    
                end
                
            end
            
            %             cd /Users/lim_clmn/Desktop/GL/timeseries/intersub_conn
            %             eval(['save CorrMatrix_' subjects{i} '_fb temp_corrMat_fb;'])
            %             eval(['save CorrMatrix_' subjects{i} '_nf temp_corrMat_nf;'])
            %             eval(['save CorrVec_' subjects{i} '_fb temp_corrVec_fb;'])
            %             eval(['save CorrVec_' subjects{i} '_nf temp_corrVec_nf;'])
            temp_mat_fb(:,:,ct) = temp_corrMat_fb;
            temp_mat_nf(:,:,ct) = temp_corrMat_nf;
            

            ct
        else
        end
    end
end
% for performing t-test of edges (corr. coefficients)
cd(outputdir)
avgMat_fb = mean(temp_mat_fb,3);
avgMat_nf = mean(temp_mat_nf,3);

figure;subplot(121);imagesc(avgMat_fb, [-0.1 0.1]);subplot(122);imagesc(avgMat_nf, [-0.1 0.1]);colormap(jet);
save temp_mat_fb temp_mat_fb
save temp_mat_nf temp_mat_nf
save avgMat_fb avgMat_fb
save avgMat_nf avgMat_nf
save temp_corrVec_fb temp_corrVec_fb
save temp_corrVec_nf temp_corrVec_nf
%% Averaging
allMat_fb = zeros(277,277);
allMat_nf = zeros(277,277);
for i = 1:length(subjects)
    
    
    
    
    cd /Users/lim_clmn/Desktop/GL/timeseries/result
    eval(['load CorrMatrix_' subjects{i} '_fb;'])
    eval(['load CorrMatrix_' subjects{i} '_nf;'])
    %     eval(['load CorrVec_' subjects{i} '_high temp_corrVec_high;'])
    %     eval(['load CorrVec_' subjects{i} '_low temp_corrVec_low;'])
    %     temp_corrMat_fb = atanh(temp_corrMat_fb);
    %     temp_corrMat_nf = atanh(temp_corrMat_nf);
    allMat_fb = allMat_fb+temp_corrMat_fb1;
    allMat_nf = allMat_nf+temp_corrMat_nf1;
    
    avgMat_fb = allMat_fb/length(subjects);
    avgMat_nf = allMat_nf/length(subjects);
end
save avgMat_fb avgMat_fb
save avgMat_nf avgMat_nf
figure;subplot(121);imagesc(avgMat_fb);subplot(122);imagesc(avgMat_nf);colormap(jet);
diff = avgMat_fb - avgMat_nf;
figure;subplot(121);imagesc(diff);

%% visual inspection of each participant
cd /Users/lim_clmn/Desktop/GL/timeseries/result
subjects = {'GL03', 'GL04', 'GL05', 'GL06', 'GL07', 'GL08', 'GL09', 'GL10', 'GL11', 'GL12', 'GL14', 'GL15', 'GL16', ...
    'GL17', 'GL18', 'GL19', 'GL20', 'GL21', 'GL22', 'GL24', 'GL25', 'GL26', 'GL27', 'GL29'};
for isub = 1:length(subjects)
    load(['CorrMatrix_' subjects{isub} '_fb'])
    load(['CorrMatrix_' subjects{isub} '_nf'])
    %     temp_corrMat_fb = atanh(temp_corrMat_fb);
    %     temp_corrMat_nf = atanh(temp_corrMat_nf);
    figure;subplot(121);imagesc(temp_corrMat_fb1);subplot(122);imagesc(temp_corrMat_nf1);colormap(jet);
end

%% ttest
clear;clc;
subjects = {'GL03', 'GL04', 'GL05', 'GL06', 'GL07', 'GL08', 'GL09', 'GL10', 'GL11', 'GL12', 'GL14', 'GL15', 'GL16', ...
    'GL17', 'GL18', 'GL19', 'GL20', 'GL21', 'GL22', 'GL24', 'GL25', 'GL26', 'GL27', 'GL29'};
% cd /Users/lim_clmn/Desktop/GL/roi
% [roilist roinames c] = xlsread('fan280_fullname_GL.xlsx', 1);%dlmread('rois_Yoo.txt');

cd /Users/lim_clmn/Desktop/GL/timeseries/intersub_conn
load('roinames')
alpha = 0.05;

load temp_mat_fb
load temp_mat_nf
load temp_corrVec_fb

a = alpha/length(temp_corrVec_fb);
for isub = 1:276
%     load(['CorrMatrix_' subjects{isub} '_fb'])
%     load(['CorrMatrix_' subjects{isub} '_nf'])
%     load(['CorrVec_' subjects{isub} '_fb'])
%     load(['CorrVec_' subjects{isub} '_nf'])
    temp_fb1 = temp_mat_fb(:,:,isub);
    temp_nf1 = temp_mat_nf(:,:,isub);
    temp_fb = atanh(temp_fb1);
    temp_nf = atanh(temp_nf1);
    
    rt=1;
    for r = 1:length(roinames)
        for rr = 1:length(roinames)
            if r<rr
                ttest_fb(isub,rt) = temp_fb(r,rr);
                ttest_nf(isub,rt) = temp_nf(r,rr);
                rt = rt+1;
            end
        end
    end
end
tt =0;
for j = 1:length(ttest_fb)
    [h, p, ci, stats] = ttest(ttest_fb(:,j), ttest_nf(:,j));
    if p<a
        tt = tt+1;
        sig_features{tt,1} = temp_corrVec_fb{1,j};
        sig_features{tt,2} = temp_corrVec_fb{2,j};
        sig_features{tt,3} = mean(ttest_fb(:,j));
        sig_features{tt,4} = mean(ttest_nf(:,j));
        sig_features{tt,5} = p;
        disp([temp_corrVec_fb{1,j} ' and ' temp_corrVec_fb{2,j} ' is significant: p = ' num2str(p) ])
    end

end

save sig_features sig_features




