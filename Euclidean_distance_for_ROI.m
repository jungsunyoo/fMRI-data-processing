% Getting Euclidean distance between ROIs 
% Created by Jungsun Yoo
clear;clc;
%cd rootdir
[center, txt, raw] = xlsread('OVRL_Tmap_2',1,'a1:d33');
%subjects = { };
PPI_name = 'MemRew';
RunNumber = 1;
ConName = 'HIT';
% cd(fullfile(rootdir, filesep, PPI_Extracted))
for nsub = 1:length(subjects)
    load(['PPI_Extraction_',PPI_name, '_', subjects{nsub},'_Run_',num2str(RunNumber),'_Con_',ConName,'_Phase2_Results.mat']);
    for nroi = 1:length(center)
        temp_roi(1,1) = center(nroi,1);
        temp_roi(2,1) = center(nroi,2);
        temp_roi(3,1) = center(nroi,3);
        for nvox = 1:length(PP_Interaction_Results_XYZmm_List)
            temp_vox(1,1) = PP_Interaction_Results_XYZmm_List(1,nvox);
            temp_vox(2,1) = PP_Interaction_Results_XYZmm_List(2,nvox);
            temp_vox(3,1) = PP_Interaction_Results_XYZmm_List(3,nvox);
            temp_euclidean(1,nvox) = sqrt((temp_roi(1,1)-temp_vox(1,1))^2 + (temp_roi(2,1)-temp_vox(2,1))^2 +(temp_roi(3,1)-temp_vox(3,1))^2);
        end
        temp_index = find(temp_euclidean == min(temp_euclidean));
        if length(temp_index)>1
            temp_index = temp_index(1);
        end
        euc_sub(nsub, 3*nroi-2) = PP_Interaction_Results_XYZmm_List(1,temp_index);
        euc_sub(nsub, 3*nroi-1) = PP_Interaction_Results_XYZmm_List(2,temp_index);
        euc_sub(nsub, 3*nroi) = PP_Interaction_Results_XYZmm_List(3,temp_index);

        clearvars temp_roi temp_vox temp_euclidean temp_index
    end
        clearvars PP_Interaction_Results_XYZmm_List
end

%% Sorting
for i = 1:length(center)
        final_nodes(i,1) =  euc_sub(1, 3*i-2);
        final_nodes(i,2) =  euc_sub(1, 3*i-1);
        final_nodes(i,3) =  euc_sub(1, 3*i);
end

        
