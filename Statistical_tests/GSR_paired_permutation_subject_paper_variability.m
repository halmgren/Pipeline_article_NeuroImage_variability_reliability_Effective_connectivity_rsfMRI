function GSR_paired_permutation_subject_paper_variability(SPM_dir,Work_dir)

name_ROI_def='Smith';

[ROI_list]=Define_ROIs_paper_variability(name_ROI_def);

tmp=0;
for VOI_number=1:size(ROI_list,1)
    ntwrk=ROI_list{VOI_number,1}(1:3);
    
    if VOI_number>1 && strcmp(ROI_list{VOI_number,1}(1:3),ROI_list{VOI_number-1,1}(1:3))
        ntwrk_size(tmp)=ntwrk_size(tmp)+1;
        continue
        
    else
        tmp=tmp+1;
        ntwrk_size(tmp)=1;
        ntwrk_name{tmp}=ROI_list{VOI_number,1}(1:3);
    end
end

network_number=1;

%Before GSR
procedure='Basic';
clear A_matrix;

tmp=0;
for number_dataset=1:4
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    for subject=1:number_subject
        
        clear PEB;
        tmp=tmp+1;
        
        try
            load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_A_mean_comp_GSR_' ntwrk_name{network_number} '.mat']);
%             load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/Lateralization_index_' ntwrk_name{network_number} '.mat']);
        catch
            continue;
        end
        %We don't include participants with less then 8
        %useful sessions.
        if length(PEB.Snames)<8
            tmp=tmp-1;
            continue
%         elseif posterior_probability<0.95&&posterior_probability>0.05
%             tmp=tmp-1;
%             continue
        else
            PEB_group_basic{tmp}=PEB;
%             Lat_group{tmp}=mean_diff;
        end
        
        clear PEB;
        
    end
end

sz=sqrt(length(DCM{1}.Ep));
A_matrix=zeros(16,length(DCM));

for subject=1:length(PEB_group_basic)
    A_matrix(:,subject)=full(PEB_group_basic{subject}.Ep);
end

%After GSR
procedure='GSR';
tmp=0;
for number_dataset=1:4
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    for subject=1:number_subject
        
        clear PEB;
        tmp=tmp+1;
        
        try
            load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_A_mean_comp_GSR_' ntwrk_name{network_number} '.mat']);
%             load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/Lateralization_index_' ntwrk_name{network_number} '.mat']);
        catch
            continue;
        end
        %We don't include participants with less then 8
        %useful sessions.
        if length(PEB.Snames)<8
            tmp=tmp-1;
            continue
%         elseif posterior_probability<0.95&&posterior_probability>0.05
%             tmp=tmp-1;
%             continue
        else
            PEB_group_GSR{tmp}=PEB;
%             Lat_group{tmp}=mean_diff;
        end
        
        clear PEB;
        
    end
end

sz=sqrt(length(DCM{1}.Ep));
A_matrix_GSR=zeros(16,length(DCM));
for subject=1:length(PEB_group_GSR)
    A_matrix_GSR(:,subject)=full(PEB_group_GSR{subject}.Ep);
end

%%%%%%%%%%%%%%%%
%test_statistic
%%%%%%%%%%%%%%%%
B=nchoosek(1:size(A_matrix,2),2);
for combi=1:length(B)
    tmp3=corrcoef(A_matrix(:,B(combi,1)),A_matrix(:,B(combi,2)));
    BS_correlation(combi)=tmp3(1,2);
end

B=nchoosek(1:size(A_matrix,2),2);
for combi=1:length(B)
    tmp3=corrcoef(A_matrix_GSR(:,B(combi,1)),A_matrix_GSR(:,B(combi,2)));
    BS_correlation_GSR(combi)=tmp3(1,2);
end

Test_Stat=(mean(BS_correlation_GSR)-mean(BS_correlation));
N_permutations=50000;
%%%%%%%%%%%%%%%%%%%%%
%Paired permutations
%%%%%%%%%%%%%%%%%%%%%
for permutation=1:N_permutations
    
    disp(num2str(permutation));
    perm=randi([0:1],1,17);
    
    perm_matrix_1(:,find(perm==1))=A_matrix(:,find(perm==1));
    perm_matrix_1(:,find(perm==0))=A_matrix_GSR(:,find(perm==0));
    
    perm_matrix_2(:,find(perm==0))=A_matrix(:,find(perm==0));
    perm_matrix_2(:,find(perm==1))=A_matrix_GSR(:,find(perm==1));
    
    B=nchoosek(1:size(A_matrix,2),2);
    for combi=1:length(B)
        tmp3_1=corrcoef(perm_matrix_1(:,B(combi,1)),perm_matrix_1(:,B(combi,2)));
        BS_correlation_1(combi)=tmp3_1(1,2);
    end
    
    for combi=1:length(B)
        tmp3_2=corrcoef(perm_matrix_2(:,B(combi,1)),perm_matrix_2(:,B(combi,2)));
        BS_correlation_2(combi)=tmp3_2(1,2);
    end
    
    empirical_values_prd(permutation)=mean(BS_correlation_2)-mean(BS_correlation_1);
end

p_val_prd=(sum(abs(empirical_values_prd)>=abs(Test_Stat)))/length(empirical_values_prd);

mkdir([Work_dir '/Results_paper_variability/DCM/Basic/Smith/Full_model/PEB_group']);
cd([Work_dir '/Results_paper_variability/DCM/Basic/Smith/Full_model/PEB_group']);

save(['GSR_results_subject.mat'],'empirical_values_prd','Test_Stat','p_val_prd');

end