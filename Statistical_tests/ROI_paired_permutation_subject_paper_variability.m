function ROI_paired_permutation_subject_paper_variability(SPM_dir,Work_dir)

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

procedure='ROI_Size';
clear A_matrix;

%4mm
ROI_Size=1;
tmp=0;
for number_dataset=1:4
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    for subject=1:number_subject
        
        clear PEB;
        tmp=tmp+1;
        
        try
            load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/' num2str(4*ROI_Size) '/PEB_A_mean_comp_ROIsize_' ntwrk_name{network_number} '.mat']);
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
            PEB_group_1{tmp}=PEB;
%             Lat_group{tmp}=mean_diff;
        end
        
        clear PEB;
        
    end
end

%8mm
ROI_Size=2;
tmp=0;
for number_dataset=1:4
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    for subject=1:number_subject
        
        clear PEB;
        tmp=tmp+1;
        
        try
            load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/' num2str(4*ROI_Size) '/PEB_A_mean_comp_ROIsize_' ntwrk_name{network_number} '.mat']);
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
            PEB_group_2{tmp}=PEB;
%             Lat_group{tmp}=mean_diff;
        end
        
        clear PEB;
        
    end
end

%12mm
ROI_Size=3;
tmp=0;
for number_dataset=1:4
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    for subject=1:number_subject
        
        clear PEB;
        tmp=tmp+1;
        
        try
            load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/' num2str(4*ROI_Size) '/PEB_A_mean_comp_ROIsize_' ntwrk_name{network_number} '.mat']);
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
            PEB_group_3{tmp}=PEB;
%             Lat_group{tmp}=mean_diff;
        end
        
        clear PEB;
        
    end
end

%16mm
ROI_Size=4;
tmp=0;
for number_dataset=1:4
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    for subject=1:number_subject
        
        clear PEB;
        tmp=tmp+1;
        
        try
            load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/' num2str(4*ROI_Size) '/PEB_A_mean_comp_ROIsize_' ntwrk_name{network_number} '.mat']);
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
            PEB_group_4{tmp}=PEB;
%             Lat_group{tmp}=mean_diff;
        end
        
        clear PEB;
        
    end
end


A_matrix_1=zeros(16,length(DCM));
for subject=1:length(PEB_group_1)
    A_matrix_1(:,subject)=full(PEB_group_1{subject}.Ep);
end


A_matrix_2=zeros(16,length(DCM));
for subject=1:length(PEB_group_2)
    A_matrix_2(:,subject)=full(PEB_group_2{subject}.Ep);
end


A_matrix_3=zeros(16,length(DCM));
for subject=1:length(PEB_group_3)
    A_matrix_3(:,subject)=full(PEB_group_3{subject}.Ep);
end


A_matrix_4=zeros(16,length(DCM));
for subject=1:length(PEB_group_4)
    A_matrix_4(:,subject)=full(PEB_group_4{subject}.Ep);
end

%%%%%%%%%%%%%%%%
%test_statistic
%%%%%%%%%%%%%%%%
B=nchoosek(1:size(A_matrix_1,2),2);
for combi=1:length(B)
    tmp3=corrcoef(A_matrix_1(:,B(combi,1)),A_matrix_1(:,B(combi,2)));
    BS_correlation_1(combi)=tmp3(1,2);
end

B=nchoosek(1:size(A_matrix_2,2),2);
for combi=1:length(B)
    tmp3=corrcoef(A_matrix_2(:,B(combi,1)),A_matrix_2(:,B(combi,2)));
    BS_correlation_2(combi)=tmp3(1,2);
end

B=nchoosek(1:size(A_matrix_3,2),2);
for combi=1:length(B)
    tmp3=corrcoef(A_matrix_3(:,B(combi,1)),A_matrix_3(:,B(combi,2)));
    BS_correlation_3(combi)=tmp3(1,2);
end

B=nchoosek(1:size(A_matrix_4,2),2);
for combi=1:length(B)
    tmp3=corrcoef(A_matrix_4(:,B(combi,1)),A_matrix_4(:,B(combi,2)));
    BS_correlation_4(combi)=tmp3(1,2);
end

Test_Stat_c1=(mean(BS_correlation_1)-mean(BS_correlation_2));
Test_Stat_c2=(mean(BS_correlation_2)-mean(BS_correlation_3));
Test_Stat_c3=(mean(BS_correlation_3)-mean(BS_correlation_4));

%%%%%%%%%%%%%%%%%%%%%
%Paired permutations
%%%%%%%%%%%%%%%%%%%%%
N_permutations=50000;

for permutation=1:N_permutations
    
    disp(num2str(permutation));
    perm=randi([0:1],1,17);
    
    perm_matrix_1(:,find(perm==1))=A_matrix_1(:,find(perm==1));
    perm_matrix_1(:,find(perm==0))=A_matrix_2(:,find(perm==0));
    
    perm_matrix_2(:,find(perm==0))=A_matrix_1(:,find(perm==0));
    perm_matrix_2(:,find(perm==1))=A_matrix_2(:,find(perm==1));
    
    B=nchoosek(1:size(A_matrix_1,2),2);
    for combi=1:length(B)
        tmp3_1=corrcoef(perm_matrix_1(:,B(combi,1)),perm_matrix_1(:,B(combi,2)));
        BS_correlation_1(combi)=tmp3_1(1,2);
    end
    
    for combi=1:length(B)
        tmp3_2=corrcoef(perm_matrix_2(:,B(combi,1)),perm_matrix_2(:,B(combi,2)));
        BS_correlation_2(combi)=tmp3_2(1,2);
    end
    
    empirical_values_prd_c1(permutation)=mean(BS_correlation_2)-mean(BS_correlation_1);
end

for permutation=1:N_permutations
    
    disp(num2str(permutation));
    perm=randi([0:1],1,17);
    
    perm_matrix_1(:,find(perm==1))=A_matrix_2(:,find(perm==1));
    perm_matrix_1(:,find(perm==0))=A_matrix_3(:,find(perm==0));
    
    perm_matrix_2(:,find(perm==0))=A_matrix_2(:,find(perm==0));
    perm_matrix_2(:,find(perm==1))=A_matrix_3(:,find(perm==1));
    
    B=nchoosek(1:size(A_matrix_2,2),2);
    for combi=1:length(B)
        tmp3_1=corrcoef(perm_matrix_1(:,B(combi,1)),perm_matrix_1(:,B(combi,2)));
        BS_correlation_1(combi)=tmp3_1(1,2);
    end
    
    for combi=1:length(B)
        tmp3_2=corrcoef(perm_matrix_2(:,B(combi,1)),perm_matrix_2(:,B(combi,2)));
        BS_correlation_2(combi)=tmp3_2(1,2);
    end
    
    empirical_values_prd_c2(permutation)=mean(BS_correlation_2)-mean(BS_correlation_1);
end

for permutation=1:N_permutations
    
    disp(num2str(permutation));
    perm=randi([0:1],1,17);
    
    perm_matrix_1(:,find(perm==1))=A_matrix_3(:,find(perm==1));
    perm_matrix_1(:,find(perm==0))=A_matrix_4(:,find(perm==0));
    
    perm_matrix_2(:,find(perm==0))=A_matrix_3(:,find(perm==0));
    perm_matrix_2(:,find(perm==1))=A_matrix_4(:,find(perm==1));
    
    B=nchoosek(1:size(A_matrix_3,2),2);
    for combi=1:length(B)
        tmp3_1=corrcoef(perm_matrix_1(:,B(combi,1)),perm_matrix_1(:,B(combi,2)));
        BS_correlation_1(combi)=tmp3_1(1,2);
    end
    
    for combi=1:length(B)
        tmp3_2=corrcoef(perm_matrix_2(:,B(combi,1)),perm_matrix_2(:,B(combi,2)));
        BS_correlation_2(combi)=tmp3_2(1,2);
    end
    
    empirical_values_prd_c3(permutation)=mean(BS_correlation_2)-mean(BS_correlation_1);
end

p_val_prd_c1=(sum(abs(empirical_values_prd_c1)>=abs(Test_Stat_c1)))/length(empirical_values_prd_c1);
p_val_prd_c2=(sum(abs(empirical_values_prd_c2)>=abs(Test_Stat_c2)))/length(empirical_values_prd_c2);
p_val_prd_c3=(sum(abs(empirical_values_prd_c3)>=abs(Test_Stat_c3)))/length(empirical_values_prd_c3);

[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh([p_val_prd_c1 p_val_prd_c2 p_val_prd_c3],0.05,'pdep','yes');

mkdir([Work_dir '/Results_paper_variability/DCM/Basic/Smith/Full_model/PEB_group']);
cd([Work_dir '/Results_paper_variability/DCM/Basic/Smith/Full_model/PEB_group']);

save('BS_results.mat','BS_correlation_1','BS_correlation_2','BS_correlation_3','BS_correlation_4','empirical_values_prd_c1','empirical_values_prd_c2','empirical_values_prd_c3','Test_Stat_c1','Test_Stat_c2','Test_Stat_c3','p_val_prd_c1','p_val_prd_c2','p_val_prd_c3','adj_p');

end