function BMR_paired_permutation_subject_paper_variability(SPM_dir,Work_dir)

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

%Before BMR
procedure='Basic';
clear A_matrix;
tmp=0;
for number_dataset=1:4
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    for subject=1:number_subject
        
        clear PEB;
        tmp=tmp+1;
        
        try
            load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_A_mean_' ntwrk_name{network_number} '.mat']);
            load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/Lateralization_index_' ntwrk_name{network_number} '.mat']);
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
            PEB_group{tmp}=PEB;
            Lat_group{tmp}=mean_diff;
        end
        
        clear PEB;
        
    end
end

A_matrix=zeros(16,length(DCM));
for subject=1:length(PEB_group)
    A_matrix(:,subject)=full(PEB_group{subject}.Ep);
end

%After BMR
procedure='BMR';
load([Work_dir '/Results_paper_variability/DCM/Basic/' name_ROI_def '/Full_model/PEB_group/PEB_A_mean_group_' ntwrk_name{network_number} '.mat'],'PEB','DCM');
    
A_matrix_BMR=zeros(16,length(DCM));
for subject=1:length(DCM)
    A_matrix_BMR(:,subject)=full(DCM{subject}.Ep);
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
    tmp3=corrcoef(A_matrix_BMR(:,B(combi,1)),A_matrix_BMR(:,B(combi,2)));
    BS_correlation_BMR(combi)=tmp3(1,2);
end

Test_Stat=(mean(BS_correlation_BMR)-mean(BS_correlation));

%%%%%%%%%%%%%%%%%%%%%
%Paired permutations
%%%%%%%%%%%%%%%%%%%%%
for permutation=1:50000
    disp(num2str(permutation));
    perm=randi([0:1],1,17);
    
    perm_matrix_1(:,find(perm==1))=A_matrix(:,find(perm==1));
    perm_matrix_1(:,find(perm==0))=A_matrix_BMR(:,find(perm==0));
    
    perm_matrix_2(:,find(perm==0))=A_matrix(:,find(perm==0));
    perm_matrix_2(:,find(perm==1))=A_matrix_BMR(:,find(perm==1));
    
    B=nchoosek(1:size(A_matrix,2),2);
    for combi=1:length(B)
        tmp3_1=corrcoef(perm_matrix_1(:,B(combi,1)),perm_matrix_1(:,B(combi,2)));
        BS_correlation_1(combi)=tmp3_1(1,2);
    end
    
    for combi=1:length(B)
        tmp3_2=corrcoef(perm_matrix_2(:,B(combi,1)),perm_matrix_2(:,B(combi,2)));
        BS_correlation_2(combi)=tmp3_2(1,2);
    end
    
    empirical_values_subject(permutation)=mean(BS_correlation_2)-mean(BS_correlation_1);
end

p_val_subject=sum(empirical_values_subject>=Test_Stat)/length(empirical_values_subject);
% figure; hist(empirical_values);
p_val_subject_BMR=p_val_subject;
empirical_values_subject_BMR=empirical_values_subject;


cd([Work_dir '/Results_paper_variability/DCM/Basic/Smith/Full_model/PEB_group']);
save('Emp_val_subject_BMR.mat','empirical_values_subject_BMR');
save('p_val_subject_BMR.mat','p_val_subject_BMR');
save('Effect_size_BMR.mat','Test_Stat');
end
