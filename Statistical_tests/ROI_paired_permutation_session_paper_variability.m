function ROI_paired_permutation_session_paper_variability(SPM_dir,Work_dir)

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

%Roi size 4mm
tmp1=0;
procedure='ROI_Size';
size_ROI=1;
for number_dataset=1:4
    
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    
    for subject=1:number_subject
        
        if strcmp(dataset,'DatasetGordon')&&(subject==3||subject==8||subject==9)
            continue;
        else
            if strcmp(procedure,'ROI_Size')
                cd([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_summary/DCM/' procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI)])
                
                for network_number=1:length(ntwrk_name)
                    
                    disp(['Dataset: ' dataset 'subject: ' num2str(subject)]);
                    try
                        tmp1=tmp1+1;
                        load(['GCM_' ntwrk_name{network_number} '_full_estim']);
                        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI) '/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                        
                        flag=zeros(1,length(GCM));
                        for diagn=1:length(GCM)
                            if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                                flag(diagn)=1;
                            end
                        end
                        
                        R=1:4;
                        R(size_ROI)=[];
                        
                        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*R(1)) '/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                        
                        flag_comp1=zeros(1,length(GCM));
                        for diagn=1:length(GCM)
                            if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                                flag_comp1(diagn)=1;
                            end
                        end
                        
                        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*R(2)) '/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                        
                        flag_comp2=zeros(1,length(GCM));
                        for diagn=1:length(GCM)
                            if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                                flag_comp2(diagn)=1;
                            end
                        end
                        
                        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*R(3)) '/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                        
                        flag_comp3=zeros(1,length(GCM));
                        for diagn=1:length(GCM)
                            if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                                flag_comp3(diagn)=1;
                            end
                        end
                        
                        
                        GCM(find(flag==1|flag_comp1==1|flag_comp2==1|flag_comp3==1))=[];
                        
                        GCM_group_1{tmp1}=GCM;
                    catch
                        continue
                    end
                end
            end
        end
    end
end

%Roi size 8mm
tmp1=0;
procedure='ROI_Size';
size_ROI=2;
for number_dataset=1:4
    
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    
    for subject=1:number_subject
        
        if strcmp(dataset,'DatasetGordon')&&(subject==3||subject==8||subject==9)
            continue;
        else
            if strcmp(procedure,'ROI_Size')
                cd([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_summary/DCM/' procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI)])
                
                for network_number=1:length(ntwrk_name)
                    
                    disp(['Dataset: ' dataset 'subject: ' num2str(subject)]);
                    try
                        tmp1=tmp1+1;
                        load(['GCM_' ntwrk_name{network_number} '_full_estim']);
                        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI) '/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                        
                        flag=zeros(1,length(GCM));
                        for diagn=1:length(GCM)
                            if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                                flag(diagn)=1;
                            end
                        end
                        
                        R=1:4;
                        R(size_ROI)=[];
                        
                        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*R(1)) '/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                        
                        flag_comp1=zeros(1,length(GCM));
                        for diagn=1:length(GCM)
                            if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                                flag_comp1(diagn)=1;
                            end
                        end
                        
                        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*R(2)) '/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                        
                        flag_comp2=zeros(1,length(GCM));
                        for diagn=1:length(GCM)
                            if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                                flag_comp2(diagn)=1;
                            end
                        end
                        
                        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*R(3)) '/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                        
                        flag_comp3=zeros(1,length(GCM));
                        for diagn=1:length(GCM)
                            if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                                flag_comp3(diagn)=1;
                            end
                        end
                        
                        
                        GCM(find(flag==1|flag_comp1==1|flag_comp2==1|flag_comp3==1))=[];
                        
                        GCM_group_2{tmp1}=GCM;
                    catch
                        continue
                    end
                end
            end
        end
    end
end

%Roi size 12mm
tmp1=0;
procedure='ROI_Size';
size_ROI=3;
for number_dataset=1:4
    
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    
    for subject=1:number_subject
        
        if strcmp(dataset,'DatasetGordon')&&(subject==3||subject==8||subject==9)
            continue;
        else
            if strcmp(procedure,'ROI_Size')
                cd([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_summary/DCM/' procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI)])
                
                for network_number=1:length(ntwrk_name)
                    
                    disp(['Dataset: ' dataset 'subject: ' num2str(subject)]);
                    try
                        tmp1=tmp1+1;
                        load(['GCM_' ntwrk_name{network_number} '_full_estim']);
                        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI) '/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                        
                        flag=zeros(1,length(GCM));
                        for diagn=1:length(GCM)
                            if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                                flag(diagn)=1;
                            end
                        end
                        
                        R=1:4;
                        R(size_ROI)=[];
                        
                        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*R(1)) '/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                        
                        flag_comp1=zeros(1,length(GCM));
                        for diagn=1:length(GCM)
                            if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                                flag_comp1(diagn)=1;
                            end
                        end
                        
                        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*R(2)) '/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                        
                        flag_comp2=zeros(1,length(GCM));
                        for diagn=1:length(GCM)
                            if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                                flag_comp2(diagn)=1;
                            end
                        end
                        
                        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*R(3)) '/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                        
                        flag_comp3=zeros(1,length(GCM));
                        for diagn=1:length(GCM)
                            if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                                flag_comp3(diagn)=1;
                            end
                        end
                        
                        
                        GCM(find(flag==1|flag_comp1==1|flag_comp2==1|flag_comp3==1))=[];
                        
                        GCM_group_3{tmp1}=GCM;
                    catch
                        continue
                    end
                end
            end
        end
    end
end

%Roi size 16mm
tmp1=0;
procedure='ROI_Size';
size_ROI=4;
for number_dataset=1:4
    
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    
    for subject=1:number_subject
        
        if strcmp(dataset,'DatasetGordon')&&(subject==3||subject==8||subject==9)
            continue;
        else
            if strcmp(procedure,'ROI_Size')
                cd([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_summary/DCM/' procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI)])
                
                for network_number=1:length(ntwrk_name)
                    
                    disp(['Dataset: ' dataset 'subject: ' num2str(subject)]);
                    try
                        tmp1=tmp1+1;
                        load(['GCM_' ntwrk_name{network_number} '_full_estim']);
                        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*size_ROI) '/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                        
                        flag=zeros(1,length(GCM));
                        for diagn=1:length(GCM)
                            if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                                flag(diagn)=1;
                            end
                        end
                        
                        R=1:4;
                        R(size_ROI)=[];
                        
                        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*R(1)) '/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                        
                        flag_comp1=zeros(1,length(GCM));
                        for diagn=1:length(GCM)
                            if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                                flag_comp1(diagn)=1;
                            end
                        end
                        
                        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*R(2)) '/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                        
                        flag_comp2=zeros(1,length(GCM));
                        for diagn=1:length(GCM)
                            if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                                flag_comp2(diagn)=1;
                            end
                        end
                        
                        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/' num2str(4*R(3)) '/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                        
                        flag_comp3=zeros(1,length(GCM));
                        for diagn=1:length(GCM)
                            if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                                flag_comp3(diagn)=1;
                            end
                        end
                        
                        
                        GCM(find(flag==1|flag_comp1==1|flag_comp2==1|flag_comp3==1))=[];
                        
                        GCM_group_4{tmp1}=GCM;
                    catch
                        continue
                    end
                end
            end
        end
    end
end

%Reorder matrices
for subject=1:17
    for session=1:length(GCM_group_1{subject})
        A_matrix_1{subject}(:,session)=GCM_group_1{subject}{session}.Ep.A(:);
    end
end

for subject=1:17
    for session=1:length(GCM_group_2{subject})
        A_matrix_2{subject}(:,session)=GCM_group_2{subject}{session}.Ep.A(:);
    end
end

for subject=1:17
    for session=1:length(GCM_group_3{subject})
        A_matrix_3{subject}(:,session)=GCM_group_3{subject}{session}.Ep.A(:);
    end
end

for subject=1:17
    for session=1:length(GCM_group_4{subject})
        A_matrix_4{subject}(:,session)=GCM_group_4{subject}{session}.Ep.A(:);
    end
end

%Test Stat
for subject=1:17
    B=nchoosek(1:length(GCM_group_1{subject}),2);
    for combi=1:length(B)
        corr_test_1{subject}(combi)=corr(GCM_group_1{subject}{B(combi,1)}.Ep.A(:),GCM_group_1{subject}{B(combi,2)}.Ep.A(:));
    end
    mean_corr_test_1(subject)=mean(corr_test_1{subject});
end

for subject=1:17
    B=nchoosek(1:length(GCM_group_2{subject}),2);
    for combi=1:length(B)
        corr_test_2{subject}(combi)=corr(GCM_group_2{subject}{B(combi,1)}.Ep.A(:),GCM_group_2{subject}{B(combi,2)}.Ep.A(:));
    end
    mean_corr_test_2(subject)=mean(corr_test_2{subject});
end

for subject=1:17
    B=nchoosek(1:length(GCM_group_3{subject}),2);
    for combi=1:length(B)
        corr_test_3{subject}(combi)=corr(GCM_group_3{subject}{B(combi,1)}.Ep.A(:),GCM_group_3{subject}{B(combi,2)}.Ep.A(:));
    end
    mean_corr_test_3(subject)=mean(corr_test_3{subject});
end

for subject=1:17
    B=nchoosek(1:length(GCM_group_4{subject}),2);
    for combi=1:length(B)
        corr_test_4{subject}(combi)=corr(GCM_group_4{subject}{B(combi,1)}.Ep.A(:),GCM_group_4{subject}{B(combi,2)}.Ep.A(:));
    end
    mean_corr_test_4(subject)=mean(corr_test_4{subject});
end

%%%%%%%%%%%%%%%%%%%%%
%Paired permutations
%%%%%%%%%%%%%%%%%%%%%
N_permutations=10000;

comp=1;
parfor subject=1:17
    empirical_values_session_comp1{subject}=parfor_ROI_paired_permutation_session_paper_variability(subject,N_permutations,GCM_group_1,A_matrix_1,A_matrix_2,comp);
end

comp=2;
parfor subject=1:17
    empirical_values_session_comp2{subject}=parfor_ROI_paired_permutation_session_paper_variability(subject,N_permutations,GCM_group_2,A_matrix_2,A_matrix_3,comp);
end

comp=3;
parfor subject=1:17
    empirical_values_session_comp3{subject}=parfor_ROI_paired_permutation_session_paper_variability(subject,N_permutations,GCM_group_3,A_matrix_3,A_matrix_4,comp);
end



for subject=1:17
    if mean_corr_test_2(subject)-mean_corr_test_1(subject)>0
        p_val_session_c1(subject)=(sum(empirical_values_session_comp1{subject}>=(mean_corr_test_2(subject)-mean_corr_test_1(subject)))+sum(empirical_values_session_comp1{subject}<=(mean_corr_test_1(subject)-mean_corr_test_2(subject))))/length(empirical_values_session_comp1{subject});
    else
        p_val_session_c1(subject)=(sum(empirical_values_session_comp1{subject}<=(mean_corr_test_2(subject)-mean_corr_test_1(subject)))+sum(empirical_values_session_comp1{subject}>=(mean_corr_test_1(subject)-mean_corr_test_2(subject))))/length(empirical_values_session_comp1{subject});
    end
end

for subject=1:17
    if mean_corr_test_3(subject)-mean_corr_test_2(subject)>0
        p_val_session_c2(subject)=(sum(empirical_values_session_comp2{subject}>=(mean_corr_test_3(subject)-mean_corr_test_2(subject)))+sum(empirical_values_session_comp2{subject}<=(mean_corr_test_2(subject)-mean_corr_test_3(subject))))/length(empirical_values_session_comp2{subject});
    else
        p_val_session_c2(subject)=(sum(empirical_values_session_comp2{subject}<=(mean_corr_test_3(subject)-mean_corr_test_2(subject)))+sum(empirical_values_session_comp2{subject}>=(mean_corr_test_2(subject)-mean_corr_test_3(subject))))/length(empirical_values_session_comp2{subject});
    end
end

for subject=1:17
    if mean_corr_test_4(subject)-mean_corr_test_3(subject)>0
        p_val_session_c3(subject)=(sum(empirical_values_session_comp3{subject}>=(mean_corr_test_4(subject)-mean_corr_test_3(subject)))+sum(empirical_values_session_comp3{subject}<=(mean_corr_test_3(subject)-mean_corr_test_4(subject))))/length(empirical_values_session_comp3{subject});
    else
        p_val_session_c3(subject)=(sum(empirical_values_session_comp3{subject}<=(mean_corr_test_4(subject)-mean_corr_test_3(subject)))+sum(empirical_values_session_comp3{subject}>=(mean_corr_test_3(subject)-mean_corr_test_4(subject))))/length(empirical_values_session_comp3{subject});
    end
end

%p-value cannot be zero, so equal to 1 over number of permutations
for subject=1:17
    if p_val_session_c1(subject)==0
        p_val_session_c1_adapt(subject)=1/N_permutations;
    else
        p_val_session_c1_adapt(subject)=p_val_session_c1(subject);
    end
end

for subject=1:17
    if p_val_session_c2(subject)==0
        p_val_session_c2_adapt(subject)=1/N_permutations;
    else
        p_val_session_c2_adapt(subject)=p_val_session_c2(subject);
    end
end

for subject=1:17
    if p_val_session_c3(subject)==0
        p_val_session_c3_adapt(subject)=1/N_permutations;
    else
        p_val_session_c3_adapt(subject)=p_val_session_c3(subject);
    end
end

[h_c1, crit_p_c1, adj_ci_cvrg_c1, adj_p_c1]=fdr_bh(p_val_session_c1,0.05,'pdep','yes');
[h_c1_adapt, crit_p_c1_adapt, adj_ci_cvrg_c1_adapt, adj_p_c1_adapt]=fdr_bh(p_val_session_c1_adapt,0.05,'pdep','yes');
disp(round(adj_p_c1_adapt,3));

[h_c2, crit_p_c2, adj_ci_cvrg_c2, adj_p_c2]=fdr_bh(p_val_session_c2,0.05,'pdep','yes');
[h_c2_adapt, crit_c2_p_adapt, adj_c2_ci_cvrg_adapt, adj_p_c2_adapt]=fdr_bh(p_val_session_c2_adapt,0.05,'pdep','yes');
disp(round(adj_p_c2_adapt,3));

[h_c3, crit_p_c3, adj_ci_cvrg_c3, adj_p_c3]=fdr_bh(p_val_session_c3,0.05,'pdep','yes');
[h_c3_adapt, crit_p_c3_adapt, adj_ci_cvrg_c3_adapt, adj_p_c3_adapt]=fdr_bh(p_val_session_c3_adapt,0.05,'pdep','yes');
disp(round(adj_p_c3_adapt,3));

[h_all, crit_p_all, adj_ci_cvrg_all, adj_p_all]=fdr_bh([p_val_session_c1 p_val_session_c2 p_val_session_c3],0.05,'pdep','yes');
[h_all_adapt, crit_p_all_adapt, adj_ci_cvrg_all_adapt, adj_p_all_adapt]=fdr_bh([p_val_session_c1 p_val_session_c2 p_val_session_c3],0.05,'pdep','yes');
disp(round(adj_p_all_adapt,3));

mkdir([Work_dir '/Results_paper_variability/DCM/Basic/Smith/Full_model/PEB_group/']);
cd([Work_dir '/Results_paper_variability/DCM/Basic/Smith/Full_model/PEB_group/']);

save('Emp_val_session_ROI.mat','empirical_values_session_comp1','empirical_values_session_comp2','empirical_values_session_comp3');
save('p_val_session_adapt_ROI.mat','p_val_session_c1_adapt','p_val_session_c2_adapt','p_val_session_c3_adapt');
save('p_val_session_ROI.mat','p_val_session_c1','p_val_session_c2','p_val_session_c3');
save('mean_corr_ROI.mat','mean_corr_test_1','mean_corr_test_2','mean_corr_test_3','mean_corr_test_4');
save('adj_p_adapt_ROI.mat','adj_p_c1_adapt','adj_p_c2_adapt','adj_p_c3_adapt','adj_p_all_adapt');
save('adj_p_ROI.mat','adj_p_c1','adj_p_c2','adj_p_c3','adj_p_all');

end
