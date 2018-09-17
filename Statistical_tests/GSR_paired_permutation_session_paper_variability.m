function GSR_paired_permutation_session_paper_variability(SPM_dir,Work_dir)

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

%Basic first
tmp1=0;
procedure='Basic';
for number_dataset=1:4
    
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    
    for subject=1:number_subject
        
        if strcmp(dataset,'DatasetGordon')&&(subject==3||subject==8||subject==9)
            continue;
        else
            if strcmp(procedure,'Basic')
                cd([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_summary/DCM/' procedure '/' name_ROI_def '/Full_model/'])
                
                for network_number=1:length(ntwrk_name)
                    
                    disp(['Dataset: ' dataset 'subject: ' num2str(subject)]);
                    try
                        tmp1=tmp1+1;
                        load(['GCM_' ntwrk_name{network_number} '_full_estim']);
                        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                        
                        flag=zeros(1,length(GCM));
                        for diagn=1:length(GCM)
                            if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                                flag(diagn)=1;
                            end
                        end
                        
                        if strcmp(procedure,'Basic')
                            procedure_comp='GSR';
                            load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure_comp '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                        end
                        
                        if strcmp(procedure,'GSR')
                            procedure_comp='Basic';
                            load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure_comp '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                        end
                        
                        flag_comp=zeros(1,length(GCM));
                        for diagn=1:length(GCM)
                            if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                                flag_comp(diagn)=1;
                            end
                        end
                        
                        GCM(find(flag==1|flag_comp==1))=[];
        
                        GCM_group_basic{tmp1}=GCM;
                    catch
                        continue
                    end
                end
            end
        end
    end
end

%GSR
tmp1=0;
procedure='GSR';
for number_dataset=1:4
    
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    
    for subject=1:number_subject
        
        if strcmp(dataset,'DatasetGordon')&&(subject==3||subject==8||subject==9)
            continue;
        else
            if strcmp(procedure,'GSR')
                cd([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_summary/DCM/' procedure '/' name_ROI_def '/Full_model/'])
                
                for network_number=1:length(ntwrk_name)
                    
                    disp(['Dataset: ' dataset 'subject: ' num2str(subject)]);
                    try
                        tmp1=tmp1+1;
                        load(['GCM_' ntwrk_name{network_number} '_full_estim']);
                        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                        
                        flag=zeros(1,length(GCM));
                        for diagn=1:length(GCM)
                            if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                                flag(diagn)=1;
                            end
                        end
                        
                        if strcmp(procedure,'Basic')
                            procedure_comp='GSR';
                            load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure_comp '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                        end
                        
                        if strcmp(procedure,'GSR')
                            procedure_comp='Basic';
                            load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure_comp '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                        end
                        
                        flag_comp=zeros(1,length(GCM));
                        for diagn=1:length(GCM)
                            if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                                flag_comp(diagn)=1;
                            end
                        end
                        
                        
                        GCM(find(flag==1|flag_comp==1))=[];
        
                        GCM_group_GSR{tmp1}=GCM;
                    catch
                        continue
                    end
                end
            end
        end
    end
end

%Basic
for subject=1:17
   for session=1:length(GCM_group_GSR{subject})
       A_matrix_group{subject}(:,session)=GCM_group_basic{subject}{session}.Ep.A(:);
   end
end

%GSR
for subject=1:17
   for session=1:length(GCM_group_GSR{subject})
       A_matrix_group_GSR{subject}(:,session)=GCM_group_GSR{subject}{session}.Ep.A(:);
   end
end

%Test Stat
for subject=1:17
    B=nchoosek(1:length(GCM_group_basic{subject}),2);
    for combi=1:length(B)
        corr_test_1{subject}(combi)=corr(GCM_group_basic{subject}{B(combi,1)}.Ep.A(:),GCM_group_basic{subject}{B(combi,2)}.Ep.A(:));
    end
    mean_corr_test_1(subject)=mean(corr_test_1{subject});
end

for subject=1:17
    B=nchoosek(1:length(GCM_group_basic{subject}),2);
    for combi=1:length(B)
        corr_test_2{subject}(combi)=corr(GCM_group_GSR{subject}{B(combi,1)}.Ep.A(:),GCM_group_GSR{subject}{B(combi,2)}.Ep.A(:));
    end
    mean_corr_test_2(subject)=mean(corr_test_2{subject});
end

%%%%%%%%%%%%%%%%%%%%%
%Paired permutations
%%%%%%%%%%%%%%%%%%%%%
N_permutations=10000;
for subject=1:17
        empirical_values_session{subject}=parfor_GSR_paired_permutation_session_paper_variability(subject,N_permutations,GCM_group_basic,A_matrix_group,A_matrix_group_GSR);
end

for subject=1:17
    if mean_corr_test_2(subject)-mean_corr_test_1(subject)>0
        p_val_session(subject)=(sum(empirical_values_session{subject}>=(mean_corr_test_2(subject)-mean_corr_test_1(subject)))+sum(empirical_values_session{subject}<=(mean_corr_test_1(subject)-mean_corr_test_2(subject))))/length(empirical_values_session{subject});
    else
        p_val_session(subject)=(sum(empirical_values_session{subject}<=(mean_corr_test_2(subject)-mean_corr_test_1(subject)))+sum(empirical_values_session{subject}>=(mean_corr_test_1(subject)-mean_corr_test_2(subject))))/length(empirical_values_session{subject});
    end
end

%p-value cannot be zero, so equal to 1 over number of permutations
for subject=1:17
    if p_val_session(subject)==0
        p_val_session_adapt(subject)=1/N_permutations;
    else
        p_val_session_adapt(subject)=p_val_session(subject);
    end
end

[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_val_session,0.05,'pdep','yes');
[h_adapt, crit_p_adapt, adj_ci_cvrg_adapt, adj_p_adapt]=fdr_bh(p_val_session_adapt,0.05,'pdep','yes');
disp(round(adj_p,3));

empirical_values_session_GSR=empirical_values_session;
p_val_session_adapt_GSR=p_val_session_adapt;
p_val_session_GSR=p_val_session;
adj_p_adapt_GSR=adj_p_adapt;
adj_p_GSR=adj_p;
mean_corr_test_2_GSR=mean_corr_test_2;
mean_corr_test_1_GSR=mean_corr_test_1;

mkdir([Work_dir '/Results_paper_variability/DCM/Basic/Smith/Full_model/PEB_group/']);
cd([Work_dir '/Results_paper_variability/DCM/Basic/Smith/Full_model/PEB_group/']);
save('Emp_val_session_GSR.mat','empirical_values_session_GSR');
save('p_val_session_adapt_GSR.mat','p_val_session_adapt_GSR');
save('p_val_session_GSR.mat','p_val_session_GSR');
save('adj_p_adapt_GSR.mat','adj_p_adapt_GSR');
save('adj_p_GSR.mat','adj_p_GSR');
save('mean_corr_session_GSR.mat','mean_corr_test_1_GSR','mean_corr_test_2_GSR');

end