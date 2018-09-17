function BMR_paired_permutation_session_paper_variability(SPM_dir,Work_dir)
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
                        
                        GCM(find(flag==1))=[];
                        
                        GCM_group_basic{tmp1}=GCM;
                    catch
                        continue
                    end
                end
            end
        end
    end
end

%BMR
tmp1=0;
procedure='Basic';
for number_dataset=1:4
    
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    
    for subject=1:number_subject
        
        if strcmp(dataset,'DatasetGordon')&&(subject==3||subject==8||subject==9)
            continue;
        else
            if strcmp(procedure,'Basic')
                cd([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/'])
                
                for network_number=1:length(ntwrk_name)
                    
                    disp(['Dataset: ' dataset 'subject: ' num2str(subject)]);
                    try
                        tmp1=tmp1+1;
                        load(['PEB_A_mean_' ntwrk_name{network_number} '.mat']);
%                         load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
%                         
%                         flag=zeros(1,length(GCM));
%                         for diagn=1:length(GCM)
%                             if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
%                                 flag(diagn)=1;
%                             end
%                         end
%                         
%                         GCM(find(flag==1))=[];
                        
                        GCM_group_BMR{tmp1}=DCM;
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
   for session=1:length(GCM_group_BMR{subject})
       A_matrix_group{subject}(:,session)=GCM_group_basic{subject}{session}.Ep.A(:);
   end
end

%BMR
for subject=1:17
   for session=1:length(GCM_group_BMR{subject})
       A_matrix_group_BMR{subject}(:,session)=GCM_group_BMR{subject}{session}.Ep.A(:);
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
        corr_test_2{subject}(combi)=corr(GCM_group_BMR{subject}{B(combi,1)}.Ep.A(:),GCM_group_BMR{subject}{B(combi,2)}.Ep.A(:));
    end
    mean_corr_test_2(subject)=mean(corr_test_2{subject});
end

%%%%%%%%%%%%%%%%%%%%%
%Paired permutations
%%%%%%%%%%%%%%%%%%%%%
N_permutations=10000;
for subject=1:17
    disp(['subject: ' num2str(subject) ' is being permuted']);
    empirical_values_session{subject}=parfor_BMR_paired_permutation_session_paper_variability(subject,N_permutations,GCM_group_basic,A_matrix_group,A_matrix_group_BMR);
end

for subject=1:17
    p_val_session(subject)=sum(empirical_values_session{subject}>=(mean_corr_test_2(subject)-mean_corr_test_1(subject)))/length(empirical_values_session{subject});
end

%prepare for FDR
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

mkdir([Work_dir '/Results_paper_variability/DCM/Basic/Smith/Full_model/PEB_group/']);
cd([Work_dir '/Results_paper_variability/DCM/Basic/Smith/Full_model/PEB_group/']);

empirical_values_session_BMR_ses=empirical_values_session;
p_val_session_adapt_BMR_ses=p_val_session_adapt;
adj_p_adapt_BMR_ses=adj_p_adapt;
adj_p_BMR_ses=adj_p;

save('Effect_sizes.mat','mean_corr_test_1','mean_corr_test_2');
save('Emp_val_session_BMR_ses.mat','empirical_values_session_BMR_ses');
save('p_val_session_adapt_BMR_ses.mat','p_val_session_adapt_BMR_ses');
save('p_val_session_BMR_ses.mat','p_val_session');
save('adj_p_adapt_BMR_ses.mat','adj_p_adapt_BMR_ses');
save('adj_p_BMR_ses.mat','adj_p_BMR_ses');


end
