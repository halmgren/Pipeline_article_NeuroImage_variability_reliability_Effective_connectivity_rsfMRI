function permutation_test_difference_sessions_paired_paper_variability(SPM_dir,Work_dir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Only extrinsic left-right elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(1532658);
N_permut=5000;

name_ROI_def='Smith';

%Create permutations
% for subject=1:17

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
                        %
                        %                             if strcmp(procedure,'Basic')
                        %                                 procedure_comp='GSR';
                        %                                 load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure_comp '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                        %                             end
                        %
                        %                             if strcmp(procedure,'GSR')
                        %                                 procedure_comp='Basic';
                        %                                 load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure_comp '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                        %                             end
                        %
                        %                             flag_comp=zeros(1,length(GCM));
                        %                             for diagn=1:length(GCM)
                        %                                 if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                        %                                     flag_comp(diagn)=1;
                        %                                 end
                        %                             end
                        
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

load([Work_dir '/Results_paper_variability/DCM/Basic/Smith/Full_model/PEB_group/Group_array_results_DMN.mat']);

count_left=0;
count_right=0;
No_asymmetry=[];

for subject=1:length(PEB_group)
    
    if posterior_probability_group(subject)>0.95
        if mean_diff_group(subject)>0
            count_left=count_left+1;
        end
    elseif posterior_probability_group(subject)<0.05
        if mean_diff_group(subject)<0
            count_right=count_right+1;
        end
    elseif posterior_probability_group(subject)>0.05&&posterior_probability_group(subject)<0.95
        No_asymmetry=[No_asymmetry subject];
    end
    
end

GCM_group_basic(No_asymmetry)=[];

%Basic
for subject=1:length(GCM_group_basic)
    for session=1:length(GCM_group_basic{subject})
        A_matrix_group{subject}(:,session)=GCM_group_basic{subject}{session}.Ep.A(:);
    end
end
% end

%%%%%%%%%%%%%%%%
%CHANGE SUBJECT
%%%%%%%%%%%%%%%%

for subject=1:length(GCM_group_basic)
    disp(['subject: ' num2str(subject) '; store permutations']);
    tic
    A=nan(N_permut,size(A_matrix_group{subject},2),6);
    for permut_number=1:N_permut
        for session=1:size(A_matrix_group{subject},2)
            A(permut_number,session,[1 4])=randsample([5 13],2);
            A(permut_number,session,[2 6])=randsample([7 15],2);
            A(permut_number,session,[3 5])=randsample([8 14],2);
        end
    end
    toc
    %     load([Work_dir '/Results_paper_variability/DCM/Basic/Smith/Full_model/PEB_group/Group_array_results_DMN.mat']);
    
    %make matrix of EP's
    %     EP=nan(17,16);
    %     for session=1:length(A_matrix_group{subject})
    %
    %         EP(session,:)=full(PEB_group{session}.Ep);
    %
    %     end
    
    %Make matrix to store permutations
    EP_permut=nan(N_permut,size(A_matrix_group{subject},2),16);
    for permut_number=1:N_permut
        EP_permut(permut_number,:,:)=A_matrix_group{subject}';
    end
    %
    %permute the labels
    EP_permut_tmp=EP_permut(:,:,[5 7:8 13:15]); %save extrinsic connections temporarily
    
    for permut_number=1:N_permut
        for session=1:size(A_matrix_group{subject},2)
            EP_permut(permut_number,session,[5 7:8 13:15])=squeeze(EP_permut(permut_number,session,squeeze(A(permut_number,session,:))));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %Normal delateralization
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    EP=A_matrix_group{subject}';
    EP_delat=A_matrix_group{subject}';
    
    for session=1:size(A_matrix_group{subject},2)
        if mean(A_matrix_group{subject}([13:15],session))<mean(A_matrix_group{subject}([5 7:8],session));
            EP_tmp=reshape(A_matrix_group{subject}(:,session)',4,4);
            EP_tmp([2 4],:)=EP_tmp([4 2],:);
            EP_tmp(:,[2 4])=EP_tmp(:,[4 2]);
            EP_delat(session,:)=spm_vec(EP_tmp);
            clear EP_tmp;
        end
    end
    
    B=nchoosek(1:size(A_matrix_group{subject},2),2);
    for combi=1:length(B)
        tmp3=corrcoef(EP_delat(B(combi,1),:),EP_delat(B(combi,2),:));
        BS_correlation_delat_1(combi)=tmp3(1,2);
    end
    
    B=nchoosek(1:size(A_matrix_group{subject},2),2);
    for combi=1:length(B)
        tmp3=corrcoef(EP(B(combi,1),:),EP(B(combi,2),:));
        BS_correlation_1(combi)=tmp3(1,2);
    end
    
    %average de/increase after delateralization
    Effect_lat_basic(subject)=mean(BS_correlation_delat_1,2)-mean(BS_correlation_1,2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Permuted delateralization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['subject: ' num2str(subject) '; permute labels']);
    tic
    EP_permut_delat=EP_permut;
    for permut_number=1:N_permut
        for session=1:size(A_matrix_group{subject},2)
            
            if mean(EP_permut(permut_number,session,[13:15]))<mean(EP_permut(permut_number,session,[5 7:8]));
                EP_permut_tmp2=reshape(EP_permut(permut_number,session,:),4,4);
                EP_permut_tmp2([2 4],:)=EP_permut_tmp2([4 2],:);
                EP_permut_tmp2(:,[2 4])=EP_permut_tmp2(:,[4 2]);
                EP_permut_delat(permut_number,session,:)=spm_vec(EP_permut_tmp2);
                clear EP_permut_tmp2;
            end
            
        end
    end
    toc
    
    disp(['subject: ' num2str(subject) '; calc delateralized correlations']);
    
    tic
    %delat
    B=nchoosek(1:size(A_matrix_group{subject},2),2);
    for permut_number=1:N_permut
        for combi=1:length(B)
            tmp3=corrcoef(EP_permut_delat(permut_number,B(combi,1),:),EP_permut_delat(permut_number,B(combi,2),:));
            BS_correlation_delat_permut(permut_number,combi)=tmp3(1,2);
        end
    end
    toc
    
    disp(['subject: ' num2str(subject) '; calc asymmetric correlations']);
    
    tic
    %orig
    B=nchoosek(1:size(A_matrix_group{subject},2),2);
    for permut_number=1:N_permut
        for combi=1:length(B)
            tmp3=corrcoef(EP_permut(permut_number,B(combi,1),:),EP_permut(permut_number,B(combi,2),:));
            BS_correlation_permut(permut_number,combi)=tmp3(1,2);
        end
    end
    toc
    
    Effect_lat_permut{subject}=mean(BS_correlation_delat_permut,2)-mean(BS_correlation_permut,2);
    
    %Test statistic
    p_value1(subject)=sum(Effect_lat_permut{subject}>=(mean(BS_correlation_delat_1)-mean(BS_correlation_1)))/N_permut;
    
    clearvars -except p_value1 subject N_permut A_matrix_group Effect_lat_basic GCM_group_basic Effect_lat_permut
end

for subject=1:length(GCM_group_basic)
    if p_value1(subject)==0
        p_value1_adapt(subject)=1/N_permut;
    else
        p_value1_adapt(subject)=p_value1(subject);
    end
end

[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_value1,0.05,'pdep','yes');
[h_adapt, crit_p_adapt, adj_ci_cvrg_adapt, adj_p_adapt]=fdr_bh(p_value1_adapt,0.05,'pdep','yes');

mkdir([Work_dir '/Results_paper_variability/DCM/Basic/Smith/Full_model/PEB_group']);
cd([Work_dir '/Results_paper_variability/DCM/Basic/Smith/Full_model/PEB_group']);

save('Effect_size_session_lat.mat','Effect_lat_basic');
save('Empirical_values.mat','Effect_lat_permut');
save('p_value_adapt.mat','p_value1_adapt');
save('p_value.mat','p_value1');
save('adj_p_adapt.mat','adj_p_adapt');
save('adj_p.mat','adj_p');
end