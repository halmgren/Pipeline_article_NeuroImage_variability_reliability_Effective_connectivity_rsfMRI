function [empirical_values_session_subj]=parfor_GSR_paired_permutation_session_paper_variability(subject,N_permutations,GCM_group_basic,A_matrix_group,A_matrix_group_GSR)

for permutation=1:N_permutations
    perm=randi([0:1],1,length(GCM_group_basic{subject}));
    
    perm_matrix_session_1(:,find(perm==1))=A_matrix_group{subject}(:,find(perm==1));
    perm_matrix_session_1(:,find(perm==0))=A_matrix_group_GSR{subject}(:,find(perm==0));
    
    perm_matrix_session_2(:,find(perm==0))=A_matrix_group{subject}(:,find(perm==0));
    perm_matrix_session_2(:,find(perm==1))=A_matrix_group_GSR{subject}(:,find(perm==1));
    
    B=nchoosek(1:size(A_matrix_group{subject},2),2);
    for combi=1:length(B)
        tmp3_1=corrcoef(perm_matrix_session_1(:,B(combi,1)),perm_matrix_session_1(:,B(combi,2)));
        BS_correlation_session_1(combi)=tmp3_1(1,2);
    end
    
    for combi=1:length(B)
        tmp3_2=corrcoef(perm_matrix_session_2(:,B(combi,1)),perm_matrix_session_2(:,B(combi,2)));
        BS_correlation_session_2(combi)=tmp3_2(1,2);
    end
    
    empirical_values_session_subj(permutation)=mean(BS_correlation_session_2)-mean(BS_correlation_session_1);
end

end