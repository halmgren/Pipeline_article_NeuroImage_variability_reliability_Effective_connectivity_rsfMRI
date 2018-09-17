function [empirical_values_session_subj]=parfor_ROI_paired_permutation_session_paper_variability(subject,N_permutations,GCM_group_basic,A_matrix_c1,A_matrix_c2,comp)

%Compare 4mm to 8mm
if comp==1
    for permutation=1:N_permutations
        disp(permutation);
        perm=randi([0:1],1,length(GCM_group_basic{subject}));
        
        perm_matrix_session_1(:,find(perm==1))=A_matrix_c1{subject}(:,find(perm==1));
        perm_matrix_session_1(:,find(perm==0))=A_matrix_c2{subject}(:,find(perm==0));
        
        perm_matrix_session_2(:,find(perm==0))=A_matrix_c1{subject}(:,find(perm==0));
        perm_matrix_session_2(:,find(perm==1))=A_matrix_c2{subject}(:,find(perm==1));
        
        B=nchoosek(1:size(A_matrix_c1{subject},2),2);
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
    
    %Compare 8mm to 12mm
elseif comp==2
    for permutation=1:N_permutations
        disp(permutation);
        perm=randi([0:1],1,length(GCM_group_basic{subject}));
        
        perm_matrix_session_1(:,find(perm==1))=A_matrix_c1{subject}(:,find(perm==1));
        perm_matrix_session_1(:,find(perm==0))=A_matrix_c2{subject}(:,find(perm==0));
        
        perm_matrix_session_2(:,find(perm==0))=A_matrix_c1{subject}(:,find(perm==0));
        perm_matrix_session_2(:,find(perm==1))=A_matrix_c2{subject}(:,find(perm==1));
        
        B=nchoosek(1:size(A_matrix_c1{subject},2),2);
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
    
    %Compare 12mm to 16mm
elseif comp==3
    for permutation=1:N_permutations
        disp(permutation);
        perm=randi([0:1],1,length(GCM_group_basic{subject}));
        
        perm_matrix_session_1(:,find(perm==1))=A_matrix_c1{subject}(:,find(perm==1));
        perm_matrix_session_1(:,find(perm==0))=A_matrix_c2{subject}(:,find(perm==0));
        
        perm_matrix_session_2(:,find(perm==0))=A_matrix_c1{subject}(:,find(perm==0));
        perm_matrix_session_2(:,find(perm==1))=A_matrix_c2{subject}(:,find(perm==1));
        
        B=nchoosek(1:size(A_matrix_c1{subject},2),2);
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
end