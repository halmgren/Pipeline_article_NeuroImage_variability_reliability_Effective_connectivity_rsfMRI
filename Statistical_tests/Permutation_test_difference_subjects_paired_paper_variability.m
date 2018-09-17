function Permutation_test_difference_subjects_paired_paper_variability(SPM_dir,Work_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Only extrinsic left-right elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(1532658);
N_permut=5000;

%Create permutations
A=nan(N_permut,17,6);
for permut_number=1:N_permut
    for subject=1:17
        A(permut_number,subject,[1 4])=randsample([5 13],2);
        A(permut_number,subject,[2 6])=randsample([7 15],2);
        A(permut_number,subject,[3 5])=randsample([8 14],2);
    end
end

load([Work_dir '/Results_paper_variability/DCM/Basic/Smith/Full_model/PEB_group/Group_array_results_DMN.mat']);

%make matrix of EP's
EP=nan(17,16);
for subject=1:17
    
    EP(subject,:)=full(PEB_group{subject}.Ep);
    
end

%Make matrix to store permutations
EP_permut=nan(N_permut,17,16);
for permut_number=1:N_permut
    EP_permut(permut_number,:,:)=EP;
end

%permute the labels
% EP_permut_tmp=EP_permut(:,:,[5 7:8 13:15]); %save extrinsic connections temporarily

for permut_number=1:N_permut
    for subject=1:17
        EP_permut(permut_number,subject,[5 7:8 13:15])=squeeze(EP_permut(permut_number,subject,squeeze(A(permut_number,subject,:))));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Normal delateralization
%%%%%%%%%%%%%%%%%%%%%%%%%%
EP_delat=EP;

for subject=1:17
    if mean(EP(subject,[13:15]))<mean(EP(subject,[5 7:8]));
        EP_tmp=reshape(EP(subject,:),4,4);
        EP_tmp([2 4],:)=EP_tmp([4 2],:);
        EP_tmp(:,[2 4])=EP_tmp(:,[4 2]);
        EP_delat(subject,:)=spm_vec(EP_tmp);
        clear EP_tmp;
    end
end

B=nchoosek(1:17,2);
for combi=1:length(B)
    tmp3=corrcoef(EP_delat(B(combi,1),:),EP_delat(B(combi,2),:));
    BS_correlation_delat_1(combi)=tmp3(1,2);
end

B=nchoosek(1:17,2);
for combi=1:length(B)
    tmp3=corrcoef(EP(B(combi,1),:),EP(B(combi,2),:));
    BS_correlation_1(combi)=tmp3(1,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Permuted delateralization
%%%%%%%%%%%%%%%%%%%%%%%%%%%
EP_permut_delat=EP_permut;
for permut_number=1:N_permut
    for subject=1:17
        
        if mean(EP_permut(permut_number,subject,[13:15]))<mean(EP_permut(permut_number,subject,[5 7:8]));
            EP_permut_tmp2=reshape(EP_permut(permut_number,subject,:),4,4);
            EP_permut_tmp2([2 4],:)=EP_permut_tmp2([4 2],:);
            EP_permut_tmp2(:,[2 4])=EP_permut_tmp2(:,[4 2]);
            EP_permut_delat(permut_number,subject,:)=spm_vec(EP_permut_tmp2);
            clear EP_permut_tmp2;
        end
        
    end
end

%delat
B=nchoosek(1:17,2);
for permut_number=1:N_permut
    for combi=1:length(B)
        tmp3=corrcoef(EP_permut_delat(permut_number,B(combi,1),:),EP_permut_delat(permut_number,B(combi,2),:));
        BS_correlation_delat_permut(permut_number,combi)=tmp3(1,2);
    end
end

%orig
B=nchoosek(1:17,2);
for permut_number=1:N_permut
    for combi=1:length(B)
        tmp3=corrcoef(EP_permut(permut_number,B(combi,1),:),EP_permut(permut_number,B(combi,2),:));
        BS_correlation_permut(permut_number,combi)=tmp3(1,2);
    end
end

Effect_lat_permut=mean(BS_correlation_delat_permut,2)-mean(BS_correlation_permut,2);

%Test statistic
p_value1=sum(Effect_lat_permut>=(mean(BS_correlation_delat_1)-mean(BS_correlation_1)))/N_permut;

% mkdir([Work_dir '/results_perm_tests/Asymmetry/subject']);
% cd([Work_dir '/results_perm_tests/Asymmetry/subject']);

% save('Empirical_values.mat','Effect_lat_permut');
save([Work_dir '/Results_paper_variability/DCM/Basic/Smith/Full_model/PEB_group/p_value_BS_reordered.mat'],'p_value1');

end
