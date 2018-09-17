function Table_lat_inf(SPM_dir,Work_dir)

load([Work_dir '/Results_paper_variability/DCM/Basic/Smith/Full_model/PEB_group/Group_array_results_DMN.mat'],'PEB_group','mean_diff_group','var_of_sum_group','posterior_probability_group');
mean_diff_group_subject_level=mean_diff_group;
clear PEB_group var_of_sum_group posterior_probability_group mean_diff_group;

tmp=0;

for number_dataset=1:4
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    for subject=1:number_subject
        tmp=tmp+1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %These specific subjects were excluded
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(dataset,'DatasetGordon')&&(subject==3||subject==8||subject==9)
            tmp=tmp-1;
            continue
        end
        
        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/Basic/Smith/Full_model/Lateral_index_individ_DMN.mat']); %load lateralization
        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/Basic/Smith/Full_model/QC/Above_treshold_marks_DMN.mat']); %load QC's
        
        flag=zeros(1,length(mean_diff));
        for diagn=1:length(mean_diff)
            if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                flag(diagn)=1;
            end
        end
        
        mean_diff(flag==1)=[];
        posterior_probability(flag==1)=[];
        flag2=zeros(1,length(mean_diff));
        
        for diagn=1:length(mean_diff)
            if posterior_probability(diagn)<0.95&&posterior_probability(diagn)>0.05
                flag2(diagn)=1;
            end
        end
        
        mean_diff_group{tmp}=mean_diff;
        %stab4(tmp)=sum(mean_diff(flag2==0)<0)/length(mean_diff(flag2==0)); 
        
        if mean_diff_group_subject_level(tmp)<0
            stab6(tmp)=sum(mean_diff(flag2==0)<0)/length(mean_diff); %%%%%%%this one was used oin calculations of perc across all subjects
        else
            stab6(tmp)=sum(mean_diff(flag2==0)>0)/length(mean_diff);
        end
        
        %stab5(tmp)=length(mean_diff(flag2==0))/length(mean_diff);
        N_left_sessions(tmp)=sum(mean_diff(flag2==0)>0);
        N_used_ses(tmp)=length(mean_diff(flag2==0));
        N_tot_ses(tmp)=length(mean_diff);
        Non_sign{tmp,:}=flag2;
        if number_dataset==1
            rownam{tmp}=['S' num2str(subject)];
        elseif number_dataset==2
            rownam{tmp}=['S' num2str(subject+8)];
        elseif number_dataset==3
            rownam{tmp}=['S' num2str(subject+9)];
        elseif number_dataset==4
            rownam{tmp}=['S' num2str(subject+10)];
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

%%%%%%%%%%%%%%%%%%%%%%%%
%average same asymmetry
%%%%%%%%%%%%%%%%%%%%%%%%

N_left_sessions(No_asymmetry)=[];
N_used_ses(No_asymmetry)=[];
N_tot_ses(No_asymmetry)=[];
rownam(No_asymmetry)=[];

%stab4(No_asymmetry)=[];
%stab5(No_asymmetry)=[];
stab6(No_asymmetry)=[];

%on average, 71% of sessions...
mean(stab6);

%In total 94% if sessions showed significant lateralization
Prop_lateral=sum(N_used_ses)/sum(N_tot_ses);

N_p_values=myBinomTest(N_left_sessions,N_used_ses,0.5,'two');
disp(round(N_p_values,4));

[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(N_p_values,0.05,'pdep','yes');
disp(round(adj_p,3));


tmp_mean_diff_group_subject_level=mean_diff_group_subject_level;
tmp_mean_diff_group_subject_level(No_asymmetry)=[];

for subj=1:length(N_left_sessions)
    if tmp_mean_diff_group_subject_level(subj)<0load([Work_dir '/Results_paper_variability/DCM/Basic/Smith/Full_model/PEB_group/Group_array_results_DMN.mat'],'PEB_group','mean_diff_group','var_of_sum_group','posterior_probability_group');
mean_diff_group_subject_level=mean_diff_group;
clear PEB_group var_of_sum_group posterior_probability_group mean_diff_group;

tmp=0;

for number_dataset=1:4
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    for subject=1:number_subject
        tmp=tmp+1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %These specific subjects were excluded
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(dataset,'DatasetGordon')&&(subject==3||subject==8||subject==9)
            tmp=tmp-1;
            continue
        end
        
        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/Basic/Smith/Full_model/Lateral_index_individ_DMN.mat']); %load lateralization
        load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/Basic/Smith/Full_model/QC/Above_treshold_marks_DMN.mat']); %load QC's
        
        flag=zeros(1,length(mean_diff));
        for diagn=1:length(mean_diff)
            if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                flag(diagn)=1;
            end
        end
        
        mean_diff(flag==1)=[];
        posterior_probability(flag==1)=[];
        flag2=zeros(1,length(mean_diff));
        
        for diagn=1:length(mean_diff)
            if posterior_probability(diagn)<0.95&&posterior_probability(diagn)>0.05
                flag2(diagn)=1;
            end
        end
        
        mean_diff_group{tmp}=mean_diff;
        %stab4(tmp)=sum(mean_diff(flag2==0)<0)/length(mean_diff(flag2==0)); 
        
        if mean_diff_group_subject_level(tmp)<0
            stab6(tmp)=sum(mean_diff(flag2==0)<0)/length(mean_diff); %%%%%%%this one was used oin calculations of perc across all subjects
        else
            stab6(tmp)=sum(mean_diff(flag2==0)>0)/length(mean_diff);
        end
        
        %stab5(tmp)=length(mean_diff(flag2==0))/length(mean_diff);
        N_left_sessions(tmp)=sum(mean_diff(flag2==0)>0);
        N_used_ses(tmp)=length(mean_diff(flag2==0));
        N_tot_ses(tmp)=length(mean_diff);
        Non_sign{tmp,:}=flag2;
        if number_dataset==1
            rownam{tmp}=['S' num2str(subject)];
        elseif number_dataset==2
            rownam{tmp}=['S' num2str(subject+8)];
        elseif number_dataset==3
            rownam{tmp}=['S' num2str(subject+9)];
        elseif number_dataset==4
            rownam{tmp}=['S' num2str(subject+10)];
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

%%%%%%%%%%%%%%%%%%%%%%%%
%average same asymmetry
%%%%%%%%%%%%%%%%%%%%%%%%

N_left_sessions(No_asymmetry)=[];
N_used_ses(No_asymmetry)=[];
N_tot_ses(No_asymmetry)=[];
rownam(No_asymmetry)=[];

%stab4(No_asymmetry)=[];
%stab5(No_asymmetry)=[];
stab6(No_asymmetry)=[];

%on average, 71% of sessions...
mean(stab6);

%In total 94% if sessions showed significant lateralization
Prop_lateral=sum(N_used_ses)/sum(N_tot_ses);

N_p_values=myBinomTest(N_left_sessions,N_used_ses,0.5,'two');
disp(round(N_p_values,4));

[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(N_p_values,0.05,'pdep','yes');
disp(round(adj_p,3));


tmp_mean_diff_group_subject_level=mean_diff_group_subject_level;
tmp_mean_diff_group_subject_level(No_asymmetry)=[];

for subj=1:length(N_left_sessions)
    if tmp_mean_diff_group_subject_level(subj)<0
        N_ses_equal_to_subj(subj)=N_used_ses(subj)-N_left_sessions(subj);
    else
        N_ses_equal_to_subj(subj)=N_left_sessions(subj);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_tot_ses=N_tot_ses([9:10 1:8 11:end]);
N_used_ses=N_used_ses([9:10 1:8 11:end]);
Prop_used_ses=round(N_used_ses./N_tot_ses,2)*100;
N_ses_equal_to_subj=N_ses_equal_to_subj([9:10 1:8 11:end]);
Effect_size=round(N_ses_equal_to_subj./N_used_ses,3)*100;
N_p_values=round(N_p_values([9:10 1:8 11:end]),3);
adj_p=round(adj_p([9:10 1:8 11:end]),3);


varnam={'Number_Asym_Sessions','prop_used_sessions','Effect_size', 'Unc_p_value','FDR_cor_p_value'}';
T=table(N_used_ses',Prop_used_ses',Effect_size',N_p_values',adj_p','RowNames',rownam','VariableNames',varnam);
        N_ses_equal_to_subj(subj)=N_used_ses(subj)-N_left_sessions(subj);
    else
        N_ses_equal_to_subj(subj)=N_left_sessions(subj);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_tot_ses=N_tot_ses([9:10 1:8 11:end]);
N_used_ses=N_used_ses([9:10 1:8 11:end]);
Prop_used_ses=round(N_used_ses./N_tot_ses,2)*100;
N_ses_equal_to_subj=N_ses_equal_to_subj([9:10 1:8 11:end]);
Effect_size=round(N_ses_equal_to_subj./N_used_ses,3)*100;
N_p_values=round(N_p_values([9:10 1:8 11:end]),3);
adj_p=round(adj_p([9:10 1:8 11:end]),3);


varnam={'Number_Asym_Sessions','prop_used_sessions','Effect_size', 'Unc_p_value','FDR_cor_p_value'}';
T=table(N_used_ses',Prop_used_ses',Effect_size',N_p_values',adj_p','RowNames',rownam','VariableNames',varnam);

cd([Work_dir '/Figures_paper_variability/Supplementary_Figures/']);

writetable(T,'Table_2','WriteRowNames',true);
disp(T);

end