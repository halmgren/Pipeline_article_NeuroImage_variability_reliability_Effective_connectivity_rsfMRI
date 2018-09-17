function table_PCA_paper_variability(SPM_dir,Work_dir)

%out of convenience, subjects are numbered from 1 to 17 (iso 20) the same order was used
load([Work_dir '/Results_paper_variability/DCM/Basic/Smith/Full_model/PEB_group/PCA_subject_spec.mat'],'coeff','score','latent','tsquared','explained','mu');
for subject=1:length(coeff)
    first_coeff_group(:,subject)=coeff{subject}(:,1);
    eval(sprintf('S%d = [coeff{subject}(:,1); explained{subject}(1)]', subject));
    %A(:,subject)=[coeff{subject}(:,1); explained{subject}(1)]';
    if subject<9
        A(:,subject+2)=[coeff{subject}(:,1); explained{subject}(1)]';
    elseif subject==9||subject==10
        A(:,subject-8)=[coeff{subject}(:,1); explained{subject}(1)]';
    else
        A(:,subject)=[coeff{subject}(:,1); explained{subject}(1)]';
    end
        
    varname{subject}=['S' num2str(subject)];
end

rownam={'con1','con2','con3','con4','con5','con6','con7','con8','con9','con10','con11','con12','con13','con14','con15','con16','expl_var'}';

T=array2table(round(A,2),'RowNames',rownam,'VariableNames',varname);

disp(T);

cd([Work_dir '/Figures_paper_variability/Supplementary_Figures/']);

writetable(T,'PCA_subject_spec.txt','WriteRowNames',true);
end