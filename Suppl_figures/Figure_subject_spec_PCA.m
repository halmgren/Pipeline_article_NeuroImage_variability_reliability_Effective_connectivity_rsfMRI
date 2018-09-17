function Figure_subject_spec_PCA(SPM_dir,Work_dir)

    load([Work_dir '/Results_paper_variability/DCM/Basic/Smith/Full_model/PEB_group/PCA_subject_spec.mat']);
    PCA_coef=[coeff{1}(:,1), coeff{2}(:,1),coeff{3}(:,1),  coeff{4}(:,1), coeff{5}(:,1), coeff{6}(:,1), coeff{7}(:,1), coeff{8}(:,1), coeff{9}(:,1), coeff{10}(:,1), coeff{11}(:,1), coeff{12}(:,1), coeff{13}(:,1), coeff{14}(:,1), coeff{15}(:,1), coeff{16}(:,1), coeff{17}(:,1)];
    PCA_explained=[explained{1}(1) explained{2}(1) explained{3}(1) explained{4}(1) explained{5}(1) explained{6}(1) explained{7}(1) explained{8}(1) explained{9}(1) explained{10}(1) explained{11}(1) explained{12}(1) explained{13}(1) explained{14}(1) explained{15}(1) explained{16}(1) explained{17}(1)];
    cd([Work_dir '/Figures_paper_variability/']);
    mkdir('Supplementary_Figures');
    
    PCA_coef=round(PCA_coef,2);
    PCA_explained=round(PCA_explained,1);
    PCA_info=[PCA_coef;PCA_explained];
    dlmwrite([Work_dir '/Figures_paper_variability/Supplementary_Figures/Table_PCA_Subject_spec.txt'],PCA_info,'delimiter',' ');

end