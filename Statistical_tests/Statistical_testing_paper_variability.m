function Statistical_testing_paper_variability(SPM_dir,Work_dir)
    
    %Lateralization permutation test
    Permutation_test_difference_subjects_paired_paper_variability(SPM_dir,Work_dir)
    permutation_test_difference_sessions_paired_paper_variability(SPM_dir,Work_dir)
    
    %BMR
    BMR_paired_permutation_subject_paper_variability(SPM_dir,Work_dir)
    BMR_paired_permutation_session_paper_variability(SPM_dir,Work_dir)
    
    %ROI size
    ROI_paired_permutation_subject_paper_variability(SPM_dir,Work_dir)
    ROI_paired_permutation_session_paper_variability(SPM_dir,Work_dir)
    
    %GSR
    GSR_paired_permutation_subject_paper_variability(SPM_dir,Work_dir)
    GSR_paired_permutation_session_paper_variability(SPM_dir,Work_dir)
end