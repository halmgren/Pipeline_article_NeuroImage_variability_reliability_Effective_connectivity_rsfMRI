# Pipeline_article_NeuroImage_variability_reliability_Effective_connectivity_rsfMRI

This repository contains the code that was used for the analyses of the article concerning variability and reliability of effective connectivity estimated with DCM for resting state fMRI (https://doi.org/10.1016/j.neuroimage.2018.08.053)

The main pipeline is located in the 'Main_code'-folder. The instructions for reproduction (and data availability) are included in the 'Pipeline_paper_variability.m'-file, and the pipeline can be applied by executing the same file. 

The input files are included in the four datasets ('Myconnectome', 'Midnight Scan Club', 'Day2day' and 'Kirby') and should be converted to BIDS format (if not downloaded in this format). The code does all analyses in the paper, and creates figures and the result section. For some figures, other toolboxes should be downloaded, which are described in the 'Pipeline_paper_variability.m'-file. Differences in software and hardware versions can have an effect on results and inference.

Special thanks to Soroor Shafeizadegan (Isfahan University of Technology) for updating the code to a newer revision of SPM12 (r7487). For use of the code with this SPM12 revision, uncomment lines 120 and 128, and comment lines 122 and 130 in 'Extract_regressors_paper_variability.m'. In addition, uncomment lines 236 and 363, and comment lines 238 and 366 in 'Extract_timeseries_paper_variability.m'

This code is provided without any warranty. 

For feedback and comments, please contact Hannes Almgren (Hannes.Almgren@ugent.be)
