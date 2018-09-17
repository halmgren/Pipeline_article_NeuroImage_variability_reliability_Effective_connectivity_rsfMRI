# Pipeline_article_NeuroImage_variability_reliability_Effective_connectivity_rsfMRI

This repository contains the code that was used for the analyses of the article concerning variability and reliability of effective connectivity estimated with DCM for resting state fMRI (https://doi.org/10.1016/j.neuroimage.2018.08.053)

The main pipeline is located in the 'Main_code'-folder. The instructions for reproduction (and data availability) are included in the 'Pipeline_paper_variability.m'-file, and the pipeline can be applied by executing the same file. 

The input files are included in the four datasets ('Myconnectome', 'Midnight Scan Club', 'Day2day' and 'Kirby') and should be converted to BIDS format (if not downloaded in this format). The code does all analyses in the paper, and creates figures and the result section. For some figures, other toolboxes should be downloaded, which are described in the 'Pipeline_paper_variability.m'-file. Differences in software and hardware versions can have an effect on results and inference.

This code is provided without any warranty. 

For feedback and comments, please contact Hannes Almgren (Hannes.Almgren@ugent.be)
