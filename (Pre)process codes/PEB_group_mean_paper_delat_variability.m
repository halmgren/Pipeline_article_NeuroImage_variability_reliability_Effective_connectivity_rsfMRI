function PEB_group_mean_paper_delat_variability(SPM_dir,Work_dir)

    
all_ROI_defs={'Smith'};
all_procedure_names={'Basic'};

for number_ROI_def=1:length(all_ROI_defs)
    
    name_ROI_def=all_ROI_defs{number_ROI_def};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Give regions name and coördinates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [ROI_list]=Define_ROIs_paper_variability(name_ROI_def);
    
    tmp1=0;
    
    for VOI_number=1:size(ROI_list,1)
        ntwrk=ROI_list{VOI_number,1}(1:3);
        
        if VOI_number>1 && strcmp(ROI_list{VOI_number,1}(1:3),ROI_list{VOI_number-1,1}(1:3))
            ntwrk_size(tmp1)=ntwrk_size(tmp1)+1;
            continue
            
        else
            tmp1=tmp1+1;
            ntwrk_size(tmp1)=1;
            ntwrk_name{tmp1}=ROI_list{VOI_number,1}(1:3);
        end
    end
    
    for number_procedure=1 %Basic analyses were used here
        procedure=all_procedure_names{number_procedure};
        
        if strcmp(procedure,'Basic')
            mkdir([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/']);
            
            for network_number=1:length(ntwrk_name)
                
                tmp=0;
                
                disp(network_number);
                for number_dataset=1:4
                    
                    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
                    for subject=1:number_subject
                        clear PEB;
                        tmp=tmp+1;
                        
                        try
                            load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_A_mean_' ntwrk_name{network_number} '.mat']);
                            load([Work_dir '/' dataset '/sub-' sprintf('%02d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/Lateralization_index_' ntwrk_name{network_number} '.mat']);
                        catch ME
                            tmp=tmp-1;
                            disp('PEB group')
                            disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
%                             pause;
                            if strcmp(ME.identifier,'MATLAB:load:couldNotReadFile')
                                continue;
                            else
                                save([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/ERROR_PEB_A_mean_' dataset '_subject_' num2str(subject) '_' ntwrk_name{network_number} '.mat'],'ME');
                                continue;
                            end
                        end
                        %We don't include participants with less then 8
                        %useful sessions.
                        if length(PEB.Snames)<8
                            tmp=tmp-1;
                            continue
                        elseif posterior_probability<0.95&&posterior_probability>0.05
                            tmp=tmp-1;
                            continue
                        else
                            if mean_diff>0
                                PEB_group{tmp}=PEB;
                            elseif mean_diff<0
                                EP=full(vec2mat(PEB.Ep(1:16),4)');
                                EP([2 4],:)=EP([4 2],:);
                                EP(:,[2 4])=EP(:,[4 2]);
                                
                                PEB.Ep=sparse(spm_vec(EP));
                                
                                EH=full(vec2mat(PEB.Eh(1:16),4)');
                                EH([2 4],:)=EH([4 2],:);
                                EH(:,[2 4])=EH(:,[4 2]);
                                
                                PEB.Cp(:,:)=PEB.Cp(:,[1 4 3 2 13 16 15 14 9 12 11 10 5 8 7 6]);
                                PEB.Cp(:,:)=PEB.Cp([1 4 3 2 13 16 15 14 9 12 11 10 5 8 7 6],:);
                                 
                                PEB.Ch(:,:)=PEB.Ch(:,[1 4 3 2 13 16 15 14 9 12 11 10 5 8 7 6]);
                                PEB.Ch(:,:)=PEB.Ch([1 4 3 2 13 16 15 14 9 12 11 10 5 8 7 6],:);
                                
                                PEB.Ce(:,:)=PEB.Ce(:,[1 4 3 2 13 16 15 14 9 12 11 10 5 8 7 6]);
                                PEB.Ce(:,:)=PEB.Ce([1 4 3 2 13 16 15 14 9 12 11 10 5 8 7 6],:);
                                
                                PEB_group{tmp}=PEB;
                            end
                        end
                        
                        clear PEB mean_diff posterior_probability;
                    end
                end
                
                
                
                M = struct();
                M.alpha = 1;
                M.beta  = 16;
                M.hE    = 0;
                M.hC    = 1/16;
                M.Q     = 'all';
                M.X=ones(length(PEB_group),1);
                
                [PEB DCM]=spm_dcm_peb_of_peb(PEB_group',M,'A');
                save([Work_dir '/Results_paper_variability/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_A_mean_group_delat_' ntwrk_name{network_number} '.mat'],'PEB','DCM','PEB_group');
                
                clear PEB_group;
            end
        end
    end
end

end