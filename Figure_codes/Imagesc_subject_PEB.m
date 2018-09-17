function Imagesc_subject_PEB(SPM_dir,Work_dir)
%%%%%%%%%%%%%%%%%
%GROUP_LEVEL PEB
%%%%%%%%%%%%%%%%%
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
    
    for number_procedure=1:length(all_procedure_names)
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
                        else
                            PEB_group{tmp}=PEB;
                        end
                        
                        clear PEB;
                    end
                end
            end
        end
    end
end

for subj=1:length(PEB_group)
    
    load([Work_dir '/DatasetKuehn/RestingStatefMRI/Subj_01_summary/DCM/Basic/Smith/Full_model/GCM_DMN_full_estim.mat']);
    DCM=GCM;
    
    PEB=PEB_group{subj};
    figure('units','normalized','outerposition',[0 0 1 1]);
    
    ci=spm_invNcdf(1-0.05);
    EP=full(vec2mat(PEB.Ep(1:16),4)');
    CP=diag(PEB.Cp);
    CP=full(vec2mat(CP(1:16),4)');
    sgn=sign(EP-ci*sqrt(CP)).*sign(EP+ci*sqrt(CP));
    
    A_matrix=full(vec2mat(PEB.Ep(1:16),4)');
    
    A_matrix(sgn==-1)=NaN;
    
    h=imagesc(A_matrix);
    set(h,'alphadata',~isnan(A_matrix))
    
    colorbar;
    
    axis square;
    xlabel('\textbf{\underline{From}}','FontSize',40,'Fontweight','bold','Interpreter','latex'); ylabel('\textbf{\underline{To}}','FontSize',40,'Fontweight','bold','Interpreter','latex');
    set(gca,'XAxisLocation','top','Tag','connectivity');
    
    title_str = ['Connectivity'];
    
    title(title_str,'FontSize',46);
    
    
    for region=1:length({DCM{1}.xY.name})
        regions(region)={DCM{1}.xY(region).name(5:end)};
    end
    
    if size(A_matrix,1) == size(A_matrix,2) && ...
            size(A_matrix,1) == length({DCM{1}.xY.name})
        set(gca,'YTickLabel',regions,...
            'YTick',1:length({DCM{1}.xY.name}),...
            'XTickLabel',regions,'fontweight','bold','fontsize',40,'XTick',...
            1:length({DCM{1}.xY.name}),'TickLabelInterpreter', 'none');
    end
    
    for side1=1:4
        for side2=1:4
            if ~isnan(A_matrix(side2,side1))
                text(side1,side2,num2str(round(A_matrix(side2,side1),2)),'HorizontalAlignment','center','Fontweight','bold','Fontsize',40);
                
            elseif isnan(A_matrix(side2,side1))
                text(side1,side2,'ns','HorizontalAlignment','center','Fontweight','bold','Fontsize',40);
            end
        end
    end
    
    %make sure that black values are seen on color background
    cmap_custom=parula;
    cmap_custom(1:4,:)=[];
    colormap(cmap_custom);
    
    set(gcf,'PaperPositionMode','auto');
    saveas(gcf,[Work_dir '/Figures_paper_variability/Supplementary_Figures/Subject_level_PEB_Ep_subject_' num2str(subj) '.bmp']);
    
    close;
    clear DCM PEB;
end
end