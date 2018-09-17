function Sign_WS_reliability_network_all(SPM_dir,Work_dir)

procedure='Basic';
name_ROI_def='Smith';


[ROI_list]=Define_ROIs_paper_variability(name_ROI_def);

tmp=0;
n=0;

for VOI_number=1:size(ROI_list,1)
    ntwrk=ROI_list{VOI_number,1}(1:3);
    
    if VOI_number>1 && strcmp(ROI_list{VOI_number,1}(1:3),ROI_list{VOI_number-1,1}(1:3))
        n=n+1;
        ntwrk_size(tmp)=ntwrk_size(tmp)+1;
        ntwrk_VOI_names{n,tmp}=ROI_list{VOI_number,1}(5:end);
        ntwrk_VOI_coord{n,tmp}=ROI_list{VOI_number,2};
        continue
        
    else
        n=1;
        tmp=tmp+1;
        ntwrk_size(tmp)=1;
        ntwrk_name{tmp}=ROI_list{VOI_number,1}(1:3);
        ntwrk_VOI_names{n,tmp}=ROI_list{VOI_number,1}(5:end);
        ntwrk_VOI_coord{n,tmp}=ROI_list{VOI_number,2};
    end
end

stability_threshold=0.75;

%%%%%%%%%%%%%%%%%%%%
%Figure of networks
%%%%%%%%%%%%%%%%%%%%
tmp=0;
for number_dataset=1:4
    
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_variability(number_dataset);
    for subject=1:number_subject
        clear PEB;
        tmp=tmp+1;
        F1=figure;
        if strcmp(dataset,'DatasetGordon')&&(subject==3||subject==8||subject==9)
            continue;
        end
        
        Plot_sign_WS_reliability_network_all(dataset,subject,ntwrk_size,ntwrk_VOI_names,ntwrk_VOI_coord,ntwrk_name{1},name_ROI_def,stability_threshold,procedure,SPM_dir,Work_dir);
    end
    close all;
end


end

