function IL_plot_elecs(subject_block,comparison)
ERP_HFA_TF = 'elec_TF';
to_comp = find(cellfun(@iscell,comparison));

% preallocation
sbjs_elecs = cell(size(subject_block,1),2);all_classes=[];
for sb = 1:size(subject_block,1)
    sbj_ID = subject_block{sb,1};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    save_folder = fullfile(Sbj_Metadata.results,'power_peaks');
    curr_block= subject_block{sb,2};
    load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']))
    
    % preallocation
    sbjs_elecs{sb,1} = Sbj_Metadata.fsname;
    sbjs_elecs{sb,2} = get_info_goodchans(info);
    classes=cell(length(sbjs_elecs{sb,2}),1);
    
    % loop the electrodes
    for el = 1:length(sbjs_elecs{sb,2})
        elec = info.channelinfo.Label{el};
        load(fullfile(save_folder,[elec , '_' curr_block '_' ERP_HFA_TF '_powerstats.mat']),'ranksumresults')
        
        classes(el,:) = IL_ranksum_class(ranksumresults,comparison);
    end
    
    all_classes = [all_classes;classes];
    
end

elecClassColors = IL_classColors(all_classes,comparison);






[AllSubElecNames,AllSubElecCoords] = gather_allelecsinfo(sbjs_elecs,'FSAVERAGE');

color_elecs_wData('fsaverage',all_classes,AllSubElecNames,AllSubElecCoords,'discrete',[3 18],'t-values','inferno',edgeColoredElecs)


end