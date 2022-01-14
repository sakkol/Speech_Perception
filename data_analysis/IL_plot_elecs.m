function IL_plot_elecs(subject_block)
ERP_HFA_TF = 'elec_TF';

for sb = 1:height(subject_block)
    sbj_ID = subject_block{sb,1};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    save_folder = fullfile(Sbj_Metadata.results,'power_peaks');
    curr_block= subject_block{sb,2};
    load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']))
    for el = 1:length(info.channelinfo.Label)
            elec = info.channelinfo.Label{el};            
            % save wlt and data
            load(fullfile(save_folder,[elec , '_' curr_block '_' ERP_HFA_TF '_powerstats.mat']),'ranksumresults')
            
%             structfun
    end
end



end
