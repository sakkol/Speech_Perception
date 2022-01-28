function IL_plot_elecs(subject_block,comparison)
ERP_HFA_TF = 'elec_TF';
% to_comp = find(cellfun(@iscell,comparison));
data_root = '/media/sakkol/HDD1/HBML/';
project_name = 'IsochronousListening';

% preallocation
sbjs_elecs = cell(size(subject_block,1),2);
all_classes=[];
for sb = 1:size(subject_block,1)
    sbj_ID = subject_block{sb,1};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    save_folder = fullfile(Sbj_Metadata.results,'power_peaks');
    curr_block= subject_block{sb,2};
    load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']),'info')
    
    % preallocation
    sbjs_elecs{sb,1} = Sbj_Metadata.fsname;
    sbjs_elecs{sb,2} = get_info_goodchans(info);
    classes=repmat({{'Non-resp/Other'}},length(sbjs_elecs{sb,2}),1);
    
    % loop the electrodes to get significant classifications
    for el = 1:length(sbjs_elecs{sb,2})
        elec = info.channelinfo.Label{el};
        load(fullfile(save_folder,[elec , '_' curr_block '_' ERP_HFA_TF '_powerstats.mat']),'ranksumresults')
        t = IL_ranksum_SignClass(ranksumresults);
        if ~isempty(t),classes{el} = t;end
    end
    
    % remove unresp/other electrodes
    nons = ~cellfun(@(x)isequal(x,1),cellfun(@(x)strcmp(x,'Non-resp/Other'),classes,'UniformOutput',0));
    classes = classes(nons);
    sbjs_elecs{sb,2} = sbjs_elecs{sb,2}(nons);
    
    all_classes = [all_classes;classes];
end

% assign color to each electrode
[elec_colors,elec_classes] = IL_classColors(all_classes,comparison);
[AllSubElecNames,AllSubElecCoords] = gather_allelecsinfo(sbjs_elecs,'FSAVERAGE');
nons = ~cellfun(@(x)isequal(x,1),cellfun(@(x)strcmp(x,'Non-resp/Other'),elec_classes,'UniformOutput',0));
elec_colors = elec_colors(nons,:);
elec_classes = elec_classes(nons);
AllSubElecNames = AllSubElecNames(nons);
AllSubElecCoords = AllSubElecCoords(nons,:);

IL_color_elecs_wData('fsaverage',elec_classes,AllSubElecNames,AllSubElecCoords,'discrete',[],'Classes',elec_colors,[])
totitle = char(join(string(cellfun(@(x)join(x,'&'),comparison,'UniformOutput',0)),' : '));
text(gca,.5,1.07,['Any combiation of ' totitle],'Units','normalized','FontSize',18,'FontWeight','bold','HorizontalAlignment','center')

% save
toname = char(join(string(cellfun(@(x)join(x,'&'),comparison,'UniformOutput',0)),'_'));
save_folder = '/media/sakkol/HDD1/HBML/PROJECTS_DATA/IsochronousListening/Collected_Results';
print(fullfile(save_folder,[toname '_signElecs.jpg']),'-djpeg','-r300')

end