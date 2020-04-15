%% Only copies plots
data_root = '/media/sakkol/HDD1/HBML/';
project_name = 'Speech_Perception';
SP_root = fullfile(data_root,'PROJECTS_DATA',project_name);
SP_res = fullfile(SP_root,'Collected_Results','Efields');
if ~exist(SP_res,'dir'),mkdir(SP_res),end

AllBlockInfo = readtable(fullfile(SP_root,[project_name '_BlockInfo.xlsx']));


sbj_IDs = unique(AllBlockInfo.sbj_ID(~ismember(AllBlockInfo.sbj_ID,'NS144_2')));
for s = 1:length(sbj_IDs)
    sbj_ID = sbj_IDs{s};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    
    for b = 1:length(Sbj_Metadata.BlockLists)
        curr_block = Sbj_Metadata.BlockLists{b};
        to_print_folder = fullfile(Sbj_Metadata.results,'Efields',curr_block);
        copyfile(fullfile(to_print_folder,[curr_block, '_efield_plot.jpg']),fullfile(SP_res,[sbj_ID '_' curr_block '_efield_plot.jpg']));
    end
    
end

%% Runs analysis and copies plots
data_root = '/media/sakkol/HDD1/HBML/';
project_name = 'Speech_Perception';
SP_root = fullfile(data_root,'PROJECTS_DATA',project_name);
SP_res = fullfile(SP_root,'Collected_Results','Efields');
if ~exist(SP_res,'dir'),mkdir(SP_res),end

aroundPeak = 1;

AllBlockInfo = readtable(fullfile(SP_root,[project_name '_BlockInfo.xlsx']));


sbj_IDs = unique(AllBlockInfo.sbj_ID(~ismember(AllBlockInfo.sbj_ID,'NS144_2'))); % Spanish version
for s = 1:length(sbj_IDs)
    sbj_ID = sbj_IDs{s};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    
    for b = 1:length(Sbj_Metadata.BlockLists)
        curr_block = Sbj_Metadata.BlockLists{b};
        SP_efields(Sbj_Metadata,curr_block,aroundPeak)
        to_print_folder = fullfile(Sbj_Metadata.results,'Efields',curr_block);
        if ~aroundPeak
            copyfile(fullfile(to_print_folder,[curr_block, '_efield_plot.jpg']),fullfile(SP_res,[sbj_ID '_' curr_block '_efield_plot.jpg']));
        else
            copyfile(fullfile(to_print_folder,[curr_block, '_efield_plot_aroundPeak.jpg']),fullfile(SP_res,[sbj_ID '_' curr_block '_efield_plot_aroundPeak.jpg']));
        end
        close all
    end
    
end