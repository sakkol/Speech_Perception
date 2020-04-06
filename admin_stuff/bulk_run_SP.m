%% Select 'subjects'
data_root = '/media/sakkol/HDD1/HBML/';
project_name = 'Speech_Perception';

AllBlockInfo = readtable(fullfile(data_root,'PROJECTS_DATA',project_name,[project_name '_BlockInfo.xlsx']));
subjects = unique(AllBlockInfo.sbj_ID);

[indx,~] = listdlg('ListString',subjects);

%% Behavioral
% create subjects in the first section
for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    
    for c = 1:size(Sbj_Metadata.BlockLists,1)
        
        curr_block = Sbj_Metadata.BlockLists{c};
        fprintf('Running for Subject:%s, Block:%s of %d\n',sbj_ID,curr_block,size(Sbj_Metadata.BlockLists,1))
        
        SP_beh_analysis(Sbj_Metadata,curr_block)
        close all
    end
end

%% Bulk run word evoked ERP-HFA plots
for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    SP_word_evoked(Sbj_Metadata)
end

%% SP_Efields
SP_root = fullfile(data_root,'PROJECTS_DATA',project_name);
SP_res = fullfile(SP_root,'Collected_Results','Efields');
if ~exist(SP_res,'dir'),mkdir(SP_res),end

AllBlockInfo = readtable(fullfile(SP_root,[project_name '_BlockInfo.xlsx']));
aroundPeak = 1;

sbj_IDs = unique(AllBlockInfo.sbj_ID(~ismember(AllBlockInfo.sbj_ID,'NS144_2')));
for s = 1:length(sbj_IDs)
    sbj_ID = sbj_IDs{s};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    
    whichblocks = AllBlockInfo.BlockList(ismember(AllBlockInfo.sbj_ID,sbj_ID) & AllBlockInfo.preproc_FU==1);
    for b = 1:length(whichblocks)
        curr_block = whichblocks{b};
        EA_efields(Sbj_Metadata,curr_block,aroundPeak)
        to_print_folder = fullfile(Sbj_Metadata.results,'Efields',curr_block);
        if ~aroundPeak
            copyfile(fullfile(to_print_folder,[curr_block, '_efield_plot.jpg']),fullfile(SP_res,[sbj_ID '_' curr_block '_efield_plot.jpg']));
        else
            copyfile(fullfile(to_print_folder,[curr_block, '_efield_plot_aroundPeak.jpg']),fullfile(SP_res,[sbj_ID '_' curr_block '_efield_plot_aroundPeak.jpg']));
        end
        close all
    end
end

