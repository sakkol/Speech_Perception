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

%% Bulk run separate fourier v2
% freq_bands = {[1 4],[4 8],[8 12]};
freq_bands = {[1 12]};
for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    SP_separate_event_fourier_v2(Sbj_Metadata,freq_bands)
    SP_control_ITPC(Sbj_Metadata)
%     SP_separate_event_fourier_v2(Sbj_Metadata)
end

%% Correcting info.events based on new info (peakRate and peakEnv)

AllBlockInfo = readtable(fullfile(data_root,'PROJECTS_DATA',project_name,[project_name '_BlockInfo.xlsx']));

sbj_IDs = unique(AllBlockInfo.sbj_ID(~ismember(AllBlockInfo.sbj_ID,'NS144_2')));
for s = 1:length(sbj_IDs)
    sbj_ID = sbj_IDs{s};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    
    whichblocks = AllBlockInfo.BlockList(ismember(AllBlockInfo.sbj_ID,sbj_ID));
    for b = 1:length(whichblocks)
        curr_block = whichblocks{b};
        fprintf('Running for pt-%s and block-%s\n',sbj_ID,curr_block)
        load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_info.mat']))
        
        for ii = 1:size(info.events,1)
            info.events.peak_info{ii} = all_info_table.peak_info{strcmp(all_info_table.sentence,info.events.sentence{ii}) ...
                & strcmp(all_info_table.rate,info.events.rate{ii})};
        end
        save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_info.mat']),'info')
    end
end

%% Bulk run separate fourier v3
% freq_bands = {[1 4],[4 8],[8 12]};
freq_bands = {[1 12]};
for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    SP_separate_event_fourier_v3(Sbj_Metadata,freq_bands)
%     SP_control_ITPC_peaks(Sbj_Metadata)
end


%% Bulk run peak events evoked ERP-HFA plots

for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    SP_event_evoked(Sbj_Metadata,[],0)
end

%% Save channel of interests

for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    
    whichblocks = AllBlockInfo.BlockList(ismember(AllBlockInfo.sbj_ID,sbj_ID));
    curr_block = whichblocks{1};
    
    load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_info.mat']))
    
    channel_OI = info.channelinfo.Label(contains(info.channelinfo.Label,{'LTs','RTs'}));
    
    save(fullfile(Sbj_Metadata.sbjDir,[Sbj_Metadata.sbj_ID, '_channel_OI.mat']),'channel_OI')
end

%% Bulk run peak events evoked ERP-HFA plots

for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    SP_control_spectrogram(Sbj_Metadata,[],0)
end
