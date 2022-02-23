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

%% run the EFields
EA_root = fullfile(data_root,'PROJECTS_DATA',project_name);
EA_res = fullfile(EA_root,'Collected_Results','Efields');
if ~exist(EA_res,'dir'),mkdir(EA_res),end

aroundPeak = 1;
elecsOI = {'allchans','goodchans','selectchans','selectgoodchans'};

AllBlockInfo = readtable(fullfile(EA_root,[project_name '_BlockInfo.xlsx']));

% subjects = {'NS148','NS148_2','NS144_2','NS150','NS170'};
subjects = {'NS174_2'};

for e = 1:length(elecsOI)
    for s = 1:length(subjects)
        sbj_ID = subjects{s};
        Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
        
        whichblocks = AllBlockInfo.BlockList(ismember(AllBlockInfo.sbj_ID,sbj_ID) & AllBlockInfo.preproc_FU==1);
%         whichblocks = AllBlockInfo.BlockList(ismember(AllBlockInfo.sbj_ID,sbj_ID) & strcmp(AllBlockInfo.conditions_code,'1,8,9') & AllBlockInfo.preproc_FU==1);
        for b = 1:length(whichblocks)
            curr_block = whichblocks{b};
            BlockInfo = makeBlockInfo(Sbj_Metadata,curr_block);
            
            EA_efields(Sbj_Metadata,curr_block,aroundPeak,{[BlockInfo.StimSide{1} 'omni'],BlockInfo.StimSide{1}},elecsOI{e})
            close all
        end
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
    SP_control_spectrogram(Sbj_Metadata,[],0,0)
    SP_control_spectrogram(Sbj_Metadata,[],0,1)
end

%% Bulk run peak events evoked ITPC plots and stats

for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    SP_POS_peakevents(Sbj_Metadata)
end

for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    SP_POS_equalpeaks(Sbj_Metadata)
end

%% Comparing HFA results from wavelet and BLP

parfor s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    SP_BLP_comp(Sbj_Metadata,[],0,0)
end

%% New sep event fourier v4

freq_bands = {[1 200]};
for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    SP_separate_event_fourier_v4(Sbj_Metadata,freq_bands,'_all')
end

%% Run the greatest plots ever

for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    SP_wordvspeak_allplot(Sbj_Metadata)
end

for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    SP_wordvspeak_allplot_withwrong(Sbj_Metadata)
end

for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    SP_wordvssyll_allplot(Sbj_Metadata)
end

for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    SP_wordvssyll_CorrVsWrng(Sbj_Metadata)
end

%% Check the acoustic differences between corr no vs wrong

for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    SP_comp_resp_acoustics(Sbj_Metadata)
end

%% Correlating the SNR and Accuracy of control condition across patients
AllBlockInfo = readtable(fullfile(data_root,'PROJECTS_DATA',project_name,[project_name '_BlockInfo.xlsx']));
sbj_block = AllBlockInfo(contains(AllBlockInfo.conditions_code,'1')&AllBlockInfo.preproc_FU==1,1:2);
sbj_block = table2cell(sbj_block);

% filter for specific patients
sbj_block = sbj_block(ismember(sbj_block(:,1),subjects),1:2);
% filter for onlt the first block
sbj_tmp = unique(sbj_block(:,1));sbj_block_tmp = cell(length(sbj_tmp),2);
for s = 1:length(sbj_tmp)
    sbj_block_tmp{s,1} = sbj_tmp{s};
    x = sbj_block(strcmp(sbj_block(:,1),sbj_tmp{s}),2);
    sbj_block_tmp{s,2} = x{1};
end
sbj_block = sbj_block_tmp;


SP_AccrSNR_corr(sbj_block);

print(fullfile(data_root,'PROJECTS_DATA',project_name,'Collected_Results','AccrSNR_corr',[strjoin(sbj_block(:,1),'_') '_allblocks.png']),'-dpng','-r300')
print(fullfile(data_root,'PROJECTS_DATA',project_name,'Collected_Results','AccrSNR_corr',[strjoin(sbj_block(:,1),'_') '_firstblocks.png']),'-dpng','-r300')

%% Noise Invariance
% Idea: looking at the electrodes that respond solely to speech vs sounds
AllBlockInfo = readtable(fullfile(data_root,'PROJECTS_DATA',project_name,[project_name '_BlockInfo.xlsx']));
sbj_block = AllBlockInfo(AllBlockInfo.manuscript==1,1:2);
sbj_block = table2cell(sbj_block);

elecsOI = {'allElecs','HG','STG',{'HG','STG'}};
for elec = 1:length(elecsOI)
    for s=1:size(sbj_block,1)
        SP_noiseInvariance(sbj_block(s,:),elecsOI{elec})
    end
    SP_noiseInvariance(sbj_block,elecsOI{elec})
end

%% entrainment - accuracy correlation for attention sentence
% Idea: looking at the electrodes that respond solely to speech vs sounds
AllBlockInfo = readtable(fullfile(data_root,'PROJECTS_DATA',project_name,[project_name '_BlockInfo.xlsx']));
sbj_block = AllBlockInfo(AllBlockInfo.manuscript==1,1:2);
sbj_block = table2cell(sbj_block);

elecsOI = {'allElecs','HG','STG',{'HG','STG'}};
for elec = 1:length(elecsOI)
    for s=1:size(sbj_block,1)
        SP_TrackingAttention_Accuracy(sbj_block(s,:),elecsOI{elec})
    end
    SP_TrackingAttention_Accuracy(sbj_block,elecsOI{elec})
end


%% entrainment - accuracy correlation for target sentence
% Idea: looking at the electrodes that respond solely to speech vs sounds
AllBlockInfo = readtable(fullfile(data_root,'PROJECTS_DATA',project_name,[project_name '_BlockInfo.xlsx']));
sbj_block = AllBlockInfo(AllBlockInfo.manuscript==1,1:2);
sbj_block = table2cell(sbj_block);

elecsOI = {'allElecs','HG','STG',{'HG','STG'}};
for elec = 1:length(elecsOI)
    for s=1:size(sbj_block,1)
        SP_TrackingTarget_Accuracy(sbj_block(s,:),elecsOI{elec})
    end
    SP_TrackingTarget_Accuracy(sbj_block,elecsOI{elec})
end

