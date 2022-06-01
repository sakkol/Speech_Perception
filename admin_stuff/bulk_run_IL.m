%% Select 'subjects'
data_root = '/media/sakkol/HDD1/HBML/';
project_name = 'IsochronousListening';

AllBlockInfo = readtable(fullfile(data_root,'PROJECTS_DATA',project_name,[project_name '_BlockInfo.xlsx']));
subjects = unique(AllBlockInfo.sbj_ID);

[indx,~] = listdlg('ListString',subjects);

% % % % % %% Run the wavelet
% % % % % for s = 1:length(indx)
% % % % %     sbj_ID = subjects{indx(s)};
% % % % %     Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
% % % % %     
% % % % %     whichblocks = AllBlockInfo.BlockList(ismember(AllBlockInfo.sbj_ID,sbj_ID) & AllBlockInfo.preproc_FU==1);
% % % % %     for b = 1:length(whichblocks)
% % % % %         curr_block = whichblocks{b};
% % % % %         %         if ~exist(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_wlt.mat']),'dir')
% % % % %         fprintf([sbj_ID,'-',curr_block,'\n'])
% % % % %         load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_avg.mat']))
% % % % %         load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']))
% % % % %         events = info.events; clear info
% % % % %         
% % % % %         pre  = 1.5; % seconds
% % % % %         post = 14; % seconds
% % % % %         [epoched_data, epoched_wlt] = master_waveletTF(ecog_avg.ftrip,events.speech_onsets,pre,post,[0.1 200]);
% % % % %         
% % % % %         % save wlt and data
% % % % %         save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_wlt_pow.mat']),'epoched_wlt','epoched_data','-v7.3');
% % % % %         clear epoched_data epoched_wlt
% % % % %         %         end
% % % % %         
% % % % %     end
% % % % % end
% % % % 
% % % % % %% Quick Plot
% % % % % % create subjects in the first section
% % % % % for s = 1:length(indx)
% % % % %     sbj_ID = subjects{indx(s)};
% % % % %     Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
% % % % %     
% % % % %     whichblocks = AllBlockInfo.BlockList(ismember(AllBlockInfo.sbj_ID,sbj_ID) & AllBlockInfo.preproc_FU==1);
% % % % %     for b = 1:length(whichblocks)
% % % % %         curr_block = whichblocks{b};
% % % % %         fprintf('%s - %s\n',sbj_ID,curr_block)
% % % % %         IL_quickPlot(Sbj_Metadata,curr_block)
% % % % %     end
% % % % %     
% % % % % end
% % % % 
% % % % % %% not wavelet, multitaper
% % % % % for s = 1:length(indx)
% % % % %     sbj_ID = subjects{indx(s)};
% % % % %     Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
% % % % %     
% % % % %     
% % % % %     whichblocks = AllBlockInfo.BlockList(ismember(AllBlockInfo.sbj_ID,sbj_ID) & AllBlockInfo.preproc_FU==1);
% % % % %     for b = 1:length(whichblocks)
% % % % %         curr_block = whichblocks{b};
% % % % %         %         if ~exist(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_wlt.mat']),'dir')
% % % % %         fprintf([sbj_ID,'-',curr_block,'\n'])
% % % % %         load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_avg.mat']))
% % % % %         load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']))
% % % % %         events = info.events; clear info
% % % % %         
% % % % %         pre  = 1.5; % seconds
% % % % %         post = 14; % seconds
% % % % %         foi = [0.2:0.1:6,50:7:200];
% % % % %         [epoched_data, epoched_wlt] = master_mtmconvolTF(ecog_avg.ftrip,events.speech_onsets,pre,post,foi,'pow');
% % % % %         
% % % % %         % save wlt and data
% % % % %         save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_mtmconvol_pow.mat']),'epoched_wlt','epoched_data','-v7.3');
% % % % %         clear epoched_data epoched_wlt
% % % % %         %         end
% % % % %         
% % % % %     end
% % % % % end

%% run TF analyses on individual electrodes using MTMCONVOL and plot the results
for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    
    whichblocks = AllBlockInfo.BlockList(ismember(AllBlockInfo.sbj_ID,sbj_ID) & AllBlockInfo.preproc_FU==1);
    for b = 1:length(whichblocks)
        curr_block = whichblocks{b};
        fprintf([sbj_ID,'-',curr_block,'\n'])
        load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_avg.mat']))
        load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']))
        events = info.events; clear info
        speech_onsets = events.speech_onsets;
        
        pre  = 6.5; % seconds
        post = 15; % seconds
        foi = [-0.1+1/6:1/12:5,50:5:200];
        parfor el = 1:length(ecog_avg.ftrip.label)
            elec = ecog_avg.ftrip.label{el};
            IL_mtmconvolTF(Sbj_Metadata,curr_block,ecog_avg.ftrip,speech_onsets,pre,post,elec,foi,'fourier');
            
            % save wlt and data
            IL_quickPlot_elec(Sbj_Metadata,curr_block,elec)
        end
    end
end


%% run peak signifcance analyses on individual electrodes
for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    
    whichblocks = AllBlockInfo.BlockList(ismember(AllBlockInfo.sbj_ID,sbj_ID) & AllBlockInfo.preproc_FU==1);
    if isempty(whichblocks),continue,end
    for b = 1:length(whichblocks)
        curr_block = whichblocks{b};
        fprintf([sbj_ID,'-',curr_block,'\n'])
        load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']))
        
        parfor el = 1:length(info.channelinfo.Label)
            elec = info.channelinfo.Label{el};            
            % save wlt and data
            IL_power_peaks(Sbj_Metadata,curr_block,elec)
        end
    end
end


%% plot peak significance analyses on common brain
subject_block={};
for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    
    whichblocks = AllBlockInfo.BlockList(ismember(AllBlockInfo.sbj_ID,sbj_ID) & AllBlockInfo.preproc_FU==1);
    if isempty(whichblocks),continue,end
    for b = 1
        curr_block = whichblocks{b};
        
        subject_block{end+1,1} = sbj_ID;
        subject_block{end,2} = curr_block;
    end
end

response_types = {'W', 'P', 'S', 'WP', 'WS', 'PS' 'WPS'};
comparisons = {'iso','4','sentence',response_types;...                              % 1
               'iso','4','scrambled',response_types;...
               'iso','4',{'sentence','scrambled'},response_types;...
               {'iso','a'},'4','sentence',response_types;...                        % 4
               'iso',{'4','3','5'},'sentence',response_types;...
               'iso',{'4','3','5'},{'sentence','scrambled'},response_types;...
               'iso',{'4','3','5'},{'sentence','scrambled'},response_types(1:3);... % 7
               'iso','4',{'sentence','scrambled'},response_types(1:3);...
               'iso','4',{'sentence','scrambled'},response_types([2,3,6]);...
               'iso','4',{'sentence','scrambled'},response_types([1]);...           % 10
               'iso','4','sentence',response_types([1]);...
               'a','4',{'sentence','scrambled'},response_types;...
               'a','4',{'sentence','scrambled'},response_types(1:3)};
% plot comparisons
for c = 1:size(comparisons,1)
    comparison=comparisons(c,:);
    IL_plot_elecs(subject_block,comparison)
    close all
end