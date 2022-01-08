%% Select 'subjects'
data_root = '/media/sakkol/HDD1/HBML/';
project_name = 'IsochronousListening';

AllBlockInfo = readtable(fullfile(data_root,'PROJECTS_DATA',project_name,[project_name '_BlockInfo.xlsx']));
subjects = unique(AllBlockInfo.sbj_ID);

[indx,~] = listdlg('ListString',subjects);

%% Run the wavelet
for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    
    
    whichblocks = AllBlockInfo.BlockList(ismember(AllBlockInfo.sbj_ID,sbj_ID) & AllBlockInfo.preproc_FU==1);
    for b = 1:length(whichblocks)
        curr_block = whichblocks{b};
        %         if ~exist(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_wlt.mat']),'dir')
        fprintf([sbj_ID,'-',curr_block,'\n'])
        load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_avg.mat']))
        load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']))
        events = info.events; clear info
        
        pre  = 1.5; % seconds
        post = 14; % seconds
        [epoched_data, epoched_wlt] = master_waveletTF(ecog_avg.ftrip,events.speech_onsets,pre,post,[0.1 200]);
        
        % save wlt and data
        save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_wlt_pow.mat']),'epoched_wlt','epoched_data','-v7.3');
        clear epoched_data epoched_wlt
        %         end
        
    end
end

%% Quick Plot
% create subjects in the first section
for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    
    whichblocks = AllBlockInfo.BlockList(ismember(AllBlockInfo.sbj_ID,sbj_ID) & AllBlockInfo.preproc_FU==1);
    for b = 1:length(whichblocks)
        curr_block = whichblocks{b};
        fprintf('%s - %s\n',sbj_ID,curr_block)
        IL_quickPlot(Sbj_Metadata,curr_block)
    end
    
end

%% not wavelet, multitaper
for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    
    
    whichblocks = AllBlockInfo.BlockList(ismember(AllBlockInfo.sbj_ID,sbj_ID) & AllBlockInfo.preproc_FU==1);
    for b = 1:length(whichblocks)
        curr_block = whichblocks{b};
        %         if ~exist(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_wlt.mat']),'dir')
        fprintf([sbj_ID,'-',curr_block,'\n'])
        load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_avg.mat']))
        load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']))
        events = info.events; clear info
        
        pre  = 1.5; % seconds
        post = 14; % seconds
        foi = [0.2:0.1:6,50:7:200];
        [epoched_data, epoched_wlt] = master_mtmconvolTF(ecog_avg.ftrip,events.speech_onsets,pre,post,foi,'pow');
        
        % save wlt and data
        save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_mtmconvol_pow.mat']),'epoched_wlt','epoched_data','-v7.3');
        clear epoched_data epoched_wlt
        %         end
        
    end
end

%% run analyses on individual electrodes: MTMCONVOL
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
        
        pre  = 1.5; % seconds
        post = 14; % seconds
        foi = [0.2:0.05:5,50:5:200];
        parfor el = 1:length(ecog_avg.ftrip.label)
            elec = ecog_avg.ftrip.label{el};
            IL_mtmconvolTF(Sbj_Metadata,curr_block,ecog_avg.ftrip,speech_onsets,pre,post,elec,foi,'fourier');
            
            % save wlt and data
            IL_quickPlot_elec(Sbj_Metadata,curr_block,elec)
        end
    end
end

%% run analyses on individual electrodes: MTMFFT
loop1={'iso','a'};
loop2={4,3,5};
loop3={'sentence','scrambled'};

for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    
    whichblocks = AllBlockInfo.BlockList(ismember(AllBlockInfo.sbj_ID,sbj_ID) & AllBlockInfo.preproc_FU==1);
    for b = 1:length(whichblocks)
        curr_block = whichblocks{b};
        fprintf([sbj_ID,'-',curr_block,'\n'])
        load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_avg.mat']))
        load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']))
        load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_trial_nos.mat']))
        events = info.events; clear info
        isofrequency = events.cond_info{1}.Frequency;
        speech_onsets = events.speech_onsets;
        
        for l1=1:2
            for l2=1:3
                for l3=1:2
                    
                    pre  = -2+(loop2{l2}/isofrequency); % seconds
                    post = (loop2{l2}*5/isofrequency)-2; % seconds
                    foi = [0.2:0.1:5,50:5:200];
                    to_save_name = [loop1{l1} '_' num2str(loop2{l2}) '_' loop3{l3}];
                    
                    for el = 1:length(ecog_avg.ftrip.label)
                        elec = ecog_avg.ftrip.label{el};
                        IL_mtmfftTF(Sbj_Metadata,curr_block,ecog_avg.ftrip,speech_onsets(trial_nos.([loop1{l1} '_' num2str(loop2{l2}) '_' loop3{l3}]))...
                            ,pre,post,elec,foi,'fourier',to_save_name);
                        
                        % save wlt and data
                        IL_quickPlot_elec_fft(Sbj_Metadata,curr_block,elec,to_save_name)
                    end
                    
                    
                end
            end
        end
        
        
        
    end
end



%% run peak signifcance analyses on individual electrodes
for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    
    whichblocks = AllBlockInfo.BlockList(ismember(AllBlockInfo.sbj_ID,sbj_ID) & AllBlockInfo.preproc_FU==1);
    for b = 1:length(whichblocks)
        curr_block = whichblocks{b};
        fprintf([sbj_ID,'-',curr_block,'\n'])
        load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']))
        save(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_trial_nos.mat']),'trial_nos')
        
        parfor el = 1:length(info.channelinfo.Label)
            elec = info.channelinfo.Label{el};            
            % save wlt and data
            IL_power_peaks(Sbj_Metadata,curr_block,elec)
        end
    end
end

%% save trial numbers in a non-table format to be read by Python FOOOF script

for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    
    whichblocks = AllBlockInfo.BlockList(ismember(AllBlockInfo.sbj_ID,sbj_ID) & AllBlockInfo.preproc_FU==1);
    for b = 1:length(whichblocks)
        curr_block = whichblocks{b};
        fprintf([sbj_ID,'-',curr_block,'\n'])
        load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']))

% 3-4-5 -- iso-a -- sentence-scrambled -- wlt-hfa-erp
loop1={'iso','a'};
loop2={4,3,5};
loop3={'sentence','scrambled'};
loop4={'wlt','hfa','data'};
trial_nos=[];
for l1=1:2
    for l2=1:3
        for l3=1:2
            trial_nos.([loop1{l1} '_' num2str(loop2{l2}) '_' loop3{l3}]) = IL_get_eventsOI(info.events, loop3{l3}, loop1{l1}, loop2{l2});
        end
    end
end
save(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_trial_nos.mat']),'trial_nos')

    end
end