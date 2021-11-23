%% Prepare
data_root = '/media/sakkol/HDD1/HBML/';
project_name = 'IsochronousListening';
sbj_ID = 'NS174_2';
Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'

% Get params directly from BlockList excel sheet
curr_block = Sbj_Metadata.BlockLists{1}
params = AllBlockInfo2params(Sbj_Metadata,curr_block)

%% Import
if ~strcmpi(params.CurrBlockInfo.EEGDAT,'edf')
    [ecog] = TDT2ecog(params);
else
    % for edfs
    params.edf = fullfile(params.directory,[params.filename, '.edf']);
    params.analog.ttl            = 'DC1';
    params.analog.audio          = 'DC2';
    % params.analog.micro          = 'C210';
    ecog=edf2ecog(params);
    tmp = ecog.ftrip.time;ecog.ftrip.time={[]};ecog.ftrip.time{1}=tmp;clear tmp
    
    % not needed, just can stay here as a reminder: data = edf2fieldtrip(params.edf)
end

%% Find bad channels using PSD
% find_bad_chans(ecog);
ecog = bad_chan_GUI([],ecog);
if size(ecog.bad_chans,2)>1,ecog.bad_chans = ecog.bad_chans';end % if bad_chans are in row order.

% Look at EEG, write bad channels to Xls-sheet
% cfg=[];
% cfg.channel = ecog.bad_chans;
% ecog_databrowser(ecog)

%% Check the good channels (write out bad channels, manually add them to ecog.bad_chans)
cfg = [];
cfg.viewmode = 'vertical';
cfg.blocksize = 20;
cfg.preproc.bsfilter       = 'yes';
cfg.preproc.bsfiltord      = 3;
cfg.preproc.bsfreq         = [59 61; 119 121; 179 181];
% cfg.preproc.bsfreq         = [59 61; 119 121; 179 181; 200 1000]; %sabina (visualisation purposes only, for noisy data)
cfg.preproc.demean         = 'yes';
good_chns = get_good_chans(ecog,2);
temp = ecog.ftrip;
temp.label = ecog.ftrip.label(good_chns);
temp.trial{1} = ecog.ftrip.trial{1}(good_chns,:);
temp.nChans = length(good_chns);
cfg.channel = temp.label(1:15);
cfg = ft_databrowser(cfg, temp);

%% If you added bad (or SOZ, spikey, out) chans to xls. Read xls in again
% this will overwrite several fields in the ecog structure
% [labelfile,ecog] = read_labelfile(params.labelfile,ecog);
% if you don't want to be asked if you'd want to overwrite.
[labelfile,ecog] = read_labelfile(params.labelfile,ecog,1); 

%% Write the bad channels seen in this block
% ecog.bad_chans = [ecog.bad_chans; {'RPs15';'RFp2';'RFp3';'RFp4';'RHs11';'RHs13';'RPc13'}]; %;'Second';'In Column Order'}];
% ecog.bad_chans = [ecog.bad_chans; {'RPs15';'RFp2';'RFp3';'RFp4';'RHs11';'RHs13';'RPc13'}]; %;'Second';'In Column Order'}];
% ecog.bad_chans = [ecog.bad_chans; {'LLa5';'LIs2';'LDh1'}];

select_bad_chans

save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog.mat']),'ecog','-v7.3');

%% Getting the onsets
% first select edf or tdt analog channel
if isfield(ecog.analog,'trial')
    analog_struct = ecog.analog;
else % for edfs, where 
    analog_struct.trial{1} = ecog.analog.ttl;
    analog_struct.fs = ecog.analog.fs;
end

% If needed: check for threshold
figure; plot(analog_struct.trial{1}(1,:));
title('Check for threshold');

% Analog2digital of the noise channel
if isfield(ecog.analog,'trial') % amplitude threshold
    thr_ampl = 0.01; % amplitude threshold
    noise_ch = analog_struct.trial{1}(1,:);
else % for edfs 
    thr_ampl = 300000;
    noise_ch = demean(analog_struct.trial{1}(1,:));
end
refract_tpts = floor(1*analog_struct.fs); % 1 second only
analog_fs = analog_struct.fs;

digital_trig_chan=analog2digital_trig(noise_ch,thr_ampl,refract_tpts,0);
trial_onsets_tpts=find(digital_trig_chan==1);

% take it from digital channel
trial_onsets_tpts = floor(ecog.digital.onset * ecog.analog.fs);

% plot events
figure('Units','normalized','Position', [0 0  1 .5]);
plot(1/analog_fs:1/analog_fs:length(noise_ch)/analog_fs,noise_ch); hold on
for i=1:length(trial_onsets_tpts)
    plot([trial_onsets_tpts(i)/analog_fs trial_onsets_tpts(i)/analog_fs],ylim)
end
title([num2str(length(trial_onsets_tpts)) ' onsets found']);
xlabel('Time (s)')
xlim([1/analog_fs length(noise_ch)/analog_fs])
print(fullfile(Sbj_Metadata.iEEG_data,curr_block,'PICS','events.jpg'),'-djpeg','-r300')

%% Syncing mic and EEG
tmp = audioinfo(fullfile(Sbj_Metadata.behavioral_root,curr_block,[curr_block '_MicRecording.MP3']));
[MicRec,MicFs]=audioread(fullfile(Sbj_Metadata.behavioral_root,curr_block,[curr_block '_MicRecording.MP3']),[1 tmp.TotalSamples]);
save(fullfile(Sbj_Metadata.behavioral_root,curr_block,[curr_block '_MicRecording.mat']),'MicRec','MicFs')
tmp=load(fullfile(Sbj_Metadata.behavioral_root,curr_block,[curr_block '_MicRecording.mat']));


figure('Units','normalized','Position', [0 0  1 .8]);
subplot(2,1,1)
howlong = 300;%sec
plot((1:(howlong*tmp.MicFs))/tmp.MicFs,tmp.MicRec(1:(howlong*tmp.MicFs),:))
title('mic recording')
subplot(212)
% howlong = 120;%sec
plot((1:(howlong*analog_struct.fs))/analog_struct.fs,ecog.analog.audio(1:(howlong*ecog.analog.fs)))
title('DC channel')

clear MicFs MicRec tmp

%% Event onsets: load behavioral
clear analog_fs noise_ch digital_trig_chan refract_tpts analog_struct thr_ampl beh_data trial_dur curr_stimdur howlong
beh_data = load(fullfile(Sbj_Metadata.behavioral_root,curr_block,[curr_block '.mat']));

% enter_code = beh_data.responses{1,2};
% space_code = beh_data.responses{2,2};
% for i=1:length(beh_data.responses)
%     if beh_data.responses{i,2} == enter_code
%         beh_data.responses{i,1} = 'Enter';
%     elseif beh_data.responses{i,2} == space_code
%         beh_data.responses{i,1} = 'space';
%     else
%         beh_data.responses{i,1} = 'other';
%     end
% end

% Create each event point
trial_onsets = (trial_onsets_tpts/ecog.analog.fs);
trial_durs = zeros(length(trial_onsets_tpts),1);
accuracy = zeros(length(trial_onsets_tpts),1);
% Responses were Enter for absent, space for present words
for t = 1:length(trial_onsets_tpts)
    % set trial duration
    trial_durs(t,1) = length(beh_data.events_table.trials{t})/44100;
    cfg = beh_data.events_table.cfgs{t};
    all_words=[];
    for p=2:6
        curr_part = cfg.(['part' num2str(p)]);
        curr_partfields = fieldnames(curr_part);
        wordfields = curr_partfields(contains(curr_partfields,'word'));
        for w=1:length(wordfields)
            all_words{end+1,1} =  cfg.(['part' num2str(p)]).(wordfields{w}){1};
        end
    end
    % if there was the word and patient responded with space, that is an accurate response
    accuracy(t,1) = any(strcmp(beh_data.words_to_catch{t},all_words)) == strcmp(beh_data.responses(t,1),'space');
end
trial_ends = trial_onsets+trial_durs;
response_time = cell2mat(beh_data.all_times.time_trial_end) - cell2mat(beh_data.all_times.actualStartTime);

% prespeech part lengths
pre_silence_length = 0.5; % in seconds
catchword_screentime = 1.25; % for NS164 and NS165, when catch word screen shows, stays 1sec and cross comes and stays for 0.25 sec before trial starts
speech_onsets = trial_onsets+repmat(pre_silence_length,size(trial_onsets));
catchword_screen = trial_onsets-repmat(catchword_screentime,size(trial_onsets));

events_table = beh_data.events_table;
% Collect trial events
event_ids = (1:size(beh_data.events_table,1))';
sbjID(1:size(beh_data.events_table,1),1) = {sbj_ID};
block(1:size(beh_data.events_table,1),1) = {curr_block};
events = beh_data.events_table;
events = [array2table(sbjID),array2table(block),array2table(event_ids),events,array2table(catchword_screen),array2table(trial_onsets),...
    array2table(speech_onsets),array2table(trial_ends),array2table(trial_durs),array2table(response_time),array2table(accuracy)];

clear block sbjID event_ids trial_onsets att_sent_onset speech_onsets trial_ends trial_durs accuracy response_time tmp...
    c t catchword_screen catchword_screentime cfg w p all_words wordfields curr_part curr_partfields trial_onsets_tpts tmp_events i pre_silence_length

%% Check and correct delays
% ttl_chan = ecog.analog.ttl;
ttl_chan = [];
% audio_chan = ecog.analog.audio;
audio_chan = ecog.analog.trial{1}(1,:);
ttlaudio_smplRate = ecog.analog.fs;

current_onsets = events.trial_onsets;
sound_files = [events.trials, repmat({44100},length(events.trials),1)];for t=1:height(sound_files),sound_files{t,1}=sound_files{t,1}(:,1);end

% new way with cross corr
ecog.analog.time = (1:length(ecog.analog.trial{1}))/ecog.analog.fs;
cor_delays=-1.2:0.001:1.2;
delays=zeros(length(current_onsets),1);
for t = 1:length(current_onsets)
    curr_rec_audio = audio_chan(nearest(ecog.analog.time,current_onsets(t)-2):nearest(ecog.analog.time,current_onsets(t)+events.trial_durs(t)+2));
    curr_soundF = resample(sound_files{t,1},round(ecog.analog.fs),sound_files{t,2});
    corr_res=[];
    for ttt = 1:length(cor_delays)
        tt=cor_delays(ttt);
        filledsound = zeros(length(curr_rec_audio),1);
        filledsound(round((2+tt)*ecog.analog.fs):round((2+tt)*ecog.analog.fs)+length(curr_soundF)-1) = curr_soundF;
        corr_res(ttt,1) = corr(filledsound,curr_rec_audio');
    end
    delays(t,1) = cor_delays(max(corr_res)==corr_res);
end
current_onsets_corr = current_onsets+delays;
% current_onsets_corr = info.events.trial_onsets;

% run the function
[delay_table] = find_delays(ttl_chan/10000000,audio_chan,ttlaudio_smplRate,current_onsets_corr,sound_files);

% add the changes
events.trial_onsets = current_onsets_corr+delay_table(:,3);
events.trial_ends = events.trial_onsets+events.trial_durs;



%% trial based rejection
fs = ecog.ftrip.fsample;
% Notch filter, demean
cfg=[];
cfg.demean         = 'yes';
cfg.bsfilter       = 'yes';
cfg.bsfiltord      = 3;
cfg.bsfreq         = [59 61; 119 121; 179 181];
cont_notched = ft_preprocessing(cfg,ecog.ftrip);

% speech onset locked
pre = 0; % seconds (prespeech part is 3.5455 seconds) [0.5sec + 3.0455sec]
post = 11.5; % seconds (longest trial is ~8.5 seconds)
trl_trg           = [];
trl_trg(:,1)      = floor( events.trial_onsets*fs - fs*pre );
trl_trg(:,2)      = floor( events.trial_onsets*fs + fs*post );
trl_trg(:,3)      = floor( -pre*fs );

% Epoch
cfg      = [];
cfg.trl  = trl_trg;
trials_all  = ft_redefinetrial(cfg,cont_notched);

cfg=[];
cfg.method = 'summary';
cfg.keepchannel = 'yes';
good_chns = get_good_chans(ecog,2);
cfg.channel = good_chns; % assuming bad channels are the same across blocks!!
cfg.metric = 'zvalue';
cfg.keeptrial = 'no';
trials_clean = ft_rejectvisual(cfg,trials_all);

clear trials_clean trials_all trl_trg good_chns fs


%% if there were bad trials, take note here and remove them from further saving
% bad_trials_idx = [3, 5]; % which numbered trials were bad
% all_idx = (1:length(trl_trg))';
% good_trials_idx = all_idx(~ismember(all_idx,bad_trials_idx));
% events = events(good_trials_idx,:);
% 
% clear bad_trials_idx all_idx good_trials_idx

%% Create info.mat file
cfg=[];
cfg.events = events;
cfg.ecog = ecog;
cfg.elecInfo = params.labelfile;
info = create_info(cfg);
save(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']),'info');

% save the auditory responsive electrodes from Auditory Localizer to here too
% [AudRespElecs] = get_AudRespElecs(Sbj_Metadata,project_name);

save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog.mat']),'ecog','-v7.3');

% Re-reference
% Average ref
ecog.ftrip = cont_notched; % nothched or not-noteched
plot_stuff=0;
ignore_szr_chans=1;
ecog_avg=ecog_avg_ref(ecog,plot_stuff,ignore_szr_chans);
save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_avg.mat']),'ecog_avg','-v7.3');

% Also BIPOLAR reference
% ecog_bp=ecog_bipolarize(ecog);
% save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_bp.mat']),'ecog_bp','-v7.3');

%% Wavelet analysis
clearvars -except Sbj_Metadata curr_block AudRespElecs
load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_avg.mat']))
load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']))
events = info.events; clear info

ecog=ecog_avg;
selected_chans = AudRespElecs;
select_bad_chans
clear ecog

% select only auditory responsive electrodes/select electrodes
cfg                = [];
cfg.channel        = selected_chans;
ecog_avg.ftrip     = ft_selectdata(cfg, ecog_avg.ftrip);

% speech onset locked
pre  = 1.5; % seconds
post = 14; % seconds
[epoched_data, epoched_wlt] = master_waveletTF(ecog_avg.ftrip,events.speech_onsets,pre,post,[0.3 200]);

% save wlt and data
save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_wlt.mat']),'epoched_wlt','epoched_data','-v7.3');

%% Prelim analysis
ElecLoc=readtable(fullfile(erase(Sbj_Metadata.data_root,'PROJECTS_DATA'),'DERIVATIVES','freesurfer','ElecLoc_master.xlsx'));
save_folder = fullfile(Sbj_Metadata.results,'iso_4word_sentenceVSscrambled');
if ~exist(save_folder,'dir'),mkdir(save_folder),end

% from fourierspectrum to powerspectrum
cfg = [];
cfg.output='abs';
cfg.keeptrials = 'yes';
epoched_wlt=ft_freqdescriptives(cfg,epoched_wlt);

iso_4_sentence = [];
iso_4_scrambled = [];
for t=1:size(events,1)
    if strcmp(events.cond_info{t}.Sentence_vs_Scrambled,'Sentence') && ...
            strcmp(events.cond_info{t}.Iso_A_chronous_Natural,'iso') && ...
            strcmp(events.cond_info{t}.Word_per_sentence,'4-word')
        iso_4_sentence(1,end+1) = t;
    elseif strcmp(events.cond_info{t}.Sentence_vs_Scrambled,'Scrambled') && ...
            strcmp(events.cond_info{t}.Iso_A_chronous_Natural,'iso') && ...
            strcmp(events.cond_info{t}.Word_per_sentence,'4-word')
        iso_4_scrambled(1,end+1) = t;
    end
end
cfg                = [];
cfg.latency        = [0 8];
cfg.trials         = iso_4_sentence;
iso_4_sentence_wlt = ft_selectdata(cfg, epoched_wlt);
cfg.trials         = iso_4_scrambled;
iso_4_scrambled_wlt= ft_selectdata(cfg, epoched_wlt);

isofrequency = events.cond_info{1}.Frequency;
for el = 1:length(epoched_wlt.label)
    figure('Position',[0 0 900 600])
    curr_label = epoched_wlt.label{el};
    for_avg=squeeze(nanmean(iso_4_sentence_wlt.powspctrm(:,el,:,:),4));
    hold on
    s1=shadedErrorBar(epoched_wlt.freq,mean(for_avg,1),stderr(for_avg),'lineprops','r');
    
    for_avg=squeeze(nanmean(iso_4_scrambled_wlt.powspctrm(:,el,:,:),4));
    hold on
    s2=shadedErrorBar(epoched_wlt.freq,mean(for_avg,1),stderr(for_avg),'lineprops','b');
    xlim([epoched_wlt.freq(1) 6])
    
    legend([s1.mainLine,s2.mainLine],{'Sentence','Scrambled'})
%     if any(stat_1click.mask(el,:,:))
%         to_title1 = '1st-click responsive';
%     else
%         to_title1 = '1st-click unresponsive';
%     end
%     if any(stat_target.mask(el,:,:))
%         to_title2 = 'target responsive';
%     else
%         to_title2 = 'target unresponsive';
%     end
    title({['PSD (mean+SEM); word frequency ' num2str(isofrequency) 'Hz'];...
        ['Elec:' curr_label ' in ' ElecLoc.Area_dtl{strcmp(ElecLoc.Subject,Sbj_Metadata.sbj_ID) & strcmp(ElecLoc.Label,curr_label)}]})
    set(gca, 'FontSize',12);
    ylims = ylim;
%     plot(xlim, [isofrequency isofrequency],'k')
    plot([isofrequency isofrequency], ylim,'k','HandleVisibility','off')
    plot([isofrequency/2 isofrequency/2], ylim,'k','HandleVisibility','off')
    plot([isofrequency/4 isofrequency/4], ylim,'k','HandleVisibility','off')
    text(isofrequency+0.1,ylims(2)-(ylims(2)-ylims(1))/3,'Word-freq','Rotation',90,'FontSize',12)
    text(isofrequency/2+0.1,ylims(2)-(ylims(2)-ylims(1))/3,'Phrase freq','Rotation',90,'FontSize',12)
    text(isofrequency/4+0.1,ylims(2)-(ylims(2)-ylims(1))/3,'Sentence freq','Rotation',90,'FontSize',12)
    
%     text(0,ylims(2)-(ylims(2)-ylims(1))/10,'1st noise period')
    
    set(gca, 'FontSize',13,'FontWeight','bold');
    xlabel('Power');
    ylabel('Frequency (Hz)');
    ylim(ylims)
    
    % Save the figure
    fprintf('\t-Saving electrode #%d %s, out of %d\n',el,curr_label,length(epoched_wlt.label))
    print(fullfile(save_folder,[curr_label , '_' curr_block '_iso_4word.jpg']),'-djpeg','-r300')
    close all
    
end


%% To show the timedomain evolution
for el = 1:length(epoched_wlt.label)
    figure('Position',[0 0 1800 600])
    curr_label = epoched_wlt.label{el};
    subplot(121)
    hold on
    for_avg=squeeze(nanmean(iso_4_sentence_wlt.powspctrm(:,el,:,:),1));
    s1=shadedErrorBar(iso_4_sentence_wlt.time,nanmean(for_avg(nearest(iso_4_sentence_wlt.freq,isofrequency),:),1),stderr(for_avg),'lineprops','r');
    for_avg=squeeze(nanmean(iso_4_sentence_wlt.powspctrm(:,el,:,:),1));
    s2=shadedErrorBar(iso_4_sentence_wlt.time,nanmean(for_avg(nearest(iso_4_sentence_wlt.freq,isofrequency/2),:),1),stderr(for_avg),'lineprops','b');
    for_avg=squeeze(nanmean(iso_4_sentence_wlt.powspctrm(:,el,:,:),1));
    s3=shadedErrorBar(iso_4_sentence_wlt.time,nanmean(for_avg(nearest(iso_4_sentence_wlt.freq,isofrequency/4),:),1),stderr(for_avg),'lineprops','g');
    xlim([iso_4_sentence_wlt.time(1) iso_4_sentence_wlt.time(end)])
    
    subplot(122)
    for_avg=squeeze(nanmean(iso_4_scrambled_wlt.powspctrm(:,el,:,:),1));
    hold on
    s1=shadedErrorBar(iso_4_scrambled_wlt.time,nanmean(for_avg(nearest(iso_4_scrambled_wlt.freq,isofrequency),:),1),stderr(for_avg),'lineprops','r');
    s2=shadedErrorBar(iso_4_scrambled_wlt.time,nanmean(for_avg(nearest(iso_4_scrambled_wlt.freq,isofrequency/2),:),1),stderr(for_avg),'lineprops','b');
    s3=shadedErrorBar(iso_4_scrambled_wlt.time,nanmean(for_avg(nearest(iso_4_scrambled_wlt.freq,isofrequency/4),:),1),stderr(for_avg),'lineprops','g');
    xlim([iso_4_sentence_wlt.time(1) iso_4_sentence_wlt.time(end)])
    
    legend([s1.mainLine,s2.mainLine],{'Sentence','Scrambled'})
%     if any(stat_1click.mask(el,:,:))
%         to_title1 = '1st-click responsive';
%     else
%         to_title1 = '1st-click unresponsive';
%     end
%     if any(stat_target.mask(el,:,:))
%         to_title2 = 'target responsive';
%     else
%         to_title2 = 'target unresponsive';
%     end
    title({['PSD (mean+SEM); word frequency ' num2str(isofrequency) 'Hz'];...
        ['Elec:' curr_label ' in ' ElecLoc.Area_dtl{strcmp(ElecLoc.Subject,Sbj_Metadata.sbj_ID) & strcmp(ElecLoc.Label,curr_label)}]})
    set(gca, 'FontSize',12);
    ylims = ylim;
%     plot(xlim, [isofrequency isofrequency],'k')
    plot([isofrequency isofrequency], ylim,'k','HandleVisibility','off')
    plot([isofrequency/2 isofrequency/2], ylim,'k','HandleVisibility','off')
    plot([isofrequency/4 isofrequency/4], ylim,'k','HandleVisibility','off')
    text(isofrequency+0.1,ylims(2)-(ylims(2)-ylims(1))/3,'Word-freq','Rotation',90,'FontSize',12)
    text(isofrequency/2+0.1,ylims(2)-(ylims(2)-ylims(1))/3,'Phrase freq','Rotation',90,'FontSize',12)
    text(isofrequency/4+0.1,ylims(2)-(ylims(2)-ylims(1))/3,'Sentence freq','Rotation',90,'FontSize',12)
    
%     text(0,ylims(2)-(ylims(2)-ylims(1))/10,'1st noise period')
    
    set(gca, 'FontSize',13,'FontWeight','bold');
    xlabel('Power');
    ylabel('Frequency (Hz)');
    ylim(ylims)
    
    % Save the figure
    fprintf('\t-Saving electrode #%d %s, out of %d\n',el,curr_label,length(epoched_wlt.label))
    print(fullfile(save_folder,[curr_label , '_' curr_block '_iso_4word_timedom.jpg']),'-djpeg','-r300')
    close all
    
end
