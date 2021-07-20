%% Prepare
data_root = '/media/sakkol/HDD1/HBML/';
project_name = 'IsochronousFreeRecall';
sbj_ID = 'NS172';
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
analog_fs = analog_struct.fs;
analog_struct.time{1} = 1/analog_fs:1/analog_fs:length(analog_struct.trial{1})/analog_fs;

% % getting events from digitaldata
% trial_onsets_tpts = ecog.digital.PtC4.onset(1:end)' * analog_fs;

% math onsets
Calc1_Op1 = trial_onsets_tpts(2:11:end)';
Calc1_plus = trial_onsets_tpts(3:11:end)';
Calc1_Op2 = trial_onsets_tpts(4:11:end)';
Calc1_equals = trial_onsets_tpts(5:11:end)';
Calc1_Op3 = trial_onsets_tpts(6:11:end)';
Calc2_Op1 = trial_onsets_tpts(7:11:end)';
Calc2_plus = trial_onsets_tpts(8:11:end)';
Calc2_Op2 = trial_onsets_tpts(9:11:end)';
Calc2_equals = trial_onsets_tpts(10:11:end)';
Calc2_Op3 = trial_onsets_tpts(11:11:end)';
FreeRecall_onset_tpts = trial_onsets_tpts(12:11:end)';
math_onsets = [array2table(Calc1_Op1),array2table(Calc1_plus),array2table(Calc1_Op2),array2table(Calc1_equals),array2table(Calc1_Op3),...
               array2table(Calc2_Op1),array2table(Calc2_plus),array2table(Calc2_Op2),array2table(Calc2_equals),array2table(Calc2_Op3)];

% plot events
figure('Units','normalized','Position', [0 0  1 .5]);
plot(1/analog_fs:1/analog_fs:length(noise_ch)/analog_fs,noise_ch); hold on
for i=1:length(trial_onsets_tpts)
    plot([trial_onsets_tpts(i)/analog_fs trial_onsets_tpts(i)/analog_fs],ylim)
end
title([num2str(length(trial_onsets_tpts)) ' onsets found']);
xlabel('Time (s)')
xlim([1/analog_fs length(noise_ch)/analog_fs])
print(fullfile(Sbj_Metadata.iEEG_data,curr_block,'PICS','_all_events.jpg'),'-djpeg','-r300')

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

%% Now create sound events
% If needed: check for threshold
figure; plot(analog_struct.time{1},demean(analog_struct.trial{1}(1,:)));
title('Check for threshold');

% Analog2digital of the noise channel
if isfield(ecog.analog,'trial') % amplitude threshold
    thr_ampl = 0.2; % amplitude threshold
    noise_ch = analog_struct.trial{1}(1,:);
else % for edfs 
    thr_ampl = 2000000;
    noise_ch = demean(analog_struct.trial{1}(1,:));
end
refract_tpts = floor(16*analog_struct.fs); % 0.1 second only

digital_trig_chan=analog2digital_trig(noise_ch,thr_ampl,refract_tpts,0);
trial_onsets_tpts=find(digital_trig_chan==1);

% plot events
figure('Units','normalized','Position', [0 0  1 .5]);
plot(1/analog_fs:1/analog_fs:length(noise_ch)/analog_fs,noise_ch); hold on
for i=1:length(trial_onsets_tpts)
    plot([trial_onsets_tpts(i)/analog_fs trial_onsets_tpts(i)/analog_fs],ylim)
end
title([num2str(length(trial_onsets_tpts)) ' onsets found']);
xlabel('Time (s)')
xlim([1/analog_fs length(noise_ch)/analog_fs])
print(fullfile(Sbj_Metadata.iEEG_data,curr_block,'PICS','_trial_onsets.jpg'),'-djpeg','-r300')

%% Create each event point
clear analog_fs noise_ch digital_trig_chan refract_tpts analog_struct thr_ampl beh_data trial_dur curr_stimdur
% trial onsets and ends
% convert event time from Accl recording to EEG recording sampling rate
% (should be in column)
trial_onsets = (trial_onsets_tpts'/ecog.analog.fs);
FreeRecall_onset = FreeRecall_onset_tpts/analog_fs;
trial_ends = trial_onsets+[diff(trial_onsets);60];

% load beh data
tmp = load(fullfile(Sbj_Metadata.behavioral_root,curr_block,[curr_block '.mat']));

for t = 1:size(tmp.events_table,1)
    % words in 
    for p=2:4
        for w=1:4
        words_info{t,1}{p-1,w} = tmp.events_table.cfgs{t}.(['part' num2str(p)]).(['word' num2str(w)]);
        end
    end
    
    Calc1_RT(t,1) = tmp.all_times.time_trial_end{t,1} - tmp.all_times.actualStartTime{t,1};
    Calc2_RT(t,1) = tmp.all_times.time_trial_end{t,2} - tmp.all_times.actualStartTime{t,1};
    
    % check if patient was accurate
    Calc1_acc(t,1) = (tmp.math.accuracy(t,1)==1 & strcmp(tmp.responses{t,1},'space')) || (tmp.math.accuracy(t,1)==0 & strcmp(tmp.responses{t,1},'Return'));
    Calc2_acc(t,1) = (tmp.math.accuracy(t,2)==1 & strcmp(tmp.responses{t,2},'space')) || (tmp.math.accuracy(t,2)==0 & strcmp(tmp.responses{t,2},'Return'));
    Calc1_first_dig(t,1) = tmp.math.first_dig(t,1);
    Calc2_first_dig(t,2) = tmp.math.first_dig(t,2);
    Calc1_stim{t,1} = tmp.math.math_stimuli(t,1:5);
    Calc2_stim{t,1} = tmp.math.math_stimuli(t,6:10);
    
    math_info{t,1} = [array2table(Calc1_acc),array2table(Calc2_acc),array2table(Calc1_first_dig),array2table(Calc2_first_dig),cell2table(Calc1_stim),cell2table(Calc2_stim)];
    
    % arrange sound stimuli
    SoundStimuli{t,1} = tmp.events_table.trials{t};
    condition_info{t,1} = tmp.events_table.cond_info{t};
    Sound_cfgs{t,1} = tmp.events_table.cfgs{t};
end

% Collect trial events
events = [cell2table(SoundStimuli),cell2table(condition_info),cell2table(Sound_cfgs),array2table(trial_onsets),math_onsets,...
    array2table(Calc1_RT),array2table(Calc2_RT),array2table(FreeRecall_onset),array2table(trial_ends),cell2table(words_info),cell2table(math_info)];

%% Check and correct delays
% ttl_chan = ecog.analog.ttl;
ttl_chan = [];
% audio_chan = ecog.analog.audio;
audio_chan = ecog.analog.trial{1}(1,:);
ttlaudio_smplRate = ecog.analog.fs;

current_onsets = events.trial_onsets;
for i=1:length(events.SoundStimuli)
    xx{i,1} = events.SoundStimuli{i,1}(:,1);
end
sound_files = [xx, repmat({44100},length(events.SoundStimuli),1)];

% new way with cross corr
ecog.analog.time = (1:length(ecog.analog.trial{1}))/ecog.analog.fs;
cor_delays=-1.2:0.001:1.2;
delays=zeros(length(current_onsets),1);
for t = 1:length(current_onsets)
    curr_rec_audio = audio_chan(nearest(ecog.analog.time,current_onsets(t)-2):nearest(ecog.analog.time,current_onsets(t)+(length(sound_files{t})/44100)+2));
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
events.trial_ends = [diff(events.trial_onsets);60];

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
pre = 0; % seconds 
post = 15; % seconds 
trl_trg           = [];
trl_trg(:,1)      = floor( events.trial_onsets*fs - fs*pre );
% trl_trg(:,2)      = floor( events.trial_ends*fs + fs*post );
trl_trg(:,2)      = floor( (events.trial_ends)*fs );
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
% bad_trials_idx = [3, 5];
% all_idx = (1:length(trl_trg))';
% good_trials_idx = all_idx(~ismember(all_idx,bad_trials_idx));
% events = events(good_trials_idx,:);

%% Create info.mat file
cfg=[];
cfg.events = events;
cfg.ecog = ecog;
cfg.elecInfo = ecog.params.labelfile;
info = create_info(cfg);
save(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']),'info');

%% 1.Save raw (unaveraged) ecog 2.average and bipolarize 3.epoch into trials
save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog.mat']),'ecog','-v7.3');

% Re-reference - epoch - save
% Average ref
ecog_avg=ecog;
ecog_avg.ftrip = cont_notched; % nothched or not-notched
plot_stuff=0;
ignore_szr_chans=1;
ecog_avg=ecog_avg_ref(ecog,plot_stuff,ignore_szr_chans);
save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_avg.mat']),'ecog_avg','-v7.3');

% Also BIPOLAR reference
ecog_bp=ecog_bipolarize(ecog);
save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_bp.mat']),'ecog_bp','-v7.3');
