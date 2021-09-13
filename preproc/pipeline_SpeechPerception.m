%% Prepare
data_root = '/media/sakkol/HDD1/HBML/';
project_name = 'Speech_Perception';
sbj_ID = 'NS174_2';
Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'

% Get params directly from BlockList excel sheet
curr_block = Sbj_Metadata.BlockLists{3}
Sbj_Metadata.params = AllBlockInfo2params(Sbj_Metadata,curr_block)

%% if response table hasn't been filled, fill it here
eventcell2responseT

%% Run quick behavioral analysis
SP_beh_analysis(Sbj_Metadata,curr_block)

%% Import
if ~strcmpi(Sbj_Metadata.params.CurrBlockInfo.EEGDAT,'edf')
    [ecog] = TDT2ecog(Sbj_Metadata.params);
    
    % bring in the stimulation output channel and save it in fieldtrip format
    data=TDTbin2mat(Sbj_Metadata.params.directory);
    if isfield(data.streams,'Stim')
        ecog.TDT_stim_output.trial{1} = data.streams.Stim.data;
        ecog.TDT_stim_output.fs = data.streams.Stim.fs;
        ecog.TDT_stim_output.time{1} = (1:length(ecog.TDT_stim_output.trial{1}))/ecog.TDT_stim_output.fs;
    elseif isfield(data.streams,'eS1r')
        ecog.TDT_stim_output.trial{1} = data.streams.eS1r.data;
        ecog.TDT_stim_output.fs = data.streams.eS1r.fs;
        ecog.TDT_stim_output.time{1} = (1:length(ecog.TDT_stim_output.trial{1}))/ecog.TDT_stim_output.fs;
    else
        error('Check stim channel')
    end
    clear data
else
    % for edfs
    Sbj_Metadata.params.edf = fullfile(Sbj_Metadata.params.directory,[Sbj_Metadata.params.filename, '.edf']);
    Sbj_Metadata.params.analog.audioL          = 'DC8';
    Sbj_Metadata.params.analog.audioR          = 'DC2';
    % params.analog.micro          = 'C210';
    ecog=edf2ecog(Sbj_Metadata.params);
%     tmp = ecog.ftrip.time;ecog.ftrip.time={[]};ecog.ftrip.time{1}=tmp;clear tmp
    
    % not needed, just can stay here as a reminder: data = edf2fieldtrip(params.edf)
end

% data=TDTbin2mat(params.directory);
% 
% save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_TDT_data.mat']),'data','-v7.3');

% not needed, just can stay here as a reminder: data = edf2fieldtrip(params.edf)

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
clear temp good_chns cfg

%% If you added bad (or SOZ, spikey, out) chans to xls. Read xls in again
% this will overwrite several fields in the ecog structure
% [labelfile,ecog] = read_labelfile(params.labelfile,ecog);
% if you don't want to be asked if you'd want to overwrite.
[labelfile,ecog] = read_labelfile(params.labelfile,ecog,1); 

%% Write the bad channels seen in this block
% ecog.bad_chans = [ecog.bad_chans; {'RPs15';'RFp2';'RFp3';'RFp4';'RHs11';'RHs13';'RPc13'}];%;'Second';'In Column Order'}];
select_bad_chans

save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog.mat']),'ecog','-v7.3');

%% Finding events
% Event onsets: load behavioral
beh_data = load(fullfile(Sbj_Metadata.behavioral_root,curr_block,[curr_block '.mat']));
tmp_events = cell2table(beh_data.events_cell);
tmp_events.Properties.VariableNames(1:end) = {'Sentence_Codes','Sentences','Cond_code','Condition','Stimuli','trial_details'};

% first select edf or tdt analog channel
if isfield(ecog.analog,'trial')
    analog_struct = ecog.analog;
else % for edfs, where 
    analog_struct.trial{1} = ecog.analog.ttl;
    analog_struct.fs = ecog.analog.fs;
end

if isfield(ecog.analog,'trial') % amplitude threshold
    % Analog2digital of the noise channel
    if strcmpi(tmp_events.trial_details{1}.LvsR, 'L')
        lr=1;
    else
        lr=2;
    end
    thr_ampl = 0.01; % amplitude threshold
    noise_ch = demean(analog_struct.trial{1}(lr,:));
else % for edfs
    thr_ampl = 50000;
    noise_ch = demean(analog_struct.trial{1}(1,:));
end
% If needed: check for threshold
figure; plot(noise_ch);
title('Check for threshold');
refract_tpts = floor(10*analog_struct.fs); % 2 second only
digital_trig_chan=analog2digital_trig(noise_ch,thr_ampl,refract_tpts,0);
analog_fs = analog_struct.fs;

% remove extra digital channel pulses based on sound length
trial_onsets_tpts=find(digital_trig_chan==1);
for i=1:length(tmp_events.Sentence_Codes)
    curr_stimdur = length(tmp_events.Stimuli{i})/24000; % in seconds; events file is 24KHz, recording may be different
    trial_dur = curr_stimdur*floor(ecog.analog.fs); % according to recording time stamps
%     trial_onsets_tpts(trial_onsets_tpts>trial_onsets_tpts(i) & trial_onsets_tpts<trial_onsets_tpts(i)+trial_dur) = [];
    trial_end_tpts(i,1) = trial_onsets_tpts(i) + trial_dur;
    trial_durs(i,1) = trial_dur;
end

% plot events
figure('Units','normalized','Position', [0 0  1 .5]);
plot(1/analog_fs:1/analog_fs:length(noise_ch)/analog_fs,noise_ch); hold on
ylims=ylim;
for i=1:length(trial_onsets_tpts)
    plot([trial_onsets_tpts(i)/analog_fs trial_onsets_tpts(i)/analog_fs],ylim)
    text(trial_onsets_tpts(i)/analog_struct.fs,4*ylims(2)/5,num2str(i))
end
title([num2str(length(trial_onsets_tpts)) ' onsets found']);
xlabel('Time (s)')
xlim([1/analog_fs length(noise_ch)/analog_fs])
print(fullfile(Sbj_Metadata.iEEG_data,curr_block,'PICS','events.jpg'),'-djpeg','-r300')



% plot events
trial_onsets = ecog.digital.PtAB.onset(1:2:end);
trial_ends = trial_onsets+[diff(trial_onsets)-1;20];
trial_durs = trial_ends-trial_onsets;

figure('Units','normalized','Position', [0 0  1 .5]);
plot(1/analog_struct.fs:1/analog_struct.fs:length(noise_ch)/analog_struct.fs,noise_ch); hold on
ylims=ylim;
for i=1:length(trial_onsets)
    plot([trial_onsets(i) trial_onsets(i)],ylim)
    text(trial_onsets(i),4*ylims(2)/5,num2str(i))
%     plot([trial_onsets_tpts(i) trial_onsets_tpts(i)],ylim)
%     text(trial_onsets_tpts(i),4*ylims(2)/5,num2str(i))
end
title([num2str(length(trial_onsets)) ' onsets found']);
xlabel('Time (s)')
xlim([1/analog_struct.fs length(noise_ch)/analog_struct.fs])
print(fullfile(Sbj_Metadata.iEEG_data,curr_block,'PICS','event_onsets.jpg'),'-djpeg','-r300')

%% Create each event point
clear analog_fs noise_ch digital_trig_chan refract_tpts analog_struct thr_ampl beh_data trial_dur curr_stimdur
% trial onsets and ends
% convert event time from Accl recording to EEG recording sampling rate
% (should be in column)
trial_onsets = (trial_onsets_tpts'/ecog.analog.fs);
trial_ends = (trial_end_tpts/ecog.analog.fs);
trial_durs = trial_ends-trial_onsets;

% prespeech part lengths (depending on if it is slow or fast)
% (this may be needed, because there may be delay before speech starts)
prespeech_noise_length = 0.5; % in seconds
curr_blockinfo = Sbj_Metadata.params.CurrBlockInfo;
if strcmp(curr_blockinfo.slowVSfast,'slow')
    prestim_att_length = 73093/24000; % in secs
elseif strcmp(curr_blockinfo.slowVSfast,'fast')
    prestim_att_length = 50848/24000; % in secs
end

att_sent_onset = trial_onsets+repmat(prespeech_noise_length,size(trial_onsets));
speech_onsets = trial_onsets+repmat(prespeech_noise_length+prestim_att_length,size(trial_onsets));

% word,syllable and phoneme onsets and gather Responses
if strcmp(curr_blockinfo.Language, 'English')
    tmp = load('English_all_info.mat');
    all_info_table = tmp.all_info_table;clear tmp
elseif strcmp(curr_blockinfo.Language, 'Spanish')
    tmp = load('Spanish_all_info.mat');
    all_info_table = tmp.peakInfotable;clear tmp
end

response_table = readtable(fullfile(Sbj_Metadata.behavioral_root, curr_block, [curr_block '_response_table.xlsx']));
% convert sub/verb/numb/adj/noun column into string if needed
for c=4:8
    if isnumeric(response_table.(response_table.Properties.VariableNames{c}))
        response_table.(response_table.Properties.VariableNames{c}) = cellstr(num2str(response_table.(response_table.Properties.VariableNames{c})));
    end
end

timing_info=all_info_table;timing_info(:,:)=[];
accuracy = response_table(:,9);
for t = 1:size(tmp_events,1)
    timing_info(t,:) = all_info_table(strcmpi(all_info_table.sentence,tmp_events.Sentences{t}) & ...
        strcmpi(all_info_table.rate,curr_blockinfo.slowVSfast),:);
        
    if any(ismember(timing_info.Properties.VariableNames,'word_info'))
    timing_info.word_info{t}(:,4)=response_table{t,4:8}';
    timing_info.word_info{t}.Properties.VariableNames(4) = {'response'};
    end
end

% Collect trial events
event_ids = (1:size(tmp_events,1))';
sbjID(1:size(tmp_events,1),1) = {sbj_ID};
block(1:size(tmp_events,1),1) = {curr_block};
events = tmp_events(:,[3,4,1,5,6]);
events = [array2table(sbjID),array2table(block),array2table(event_ids),events,array2table(trial_onsets),array2table(att_sent_onset),...
    array2table(speech_onsets),array2table(trial_ends),array2table(trial_durs),...
    timing_info,accuracy];

clear block sbjID event_ids trial_onsets att_sent_onset speech_onsets trial_ends trial_durs accuracy timing_info ...
    c t response_table all_info_table prestim_att_length trial_onsets_tpts tmp_events i prespeech_noise_length

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
pre = 4; % seconds (prespeech part is 3.5455 seconds) [0.5sec + 3.0455sec]
post = 6; % seconds (longest trial is ~8.5 seconds)
trl_trg           = [];
trl_trg(:,1)      = floor( events.speech_onsets*fs - fs*pre );
trl_trg(:,2)      = floor( events.speech_onsets*fs + fs*post );
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
cfg.elecInfo = ecog.params.labelfile;
info = create_info(cfg);
save(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']),'info');

% save the auditory responsive electrodes from Auditory Localizer to here too
[AudRespElecs] = get_AudRespElecs(Sbj_Metadata,project_name);

save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog.mat']),'ecog','-v7.3');

% Re-reference
% Average ref
ecog_avg = ecog;
ecog_avg.ftrip = cont_notched; % nothched or not-noteched
plot_stuff=0;
ignore_szr_chans=1;
ecog_avg=ecog_avg_ref(ecog_avg,plot_stuff,ignore_szr_chans);
save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_avg.mat']),'ecog_avg','-v7.3');

% Also BIPOLAR reference
ecog_bp=ecog_bipolarize(ecog);
save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_bp.mat']),'ecog_bp','-v7.3');

%% Wavelet analysis
clearvars -except Sbj_Metadata curr_block
load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_avg.mat']))
load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']))
events = info.events; clear info

% speech onset locked
pre  = 4; % seconds (prespeech part is 3.5455 seconds) [0.5sec + 3.0455sec] - to not include any NaNs in wavelet increase by 2 seconds (from 4 to 6)
post = 6; % seconds (longest trial is ~8.5 seconds)
[epoched_data, epoched_wlt] = master_waveletTF(ecog_avg.ftrip,events.speech_onsets,pre,post);

% save wlt and data
save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_wlt.mat']),'epoched_wlt','epoched_data','-v7.3');

%% Run the EFields analysis
aroundPeak = 1;
elecsOI = {'allchans','goodchans','selectchans','selectgoodchans'};

BlockInfo = makeBlockInfo(Sbj_Metadata,curr_block);
for e = 1:4
    EA_efields(Sbj_Metadata,curr_block,aroundPeak,{[BlockInfo.StimSide{1} 'omni'],BlockInfo.StimSide{1}},elecsOI{e})
    close all
end

% % save the select channels for EA_efields use
% selectchans = info.channelinfo.Label(contains(info.channelinfo.Label,{'LTp','LDp','LDa','LDh','LTx','LTi','LTc','LTs','LPc'}) & isnan(info.channelinfo.outofthebrain))
% cd(fullfile(Sbj_Metadata.sbjDir));
% save selectchans.mat selectchans