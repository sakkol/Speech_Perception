%% Prepare
data_root = '/media/sakkol/HDD1/HBML/';
project_name = 'IsochronousFreeRecall';
sbj_ID = 'NS170';
Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'

% Get params directly from BlockList excel sheet
curr_block = Sbj_Metadata.BlockLists{2}
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

clear cfg temp

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
    thr_ampl = 2000000;
    noise_ch = demean(analog_struct.trial{1}(1,:));
end

ttl_onsets = ecog.digital.PtC4.onset;

% plot events
figure('Units','normalized','Position', [0 0  1 .5]);
plot(1/analog_struct.fs:1/analog_struct.fs:length(noise_ch)/analog_struct.fs,noise_ch); hold on
ylims=ylim;
for i=1:length(ttl_onsets)
    plot([ttl_onsets(i) ttl_onsets(i)],ylim)
    text(ttl_onsets(i),4*ylims(2)/5,num2str(i))
%     plot([trial_onsets_tpts(i) trial_onsets_tpts(i)],ylim)
%     text(trial_onsets_tpts(i),4*ylims(2)/5,num2str(i))
end
title([num2str(length(ttl_onsets)) ' TTL onsets found']);
xlabel('Time (s)')
xlim([1/analog_struct.fs length(noise_ch)/analog_struct.fs])
print(fullfile(Sbj_Metadata.iEEG_data,curr_block,'PICS','_ttls.jpg'),'-djpeg','-r300')

%% Loop over trials to get events and behavioral data
% load beh data
tmp = load(fullfile(Sbj_Metadata.behavioral_root,curr_block,[curr_block '.mat']));

words_info = cell([size(tmp.events_table, 1), 1]);
math_info = cell([size(tmp.events_table, 1), 1]);
SoundStimuli = cell([size(tmp.events_table, 1), 1]);
condition_info = cell([size(tmp.events_table, 1), 1]);
Sound_cfgs = cell([size(tmp.events_table, 1), 1]);
math_onsets = cell([size(tmp.events_table, 1), 1]);
Calc1_RT = zeros(size(tmp.events_table,1), 1);
Calc2_RT = zeros(size(tmp.events_table,1), 1);
FreeRecall1_onset = zeros(size(tmp.events_table,1), 1);
FreeRecall2_onset = zeros(size(tmp.events_table,1), 1);
trial_onsets = zeros(size(tmp.events_table,1), 1);
to = 1;
for t = 1:size(tmp.events_table,1)
    % words in
    for p=2:4
        for w=1:4
            words_info{t,1}{p-1,w} = tmp.events_table.cfgs{t}.(['part' num2str(p)]).(['word' num2str(w)]);
        end
    end
    
    % arrange sound stimuli
    SoundStimuli{t,1} = tmp.events_table.trials{t};
    condition_info{t,1} = tmp.events_table.conds{t};
    Sound_cfgs{t,1} = tmp.events_table.cfgs{t};
    
    % check if calculation was run
    if isempty(tmp.math_responses{t,1})
        trial_onsets(t,1) = ttl_onsets(to);
        FreeRecall1_onset(t,1) = ttl_onsets(to+1);
        to = to + 2;
    else
        trial_onsets(t,1) = ttl_onsets(to);
        FreeRecall1_onset(t,1) = ttl_onsets(to+1);
        FreeRecall2_onset(t,1) = ttl_onsets(to+14);
        
        Calc1_RT(t,1) = ttl_onsets(to+7) - ttl_onsets(to+6);
        Calc2_RT(t,1) = ttl_onsets(to+13) - ttl_onsets(to+12);
        
        % math onsets
        Calc1_Op1 = ttl_onsets(to+2);
        Calc1_plus = ttl_onsets(to+3);
        Calc1_Op2 = ttl_onsets(to+4);
        Calc1_equals = ttl_onsets(to+5);
        Calc1_Op3 = ttl_onsets(to+6);
        Calc2_Op1 = ttl_onsets(to+8);
        Calc2_plus = ttl_onsets(to+9);
        Calc2_Op2 = ttl_onsets(to+10);
        Calc2_equals = ttl_onsets(to+11);
        Calc2_Op3 = ttl_onsets(to+12);
        math_onsets{t,1} = [array2table(Calc1_Op1),array2table(Calc1_plus),array2table(Calc1_Op2),array2table(Calc1_equals),array2table(Calc1_Op3),...
            array2table(Calc2_Op1),array2table(Calc2_plus),array2table(Calc2_Op2),array2table(Calc2_equals),array2table(Calc2_Op3)];

        
        % check if patient was accurate
        Calc1_acc = (tmp.math.accuracy(t,1)==1 && tmp.math_responses{t,1}.whichButton==1) || (tmp.math.accuracy(t,1)==0 && tmp.math_responses{t,1}.whichButton == 3);
        Calc2_acc = (tmp.math.accuracy(t,2)==1 && tmp.math_responses{t,2}.whichButton==1) || (tmp.math.accuracy(t,2)==0 && tmp.math_responses{t,2}.whichButton==3);
        Calc1_first_dig = tmp.math.first_dig(t,1);
        Calc2_first_dig = tmp.math.first_dig(t,2);
        Calc1_stim = tmp.math.math_stimuli(t,1:5);
        Calc2_stim = tmp.math.math_stimuli(t,6:10);
        
        math_info{t,1} = [array2table(Calc1_acc),array2table(Calc2_acc),array2table(Calc1_first_dig),array2table(Calc2_first_dig),cell2table(Calc1_stim),cell2table(Calc2_stim)];
        to = to + 15;
    end
end


trial_ends = trial_onsets+[diff(trial_onsets)-1;60];

% plot event onsets
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



% Collect trial events
events = [cell2table(SoundStimuli),cell2table(condition_info),cell2table(Sound_cfgs),array2table(trial_onsets),cell2table(math_onsets),...
    array2table(Calc1_RT),array2table(Calc2_RT),array2table(FreeRecall1_onset),array2table(FreeRecall2_onset),array2table(trial_ends),cell2table(words_info),cell2table(math_info)];

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
post = 35; % seconds 
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
ecog.ftrip = cont_notched; % nothched or not-notched
plot_stuff=0;
ignore_szr_chans=1;
ecog_avg=ecog_avg_ref(ecog,plot_stuff,ignore_szr_chans);
save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_avg.mat']),'ecog_avg','-v7.3');

% Also BIPOLAR reference
ecog_bp=ecog_bipolarize(ecog);
save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_bp.mat']),'ecog_bp','-v7.3');
