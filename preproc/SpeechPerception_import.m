%% Prepare
data_root = '/media/sakkol/HDD1/HBML/';
project_name = 'Speech_Perception';
sbj_ID = 'NS148';
Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'

%% Get params directly from BlockList excel sheet
curr_block = Sbj_Metadata.BlockLists{1}
params = create_Params(Sbj_Metadata,curr_block)


% Import
[ecog] = TDT2ecog(params);

%% Find bad channels using PSD
% find_bad_chans(ecog);
ecog = bad_chan_GUI([],ecog);
if size(ecog.bad_chans,2)>1,ecog.bad_chans = ecog.bad_chans';end % if bad_chans are in row order.

%% Look at EEG, write bad channels to Xls-sheet
cfg=[];
cfg.channel = ecog.bad_chans;
ecog_databrowser(ecog)

% If you added bad (or SOZ, spikey, out) chans to xls. Read xls in again
% this will overwrite several fields in the ecog structure

%% [labelfile,ecog] = read_labelfile(params.labelfile,ecog);
% if you don't want to be asked if you'd want to overwrite.
[labelfile,ecog] = read_labelfile(params.labelfile,ecog,1); 

%% Write the bad channels seen in this block
ecog.bad_chans = [ecog.bad_chans; {'First';'Second';'In Column Order'}];

save(fullfile(params.directoryOUT, [params.filename '_ecog.mat']),'ecog');

%% If needed: check for threshold
% figure; plot(ecog.analog.trial{1}(2,:));
% title('Check for threshold');

% Analog 2 digital of the noise channel
refract_tpts = floor(10*ecog.analog.fs);
thr_ampl = 0.05; % amplitude threshold
noise_ch = ecog.analog.trial{1}(2,:);
digital_trig_chan=analog2digital_trig(noise_ch,0.01,refract_tpts,0);


% Get all events
t=find(digital_trig_chan==1);

pulse_spacing = 1000;

[ev_code, ev_start_tpt]=read_event_codes(digital_trig_chan,pulse_spacing);


% plot events
figure; plot(noise_ch); hold on
for i=1:length(t)
    plot([t(i) t(i)],[-2 2])
end
plot(ev_start_tpt,noise_ch(ev_start_tpt),'*g')
title([num2str(length(ev_code)) ' bursts found, with these burst lengths: ' num2str(unique(ev_code))]);
print(fullfile(Sbj_Metadata.iEEG_data,curr_block,'PICS','events.jpg'),'-djpeg','-r300')

%% Create each event point
% load behavioral
beh_data = load(fullfile(Sbj_Metadata.behavioral_root,curr_block,[curr_block '.mat']));
tmp_table = cell2table(beh_data.events_cell);
tmp_table.Properties.VariableNames(1:end) = {'Sentence_Codes','Sentences','Cond_code','Condition','Stimuli','cfg'};


% trial onset
% convert event time from Accl recording to EEG recording sampling rate
trial_onsets = floor((ev_start_tpt'/ecog.analog.fs)*ecog.ftrip.fsample);

% prespeech end (first if it is slow or fast)
% (this may be needed, because there may be delay before speech starts)
BlockInfo = readtable(fullfile(Sbj_Metadata.project_root,[Sbj_Metadata.project_name '_BlockInfo.xlsx']));
curr_blockinfo = BlockInfo((strcmp(BlockInfo.sbj_ID,Sbj_Metadata.sbj_ID) & strcmp(BlockInfo.BlockList,curr_block)),:);
slowVSfast = BlockInfo.EEGDAT((strcmp(BlockInfo.sbj_ID,Sbj_Metadata.sbj_ID) & strcmp(BlockInfo.BlockList,curr_block)),:);
if strcmp(curr_blockinfo.slowVSfast,'slow')
    prestim_att_length = 73093/24000; % in secs
elseif strcmp(curr_blockinfo.slowVSfast,'fast')
    prestim_att_length = 50848/24000; % in secs
end

prespeech_noise_length = 0.5; % in seconds

prespeech_ends = trial_onsets+floor((prespeech_noise_length+prestim_att_length)*ecog.ftrip.fsample);

% speech onset
% if present, add delays
for i=1:length(prespeech_ends)
    curr_trial = beh_data.events_cell{i,6};
    if isfield(curr_trial,'nostim') % if control condition
        speech_onsets(i,1) = prespeech_ends(i,1);
    else
        speech_onsets(i,1) = prespeech_ends(i,1)+floor(curr_trial.delay);
    end
    
    % stim end
    stim_length(i,1) = length(beh_data.events_cell{i,5})/24000; % in secs
    trial_ends(i,1) = trial_onsets(i,1) + floor(stim_length(i,1)*ecog.ftrip.fsample);
end

%% Collect trial events
events = tmp_table(:,[3,4,1,2]);
events = [events,array2table(trial_onsets),array2table(prespeech_ends),...
    array2table(speech_onsets),array2table(trial_ends),array2table(stim_length)];

ecog.events = events;

%% Create info.mat file

% run if corr sheet doesn't include elecInfo:
% cfg=[];
% cfg.subj_folder = '/media/sakkol/HDD1/HBML/DERIVATIVES/freesurfer/NS148';
% cfg.fsaverage_dir = '/media/sakkol/HDD1/HBML/DERIVATIVES/freesurfer/fsaverage';
% cfg.FS_atlas_info = '/media/sakkol/HDD1/HBML/DERIVATIVES/freesurfer/Freesurfer_Atlas_Labels.xlsx';
% sbj_name = 'NS148';
% create_elecInfo(sbj_name, cfg)

% create info variable:
cfg=[];
cfg.events = events;
cfg.ecog = ecog;
cfg.corr_sheet = '/media/sakkol/HDD1/HBML/PROJECTS_DATA/Speech_Perception/NS148/NS148_Electrodes_Natus_TDT_correspondence.xlsx';
info = create_info(cfg);
save(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']),'info');

%% Re-reference
% Average ref
plot_stuff=0;
ignore_szr_chans=1;
ecog_avg=ecog_avg_ref(ecog,plot_stuff,ignore_szr_chans);
save(fullfile(Sbj_Metadata.iEEG_data, [curr_block '_ecog_avg.mat']),'ecog_avg','-v7.3');

% Also BIPOLAR reference
ecog_bp=ecog_bipolarize(ecog);
save(fullfile(Sbj_Metadata.iEEG_data, [curr_block '_ecog_bp.mat']),'ecog_bp','-v7.3');