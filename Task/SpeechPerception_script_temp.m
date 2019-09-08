%% Temporary Speech Perception Script
clc,clear
%% Folder related
% where the stimuli are: (also an input to event_creator)
main_stim_loc = '/home/sakkol/Documents/Speech_Perception_stim/4th_Generation/';

% Input run ID
runID = inputdlg('Please write the ID of this run!');
runID = runID{1};

% these can be similar to BlinkSuppression/EntrainSounds:
log_dir = [pwd filesep 'SpeechPerception_log'];
if ~exist(lod_dir,'dir'); mkdir(log_dir); end
save_dir = [log_dir filesep runID];
if ~exist(save_dir,'dir'); mkdir(save_dir); end

%% Adaptation
% SA: I'll put several example sounds to play, without noise. Their sounds
% are going to be ready to play. Just couple of introduction may be needed.
% 1. Intro:
%   a. 'You will be listening simple sentences composed of 5 words'
%   b. 'When you are ready, press space bar to start each sentence'
% 2. Loop for 30 stimuli (in a cell structure with same order for each
% patient)


%% Threshold:
% SA: I need to work on this
% Probably same as adaptation, with set group of sentences with 3 different
% SNR levels. If that is, we can use a similar loop as in previous section


%% Real deal
% Creating events: description of output:
% first column is Code of sentence; second is Sentence itself; third and
% fourth are Condition code and name; fifth is stimulus which will be given
% and sixth is cfg input to stim_creatorv2 (saving this to be on the safer side.)
events_cell = event_creator(main_stim_loc);

%% Loop for events_cell
% 1. Intro:
%   a. 'Now, you will listen similar sentences' stays for 3sec
%   b. 'Please try to attend and tell us what was the sentence' stays for 3sec
%   c. 'To start a new sentence, press space bar!' stays for 3sec
% 2. Loop starts:
% 3. Cross hair comes, stays as long as stimulus ends. Also there is going
% to be the code of the stimulus in the lower left corner of the screen.
% This code is the first column of events_cell. It is a note for
% experimenter to check the sheet and find the sentence.
% 4. Button press
% 5. After stim ends, 'What was the sentence?' stays for 2 sentence.
% 6. Cross hair comes. End of loop.
% 7. After loop ends, 'Thank you for your time!'

%% Save
% save events_cell
save([save_dir filesep runID 'events_cell.mat'],'events_cell')

% maybe also something like this:
% par.savetime = datestr(now,'mm.dd.yyyy.HH.MM');
% save([filename], 'par','flip','targets','events','all_stim_order', 'all_trials_events2do','all_trials_events_wait_times', 'event_code_index');

