%% Steps to run Speech Perception task

%% Folders related
clc,clear
% Where the stimuli is (the directory of where Sentences_Rate are).
main_stim_loc = '/home/sakkol/Documents/Speech_Perception_stim/4th_Generation/';

% Where you want to save the events and possibly trigger onsets
if ~exist([pwd filesep 'Speech_Perception_log'],'dir'); mkdir([pwd filesep 'Speech_Perception_log']); end
curr_dir=pwd;
log_dir = [pwd filesep 'Speech_Perception_log'];

% Run ID dialog and make dir
run_ID = inputdlg('Please write run_ID!');
if isempty(run_ID); BailOut; return; end
run_log_dir = [log_dir filesep run_ID];
if ~exist(run_log_dir,'dir'),mkdir(run_log_dir);end

%% Adaptation
% here maybe several sentences to get the sound volume?
% 1.Intro: Small intrduction saying 'You will listen simple sentences with
%    5 words, please try to attend'

% SA: I need to work on this still. Maybe just simple presentation similar
% to loop events part. Events_cell structure will be the same.


%% Threshold: different script?
% 1.Intro: Small intrduction saying 
%   a. "Now, you will listen similar sentences embedded in different noise levels'
%   b. 'Please try to focus and catch the sentence which will start after
%   'Please pay attention and remember this sentence'

% SA: I need to work on this still. 


%% Real deal: Create events
events_cell = event_creator(main_stim_loc);

%% loop events
% Event steps: 
% 1. Intro: Small intrduction saying 
%       a. 'Now, you will listen sentences with same noise level' (stays 5 seconds)
%       b. 'To start listening, push space bar?????'
% 2. Loop starts. Events can be presented with same order as in events_cell
% 3. 'Push space bar when ready' or cross hair comes and stops
% 4. Button press
% 5. Stimulus is presented
% 6. 'What was the sentence?' comes and stays for 3 seconds
% 7. After 3 seconds, cross hair comes again
% Going back step 2 until events end
% 8. Last prompt: 'Thank you for your time!'

%% Save
save([run_log_dir filesep run_ID '_events.mat'],'events_cell');
fprintf('/tEvents are saved to %s/n',[run_log_dir filesep run_ID '_events.mat'])

% maybe also something like this
%     par.savetime = datestr(now,'mm.dd.yyyy.HH.MM');
%     save([filename], 'par','flip','targets','events','all_stim_order', 'all_trials_events2do','all_trials_events_wait_times', 'event_code_index');
