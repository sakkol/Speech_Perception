%% Script to Syllable Perception
% This script combines English and Spanish versions of the task so that
% it'll be much easier to track and save changes.
%
% Serdar Akkol
% Noah Markowitz
% Stephan Bickel
% Human Brain Mapping Lab
% North Shore University Hospital
% September 2020
% Script is based on "Matrix_Speech_Task_masterScript.m"

%% Preliminary Experimental Startup

% Pre-Experiment Items
clear
commandwindow;                                                              % Typed Characters Appear in Command Window

addpath(genpath('Speech_Perception'));
% Create and Locate Directory
if ~exist([pwd filesep 'log'],'dir'); mkdir([pwd filesep 'log']); end
curr_dir=pwd;

% Experimental Setup Settings
par.textcolor = 255;                         % Text appears white
par.textsize = 60;                           % Text font size
par.cross_length = 40;                       % Size in pixels
par.cross_color = 255;
par.N_entrain_stim = 4;                      % Number of entrainment sounds prior to the probe
par.eeg_pPort = 'DFF8';                      % For parallel port
par.port_id = '/dev/ttyACM01';               % Something like '/dev/ttyACM01' for new Ubuntu laptop; or '/dev/tty.usbmodem14101' for Mac; or 'COM4' for Windows
par.rec_comp_mic = 1;                        % if wanting to record microphone from laptop
par.time = string(datetime('now'));          % save the date and time
[~, par.ComputerID] = system('hostname');    % save computer ID
par.PTB_fs = 44100;

%% Prompt for main settings
non_acceptable = 1;
dlg_title = 'Enter information';

while non_acceptable
    prompt = {'Run ID',...
        'How many syllables to present, enter multiples of 64',...
        'TTL Pulse (1 = None, 2 = MMB, 3 = Parallel Port)'};
    def = {'B1_NS001','1','1','1'};
    answer = inputdlg(prompt,dlg_title,1,def);
    
    if isempty(answer); errordlg('Exiting task. Rerun script for dialogue', 'Aborting'); return; end
    
    par.runID = answer{1};
    trial_count = str2double(answer{2});
    ttl_sender = str2double(answer{3});
    
    log_dir = fullfile(curr_dir, 'log', par.runID);
    % Put in bits that checks to make sure input is acceptable
    if exist(log_dir,'dir')
        dlg_title = 'Run ID EXIST ALREADY!';
        non_acceptable = 1;
    elseif rem(trial_count,64) ~= 0
        dlg_title = 'Invalid choice for number of trials, multiples of 64!!';
        non_acceptable = 1;
    elseif ~ismember(ttl_sender,[1, 2, 3])
        dlg_title = 'Invalid choice for TTL sender!';
        non_acceptable = 1;
    else
        non_acceptable = 0;
    end
    
end

%% Setup second part
% Create directory structure to store data
if ~exist(log_dir,'dir'),mkdir(log_dir),end                                % General directory
save_filename = fullfile(log_dir, par.runID);                              % The file that will store the end results

% Anonymous functions for TTL
if ttl_sender == 1
    send_ttl = @(pulse_code, ttl_handle) pulse_code*2;%display('No pulse sent');
    port_handle = 0;
elseif ttl_sender == 2  % Set up TTLs for MMB Trigger box
    try
        port_handle = serial(par.port_id);
        set(port_handle,'BaudRate',9600); % 9600 is recommended by Neurospec
        fopen(port_handle);
        send_ttl = @(pulse_code, ttl_handle) fwrite(ttl_handle, pulse_code);
    catch
        error('MMB No working!');
    end
elseif ttl_sender == 3 % Set TTLs for parallel port
    eeg_obj = io64;                           % obj = io32;
    eeg_status = io64(eeg_obj);               % stat = io32(obj);
    port_handle = hex2dec(par.eeg_pPort);
    send_ttl = @(pulse_code, ttl_handle) io64(eeg_obj,ttl_handle, pulse_code);
end

% send a highpitched noise before starting
questdlg('Remove the splitter cable from the computer!', ...
        'and plug splitter into Marantz input!', ...
        'ok','OK','OK');
%hpn = [rand(44000,2)*0.1;zeros(4000,2);rand(4000,2);zeros(4000,2);rand(4000,2);zeros(4000,2);rand(4000,2)];
%obj = audioplayer(hpn,par.PTB_fs);
%playblocking(obj);
questdlg('Now put splitter back in, thanks!', ...
        'Dont forget to connect mic, thanks!', ...
        'ok','OK','OK');

    
% set where the stimuli will be found and the dialogs
adaptation_intro_msg = ['You will see single syllables\nTell us what you see\n\n'...
    'When you are ready, \npress space bar to start sentences'];
threshold_intro_msg = ['Now, you will listen similar sentences,\nembedded in noise\n\n'...
    'Please try to catch the sentence coming after\n'...
    '''now catch these words''\n\n'...
    'When you are ready, \npress space bar to start each sentence'];
real_deal_msg = ['Again, you will listen similar sentences,\nembedded in noise\n\n'...
    'Please try to catch the sentence coming after\n'...
    '''now catch these words''\n\n'...
    'When you are ready, \npress space bar to start each sentence'];
end_msg = 'Thank you for your time!\n\nSaving data...';
what_sentence = 'What was the syllable?';
pass_list_msg = 'Catch: ';
pass_list_intro_msg = ['You will listen series of words\n\n'...
    'Please try to catch the word\n'...
    'in the beginning of each series\n'...
    'Press space bar to start'];
word_catch_msg = 'Was that present?';

%%
    
select_sylls = {'pa','ba','ta','da','fa','va','sa','za'};

all_sylls = repmat(select_sylls,[1,10]);
all_sylls = all_sylls(randperm(length(all_sylls)))';

words_to_catch=cell(size(all_sylls,1),1);
responses = cell(size(all_sylls,1),1);
startTime = cell(size(all_sylls,1),1);
actualStartTime = cell(size(all_sylls,1),1);
estStopTime = cell(size(all_sylls,1),1);
time_trial_end = cell(size(all_sylls,1),1);


save([save_filename '_tmp0.mat'], 'par', 'all_sylls')

startscreen_now

% Play the sounds
for trialN = 1:size(events_table,1)
    if trialN==1,WaitSecs(.5);end % in the first, give a little pause, it may be very scary at first
    
    % Choose dialog messages
    if trialN==1
        % at the start of the adaptation part
        DrawFormattedText(window, adaptation_intro_msg,'center','center',par.textcolor);
        Screen('Flip',window);
        [~, key, ~] = KbWait(-1);
        if strcmp(KbName(key), 'ESCAPE'); break; end
        % put cross hair
        Screen('DrawLines', window, cross_Coords,CrossWidth, par.cross_color, [xCenter yCenter]);
        Screen('DrawText',window,num2str(trialN),winRect(3)/20,winRect(4)*0.9,par.textcolor);
        Screen('Flip',window);
    end
    
    DrawFormattedText(window, all_sylls{trialN},'center','center',par.textcolor);
    send_ttl(255, port_handle);
    Screen('Flip',window);
    send_ttl(0, port_handle);
    WaitSecs(2); 
    
    % this time period when patient answers, then wait for input to start next trial
    [time_trial_end{trialN,1}, key, ~] = KbWait(-1);
    if strcmp(KbName(key), 'ESCAPE'); break; end
    % draw cross hair
    Screen('DrawLines', window, cross_Coords,CrossWidth, par.cross_color, [xCenter yCenter]);
    Screen('DrawText',window,num2str(trialN+1),winRect(3)/20,winRect(4)*0.9,par.textcolor);
    Screen('Flip',window);
    
    % wait half second before starting next trial
    WaitSecs(0.5);
    
    % save what is in here every 5 trials
    if mod(trialN,10)==0
        save([save_filename 'tmp' num2str(trialN/5) '.mat'], 'par', 'all_sylls')
    end
    
end

% Draw thank you msg
DrawFormattedText(window, end_msg,'center','center',par.textcolor);
Screen('Flip',window);
WaitSecs(0.5);

% gather the times
all_times=[];
all_times.startTime = startTime;startTime=[];
all_times.estStopTime = estStopTime;estStopTime=[];
all_times.actualStartTime=actualStartTime;actualStartTime=[];
all_times.time_trial_end = time_trial_end;time_trial_end=[];

% Save the info
save([save_filename '.mat'], 'par', 'all_sylls','all_times')


