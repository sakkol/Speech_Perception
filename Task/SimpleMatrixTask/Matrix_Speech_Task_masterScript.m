%% Script to run Matrix Sentence based Speech in Noise Task (Matrix Speech Task, in short)
% This script combines English and Spanish versions of the task so that
% it'll be much easier to track and save changes.
%
% Serdar Akkol
% Noah Markowitz
% Stephan Bickel
% Human Brain Mapping Lab
% North Shore University Hospital
% September 2019
% Updates:
% July 2020: added mic recording option and some cleaning.
% July 2020: English and Spanish version combined into one script. Part to set volume level done.

%% Preliminary Experimental Startup

% Pre-Experiment Items
clear all
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
par.port_id = '/dev/ttyACM0';               % Something like '/dev/ttyACM01' for new Ubuntu laptop; or '/dev/tty.usbmodem14101' for Mac; or 'COM4' for Windows
par.rec_comp_mic = 0;                        % if wanting to record microphone from laptop
par.time = string(datetime('now'));          % save the date and time
[~, par.ComputerID] = system('hostname');    % save computer ID
par.PTB_fs = 44100;

%% Prompt for main settings
non_acceptable = 1;
dlg_title = 'Enter information';

while non_acceptable
    prompt = {'Run ID',...
        'English (1) vs Spanish (2)',...
        'Slow (1) vs Fast (2)',...
        'TTL Pulse (1 = None, 2 = MMB, 3 = Parallel Port)',...
        'Adaptation+Threshold (1) vs Real Deal (2) vs No-noise Stimuli (3)'};
    def = {'B1_NS001','1', '1', '1','1'};
    answer = inputdlg(prompt,dlg_title,1,def);
    
    if isempty(answer); errordlg('Exiting task. Rerun script for dialogue', 'Aborting'); return; end
    
    par.runID = answer{1};
    EngvsSpa = str2double(answer{2});
    slowVSfast = str2double(answer{3});
    ttl_sender = str2double(answer{4});
    thresVSreal = str2double(answer{5});
    
    log_dir = fullfile(curr_dir, 'log', par.runID);
    % Put in bits that checks to make sure input is acceptable
    if exist(log_dir,'dir')
        dlg_title = 'Run ID exists already!';
        non_acceptable = 1;
    elseif ~ismember(EngvsSpa,[1, 2])
        dlg_title = 'Invalid choice for English vs Spanish!';
        non_acceptable = 1;
    elseif ~ismember(slowVSfast,[1, 2])
        dlg_title = 'Invalid choice for slow or fast for threshold!';
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
save_filename = [log_dir filesep par.runID];                               % The file that will store the end results

% set where the stimuli will be found and the dialogs
if EngvsSpa == 1
    main_stim_loc = fullfile(curr_dir,'English_main_stim_loc');
    adaptation_intro_msg = ['You will be listening simple sentences,\ncomposed of 5 words\n\n'...
        'When you are ready, \npress space bar to start each sentence'];
    threshold_intro_msg = ['Now, you will listen similar sentences,\nembedded in noise\n\n'...
        'Please try to catch the sentence coming after\n'...
        '''Please pay attention and remember this sentence''\n\n'...
        'When you are ready, \npress space bar to start each sentence'];
    real_deal_msg = ['Now, you will listen similar sentences.\n\n'...
        'Please try to attend and tell us what is the sentence!\n\n'...
        'To start a new sentence, press space bar.'];
    end_msg = 'Thank you for your time!\n\nSaving data...';
    what_sentence = 'What was the sentence?';
elseif EngvsSpa == 2
    main_stim_loc = fullfile(curr_dir,'Spanish_main_stim_loc');
    adaptation_intro_msg = ['Escucharás oraciones simples,\ncompuestas de 5 palabras\n\n'...
        'Cuando esté listo,\npresione la barra espaciadora\npara comenzar cada oración'];
    threshold_intro_msg = ['Ahora, escuchará oraciones similares,\nincrustadas en el ruido\n\n'...
        'Por favor, trata de entender la frase que viene después\n'...
        '''Presta atención y recuerda esta oración.''\n\n'...
        'Cuando esté listo,\npresione la barra espaciadora\npara comenzar cada oración'];
    real_deal_msg = ['Ahora, escuchará oraciones similares.\n\n'...
        'Por favor, intente asistir y díganos cuál es la oración.\n\n'...
        'Para comenzar una nueva oración, presione la barra espaciadora.'];
    end_msg = '¡Gracias por tu tiempo!\n\nGuardando datos...';
    what_sentence = '¿Cuál fue la oración?';
end

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

%% Setup Psychtoolbox and Screen
try

if thresVSreal == 1

% Load all the stimuli for the 'Adaptation Phase' from directory
% called 'adaptation_sounds'
if slowVSfast == 1
    adaption_sounds_dir = fullfile(main_stim_loc, 'adaptation_sentences_slow');
elseif slowVSfast == 2
    adaption_sounds_dir = fullfile(main_stim_loc, 'adaptation_sounds_fast');
end
list_adaptation_sound_files = dir([adaption_sounds_dir filesep '*.wav']);
adaptation_sounds = cell(length(list_adaptation_sound_files),2);
for ss = 1:length(adaptation_sounds)
    adaptation_sounds{ss,1} = list_adaptation_sound_files(ss).name;
    [y, ~] = audioread(fullfile(adaption_sounds_dir, list_adaptation_sound_files(ss).name));
    adaptation_sounds{ss,2} = [y,y]'; % Because it is 1 column
    if par.PTB_fs ~= 24000 % 24000 is what audio is in
        adaptation_sounds{ss,2} = resample(adaptation_sounds{ss,2}',par.PTB_fs,24000)';
    end
end
% Load all the stimuli for the 'Threshold Phase' from directory called
% 'threshold_sounds' 
if slowVSfast == 1
    threshold_sounds_dir = fullfile(main_stim_loc, 'threshold_sounds_slow');
elseif slowVSfast == 2
    threshold_sounds_dir = fullfile(main_stim_loc, 'threshold_sounds_fast');
end
list_threshold_sound_files = dir([threshold_sounds_dir filesep '*.wav']);
threshold_sounds = cell(length(list_threshold_sound_files),2);
for ss = 1:length(threshold_sounds)
    threshold_sounds{ss,1} = list_threshold_sound_files(ss).name;
    [y, ~] = audioread(fullfile(threshold_sounds_dir, list_threshold_sound_files(ss).name));
    threshold_sounds{ss,2} = y';
    if par.PTB_fs ~= 24000 % 24000 is what audio is in
        threshold_sounds{ss,2} = resample(threshold_sounds{ss,2}',par.PTB_fs,24000)';
    end
end

%% SOUND VOLUME SETTING
sound_good=0;
while sound_good~=1
    obj = audioplayer(threshold_sounds{ss,2},par.PTB_fs);
    playblocking(obj);
    sound_good = input('Is the sound level good (1) or bad (2)?');
end

%% Adaptation
% There'll be several example sounds without noise. Sounds are ready to
% play in the adaptation_sounds folder. There will be 2 lines saying:
% 1. Intro:
%   a. 'You will be listening simple sentences composed of 5 words'
%   b. 'When you are ready, press space bar to start each sentence'
% 2. Loop for 30 stimuli (in a cell structure with same order for each
% patient)

startscreen_now

% put dialog message to start adaptation part
DrawFormattedText(window, adaptation_intro_msg,'center','center',par.textcolor);
Screen('Flip',window);
KbWait(-1);

% put cross hair
Screen('DrawLines', window, cross_Coords,CrossWidth, par.cross_color, [xCenter yCenter]);
Screen('Flip',window);

% Play the sounds
for trialN = [1,4,7,10,13,16,19,22,25,28,30,11,14,17,20,23,26,2,6,29,27,24,21,18,15,12,9,3,8,5] % 1:length(adaptation_sounds)
    if trialN==1,WaitSecs(.5);end % in the first, give a little pause, it may be very scary at first
    PsychPortAudio('FillBuffer',pahandle, adaptation_sounds{trialN,2});
    startTime{trialN,1} = PsychPortAudio('Start', pahandle, repetitions, startCue, waitforDeviceStart);
    send_ttl(255, port_handle);
    [actualStartTime{trialN,1}, ~,~,estStopTime{trialN,1}] = PsychPortAudio('Stop',pahandle,1,1);
    send_ttl(0, port_handle); 
    WaitSecs(0.5); % Wait until audio ends and then wait another 0.5 sec
    [time_trial_end{trialN,1}, key, ~] = KbWait(-1);
    if strcmp(KbName(key), 'ESCAPE'); break; end
end
all_times_adapt=[];
all_times_adapt.startTime = startTime;startTime=[];
all_times_adapt.estStopTime = estStopTime;estStopTime=[];
all_times_adapt.actualStartTime=actualStartTime;actualStartTime=[];
all_times_adapt.time_trial_end = time_trial_end;time_trial_end=[];

%% Threshold:
% There are set group of sentences to play in threshold_sounds folder.
% However this time, these sounds will be embedded in noise of 3 different
% levels. 

% put the dialog message
DrawFormattedText(window, threshold_intro_msg,'center','center',par.textcolor);
Screen('Flip',window);
KbWait(-1);

% draw cross hair
Screen('DrawLines', window, cross_Coords,CrossWidth, par.cross_color, [xCenter yCenter]);
Screen('Flip',window);

% Looping for threshold audio
for trialN = 1:length(threshold_sounds)
    PsychPortAudio('FillBuffer',pahandle, threshold_sounds{trialN,2});
    startTime{trialN,1} = PsychPortAudio('Start', pahandle, repetitions, startCue, waitforDeviceStart);
    send_ttl(255, port_handle);
    [actualStartTime{trialN,1}, ~,~,estStopTime{trialN,1}] = PsychPortAudio('Stop',pahandle,1,1);
    send_ttl(0, port_handle);
    if par.rec_comp_mic
        [recordedaudio{trialN,1}, freq{trialN,1}] = audio_capture(pahandle);
    else
        recordedaudio=[];freq=[];
    end
    WaitSecs(0.25); % Wait until audio and/or utterance ends and then wait another 0.25sec
    [time_trial_end{trialN,1}, key, ~] = KbWait(-1);
    if strcmp(KbName(key), 'ESCAPE'); break; end
end

% Draw thank you msg
DrawFormattedText(window, end_msg,'center','center',par.textcolor);
Screen('Flip',window);
WaitSecs(2);

all_times_thresh=[];
all_times_thresh.startTime = startTime;
all_times_thresh.estStopTime = estStopTime;
all_times_thresh.actualStartTime=actualStartTime;
all_times_thresh.time_trial_end = time_trial_end;
save([save_filename '.mat'], 'par', 'EngvsSpa', 'threshold_sounds', 'slowVSfast', 'recordedaudio', 'freq','all_times_adapt','all_times_thresh')


%% Do Real deal
elseif thresVSreal == 2 || thresVSreal == 3
% Creating events: description of output:
% first column is Code of sentence; second is Sentence itself; third and
% fourth are Condition code and name; fifth is stimulus which will be given
% and sixth is cfg input to stim_creatorv2 (saving this to be on the safer side.)
if thresVSreal == 2
    events_cell = event_creator(main_stim_loc,slowVSfast,[]);
elseif thresVSreal == 3
    events_cell = nonoise_event_creator(main_stim_loc,slowVSfast);
end
save([save_filename '0.mat'], 'par','events_cell')

%% Loop for events_cell (=trials)
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
% 5. After stim ends, 'What was the sentence?' screen comes. If no mic
% recording needed, it stays for 2 seconds. If mic recording, until
% response is finished.
% 6. Cross hair comes. End of loop.
% 7. After loop ends, 'Thank you for your time!'

startscreen_now;

DrawFormattedText(window, real_deal_msg,'center','center',par.textcolor);
Screen('Flip',window);
[~, key, ~] = KbWait(-1);
if strcmp(KbName(key), 'ESCAPE'); return; end
    
% Start loop through events_cell
for trialN = 1:length(events_cell)
    % Set sound and code for this trial
    this_trial_sounds = events_cell{trialN,5};
    this_trial_code = events_cell{trialN,1};
    if par.PTB_fs ~= 24000 % 24000 is what audio is in
        this_trial_sounds = resample(this_trial_sounds,par.PTB_fs,24000);
    end
    
    % Create cross and write code to screen for this trial
    Screen('DrawLines', window, cross_Coords,CrossWidth, par.cross_color, [xCenter yCenter]);
    Screen('DrawText',window,num2str(this_trial_code),winRect(3)/20,winRect(4)*0.9,par.textcolor);
    Screen('Flip', window);
    
    % Play sound for this trial
    PsychPortAudio('FillBuffer',pahandle, this_trial_sounds');
    %WaitSecs(0.5);
    startTime{trialN,1} = PsychPortAudio('Start', pahandle, repetitions, startCue, waitforDeviceStart);
    send_ttl(255, port_handle);
    [actualStartTime{trialN,1}, ~,~,estStopTime{trialN,1}] = PsychPortAudio('Stop',pahandle,1,1);
    send_ttl(0, port_handle);
    WaitSecs(0.25); % Wait until audio ends and then wait another 0.25sec
    
    % Prompt for answer from participant. Show message for 2sec
    DrawFormattedText(window, what_sentence,'center','center',par.textcolor);
    Screen('DrawText',window,num2str(this_trial_code),winRect(3)/20,winRect(4)*0.9,par.textcolor);
    Screen('Flip',window);
    if par.rec_comp_mic
        [recordedaudio{trialN,1}, freq{trialN,1}] = audio_capture(pahandle);
    else
        recordedaudio=[];freq=[];
        WaitSecs(1.75);
    end
    WaitSecs(0.25);
    
    % Go back to cross and wait for button press before proceeding to the next trial
    Screen('DrawLines', window, cross_Coords,CrossWidth, par.cross_color, [xCenter yCenter]);
    Screen('DrawText',window,num2str(this_trial_code),winRect(3)/20,winRect(4)*0.9,par.textcolor);
    Screen('Flip', window);
    [time_trial_end{trialN,1}, key, ~] = KbWait(-1);
    if strcmp(KbName(key), 'ESCAPE'); break; end
end

all_times_real=[];
all_times_real.startTime = startTime;
all_times_real.estStopTime = estStopTime;
all_times_real.actualStartTime=actualStartTime;
all_times_real.time_trial_end = time_trial_end;

% Draw thank you msg
DrawFormattedText(window, end_msg,'center','center',par.textcolor);
Screen('Flip',window);
WaitSecs(2);
% Save the info
save([save_filename '.mat'], 'par','events_cell','recordedaudio','freq','all_times_real')

end
%% End of task

% If MMB trigger box was used, close it
if ttl_sender == 2
    fclose(port_handle);
    delete(port_handle);
    clear port_handle
end

% Close Experiment
ListenChar(0);                                                              % Characters Show in Command Window
ShowCursor();                                                               % Shows Cursor
PsychPortAudio('Close', pahandle);                                          % Close the audio device
Screen('CloseAll');                                                         % Close PsychToolbox Screen
sca

catch
    ListenChar(0);                                                          % Characters Show in Command Window
    ShowCursor();                                                           % Shows Cursor
    PsychPortAudio('Close', pahandle);                                      % Close the audio device
    Screen('CloseAll');                                                     % Close PsychToolbox Screen
    sca
    ple                                                                     % Print Last Error
    
    % If MMB trigger box was used, close it
    if ttl_sender == 2
        fclose(port_handle);
        delete(port_handle);
        clear port_handle
    end
end
