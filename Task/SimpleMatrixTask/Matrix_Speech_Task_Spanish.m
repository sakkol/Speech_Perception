%% Task Script of Spanish version of Matrix Sentence Test based Speech Perception Task

% Serdar Akkol
% Noah Markowitz
% Stephan Bickel
% Human Brain Mapping Lab
% North Shore University Hospital
% September 2019
% Updates:
% June 2020: added mic recording option and some cleaning
%% Preliminary Experimental Startup

% Pre-Experiment Items
clear all
commandwindow;                                                              % Typed Characters Appear in Command Window 

addpath(genpath('Speech_Perception'));
% Create and Locate Directory
if ~exist([pwd filesep 'log'],'dir'); mkdir([pwd filesep 'log']); end
curr_dir=pwd;
main_stim_loc = fullfile(curr_dir,'main_stim_loc');

% Experimental Setup Settings
par.textcolor = 255;                         % Text appears white
par.textsize = 60;                           % Text font size
par.cross_length = 40;                       % Size in pixels
par.cross_color = 255;
par.N_entrain_stim = 4;                      % Number of entrainment sounds prior to the probe 
par.eeg_pPort = 'DFF8';                      % For parallel port
par.port_id = '/dev/tty.usbmodem14201';      % Something like '/dev/tty.usbmodem14101' for Mac and something like 'COM4' for Windows
par.rec_comp_mic = 1;                        % if wanting to record microphone from laptop

%% Prompt for settings
non_acceptable = 1;
dlg_title = 'Enter information';
while non_acceptable
    prompt = {'Run ID',...
        'Slow (1) vs Fast (2)',...
    'TTL Pulse (1 = None, 2 = MMB, 3 = Parallel Port)',...
    'Adaptation+Threshold (1) vs Real Deal (2)'};
    def = {'B1_NS001', '1', '1','1'};
    answer = inputdlg(prompt,dlg_title,1,def);
    
    if isempty(answer); errordlg('Exiting task. Rerun script for dialogue', 'Aborting'); return; end
    
    par.runID = answer{1};
    slowVSfast = str2double(answer{2});
    ttl_sender = str2double(answer{3});
    thresVSreal = str2double(answer{4});
    log_dir = [curr_dir filesep 'log' filesep par.runID];
    
    % Put in bits that checks to make sure input is acceptable
    if exist(log_dir,'dir')
        dlg_title = 'Run ID exists already!';
        non_acceptable = 1;
    elseif ~ismember(slowVSfast,[1, 2])
        dlg_title = 'Invalid choice for slow or fast for threshold!';
        non_acceptable = 1;
    elseif ~ismember(ttl_sender,[1, 2,3])
        dlg_title = 'Invalid choice for TTL sender!';
        non_acceptable = 1;
    else
        non_acceptable = 0;
    end
        
end

%% Setup second part

% Create directory structure to store data
mkdir(log_dir);                                                           % General directory
save_filename = [log_dir filesep par.runID];                  % The file that will store the end results
tmp_dir = fullfile(log_dir, 'tmp');
mkdir(tmp_dir);    % Temporary storage

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
        eeg_obj = io64;                                   % obj = io32;
        eeg_status = io64(eeg_obj);               % stat = io32(obj);
        port_handle = hex2dec(par.eeg_pPort); 
        send_ttl = @(pulse_code, ttl_handle) io64(eeg_obj,ttl_handle, pulse_code);
end

%% Setup Psychtoolbox and Screen
try

if thresVSreal == 1

startscreen_now

%% Adaptation
% There'll be several example sounds without noise. Sounds are ready to
% play in the adaptation_sounds folder. There will be 2 lines saying:
% 1. Intro:
%   a. 'You will be listening simple sentences composed of 5 words'
%   b. 'When you are ready, press space bar to start each sentence'
% 2. Loop for 30 stimuli (in a cell structure with same order for each
% patient)


% Load all the stimuli for the 'Adaptation Phase' from directory
% called 'adaptation_sounds'
if slowVSfast == 1
    adaption_sounds_dir = fullfile(main_stim_loc, 'adaptation_sentences_slow');
elseif slowVSfast == 2
    adaption_sounds_dir = fullfile(main_stim_loc, 'adaptation_sentences_fast');
end
list_adaptation_sound_files = dir([adaption_sounds_dir filesep '*.wav']);
adaptation_sounds = cell(length(list_adaptation_sound_files),2);
for ss = 1:length(adaptation_sounds)
    adaptation_sounds{ss,1} = list_adaptation_sound_files(ss).name;
    [y, ~] = audioread(fullfile(adaption_sounds_dir, list_adaptation_sound_files(ss).name));
    adaptation_sounds{ss,2} = [y,y]'; % Because it is 1 column
end

adaptation_intro_msg = ['Escucharás oraciones simples,\ncompuestas de 5 palabras\n\n'...
    'Cuando esté listo,\npresione la barra espaciadora\npara comenzar cada oración'];

DrawFormattedText(window, adaptation_intro_msg,'center','center',par.textcolor);
Screen('Flip',window);
KbWait(-1);

Screen('DrawLines', window, cross_Coords,CrossWidth, par.cross_color, [xCenter yCenter]);
Screen('Flip',window);

% Play the sounds
for trialN = [1,4,7,10,13,16,19,22,25,28,30,11,14,17,20,23,26,2,6,29,27,24,21,18,15,12,9,3,8,5] % 1:length(adaptation_sounds)
    if trialN==1,WaitSecs(.5);end % in the first, give a little pause, it may be very scary at first
    PsychPortAudio('FillBuffer',pahandle, adaptation_sounds{trialN,2});
    startTime = PsychPortAudio('Start', pahandle, repetitions, startCue, waitforDeviceStart);
    send_ttl(255, port_handle);
    [actualStartTime, ~,~,estStopTime] = PsychPortAudio('Stop',pahandle,1,1);
    send_ttl(0, port_handle); 
    WaitSecs(0.5); % Wait until audio ends and then wait another 1sec
    [time_trial_end, key, ~] = KbWait(-1);
    if strcmp(KbName(key), 'ESCAPE'); break; end
end

%% Threshold:
% There are set group of sentences to play in threshold_sounds folder.
% However this time, these sounds will be embedded in noise of 3 different
% levels. 

% Looping for threshold audio
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
end

threshold_intro_msg = ['Ahora, escuchará oraciones similares,\nincrustadas en el ruido\n\n'...
    'Por favor, trata de entender la frase que viene después\n'...
    '''Presta atención y recuerda esta oración.''\n\n'...
    'Cuando esté listo,\npresione la barra espaciadora\npara comenzar cada oración'];

DrawFormattedText(window, threshold_intro_msg,'center','center',par.textcolor);
Screen('Flip',window);
KbWait(-1);

Screen('DrawLines', window, cross_Coords,CrossWidth, par.cross_color, [xCenter yCenter]);
Screen('Flip',window);

% Play the sounds
for trialN = 1:length(threshold_sounds)
    PsychPortAudio('FillBuffer',pahandle, threshold_sounds{trialN,2});
    startTime = PsychPortAudio('Start', pahandle, repetitions, startCue, waitforDeviceStart);
    send_ttl(255, port_handle);
    [actualStartTime, ~,~,estStopTime] = PsychPortAudio('Stop',pahandle,1,1);
    send_ttl(0, port_handle);
    if par.rec_comp_mic
        [recordedaudio{trialN,1}, freq{trialN,1}] = audio_capture(pahandle);
    else
        recordedaudio=[];freq=[];
    end
    WaitSecs(0.25); % Wait until audio and/or utterance ends and then wait another 0.25sec
    [time_trial_end, key, ~] = KbWait(-1);
    if strcmp(KbName(key), 'ESCAPE'); break; end
end

save([save_filename '.mat'], 'par', 'threshold_sounds', 'slowVSfast','recordedaudio','freq')

%% Do the real task
elseif thresVSreal == 2

% Threshold estimation
% Input accuracies for each SNR
if slowVSfast == 1
   threshold_accr = inputdlg({'SNR = -8:','SNR = -4:','SNR = 0:'},'Please input how many words were correct for each SNR!');
data_threshold =    [...
    -8,   str2double(threshold_accr{1}),   50.0000;...
    -4,   str2double(threshold_accr{2}),   50.0000;...
    0,   str2double(threshold_accr{3}),   50.0000];
elseif slowVSfast == 2
    threshold_accr = inputdlg({'SNR = -6:','SNR = -2:','SNR = 2:'},'Please input how many words were correct for each SNR!');
data_threshold =    [...
    -6,   str2double(threshold_accr{1}),   50.0000;...
    -2,   str2double(threshold_accr{2}),   50.0000;...
    2,   str2double(threshold_accr{3}),   50.0000];
end

options             = [];       % initialize as an empty struct
options.sigmoidName = 'norm';   % normal Gaussian curve as the sigmoid
% options.expType     = '4AFC'; % ours is not forced choice, thus no input should be given
options.threshPC       = .5;

result = psignifit(data_threshold,options);
p1=plotPsych(result);
print([save_filename '_psignifit.jpg'],'-djpeg','-r300')
current_SNR = getThreshold(result,0.5)
% close p1

%% Real deal
% Creating events: description of output:
% first column is Code of sentence; second is Sentence itself; third and
% fourth are Condition code and name; fifth is stimulus which will be given
% and sixth is cfg input to stim_creatorv2 (saving this to be on the safer side.)
events_cell = event_creator(main_stim_loc,slowVSfast,current_SNR);

%% Loop for events_cell, meaning trials
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

real_deal_msg = ['Ahora, escuchará oraciones similares.\n\n'...
   'Por favor, intente asistir y díganos cuál es la oración.\n\n'...
   'Para comenzar una nueva oración, presione la barra espaciadora.'];

DrawFormattedText(window, real_deal_msg,'center','center',par.textcolor);
Screen('Flip',window);
[~, key, ~] = KbWait(-1);
if strcmp(KbName(key), 'ESCAPE'); return; end
    
% Start loop through events_cell
for trialN = 1:length(events_cell)
    % Set sound and code for this trial
    this_trial_sounds = events_cell{trialN,5};
    this_trial_code = events_cell{trialN,1};
    
    % Create cross and write code to screen for this trial
    Screen('DrawLines', window, cross_Coords,CrossWidth, par.cross_color, [xCenter yCenter]);
    Screen('DrawText',window,num2str(this_trial_code),winRect(3)/20,winRect(4)*0.9,par.textcolor);
    Screen('Flip', window);
    
    % Play sound for this trial
    PsychPortAudio('FillBuffer',pahandle, this_trial_sounds');
    %WaitSecs(0.5);
    startTime = PsychPortAudio('Start', pahandle, repetitions, startCue, waitforDeviceStart);
    send_ttl(255, port_handle);
    [actualStartTime, ~,~,estStopTime] = PsychPortAudio('Stop',pahandle,1,1);
    send_ttl(0, port_handle);
    WaitSecs(0.25); % Wait until audio ends and then wait another 0.25sec
    
    % Prompt for answer from participant. Show message for 2sec
    DrawFormattedText(window, '¿Cuál fue la oración?','center','center',par.textcolor);
    Screen('DrawText',window,num2str(this_trial_code),winRect(3)/20,winRect(4)*0.9,par.textcolor);
    Screen('Flip',window);
    if par.rec_comp_mic
        [recordedaudio{trialN,1}, freq{trialN,1}] = audio_capture(pahandle);
    else
        recordedaudio=[];freq=[];
        WaitSecs(1.75);
    end
    WaitSecs(0.25);
    
    % Go back to cross and wait for button press before proceeding to
    % next trial
    Screen('DrawLines', window, cross_Coords,CrossWidth, par.cross_color, [xCenter yCenter]);
    Screen('DrawText',window,num2str(this_trial_code),winRect(3)/20,winRect(4)*0.9,par.textcolor);
    Screen('Flip', window);
    [time_trial_end, key, ~] = KbWait(-1);
    if strcmp(KbName(key), 'ESCAPE'); break; end
end

% Save the info
save([save_filename '.mat'], 'par','result','options','events_cell','data_threshold','recordedaudio','freq')

end
%% End of task

% Draw thank you msg
end_msg = '¡Gracias por tu tiempo!';
DrawFormattedText(window, end_msg,'center','center',par.textcolor);
Screen('Flip',window);
WaitSecs(2);

% If MMB trigger box was used, close it
if ttl_sender == 2
    fclose(port_handle);
    delete(port_handle);
    clear port_handle
end

% Close Experiment
ListenChar(0);                                                               % Characters Show in Command Window
ShowCursor();                                                               % Shows Cursor
PsychPortAudio('Close', pahandle);                                 % Close the audio device
Screen('CloseAll');                                                         % Close PsychToolbox Screen
sca

catch
    
    ListenChar(0);                                                          % Characters Show in Command Window
    ShowCursor();                                                           % Shows Cursor
    %PsychPortAudio('Close', pahandle);                             % Close the audio device
    Screen('CloseAll');                                                     % Close PsychToolbox Screen
    sca
    ple                                                                     % Print Last Error
end