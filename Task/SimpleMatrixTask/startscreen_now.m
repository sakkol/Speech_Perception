PsychDefaultSetup(2);
HideCursor;
% Setup Screen and
Screen('Preference', 'SkipSyncTests', 1); % Setup Screen
screenid = max(Screen('Screens')); %Change back to maximum after finished testing; min will make the default the secondary display
[scresw, scresh]=Screen('windowSize',screenid); % Get screen resolution
center = [scresw scresh]/2; % Center of screen

% Setup Cross
xCenter = (scresw)/2;
yCenter = (scresh)/2;
CrossArmLength = par.cross_length;
CrossWidth = 6;
xCoords = [-CrossArmLength CrossArmLength 0 0];                             % Set Center of Cross, x
yCoords = [0 0 -CrossArmLength CrossArmLength];                             % Set Center of Cross, y
cross_Coords = [xCoords; yCoords];

% Retrieves color codes for black and white and grey.
black = BlackIndex(screenid);
white = WhiteIndex(screenid);
grey = (black + white) / 2;
if round(grey)==white; grey=black; end

%%% For debugging make screen small (comment this, use 2nd line for the real deal)
% [window winRect] = Screen('Openwindow', screenid, grey,[0 0 640 480]);
[window winRect] = Screen('Openwindow', screenid, grey);

% Set text parameters
Screen('TextFont',window,'Times');
Screen('TextSize',window,par.textsize);

% Setup PsychPortAudio and Sounds
InitializePsychSound(1);                             % Initialize Sound Driver
nrchannels = 2;                                      % Number of Channels
fs = par.PTB_fs;                                     % Frequency of Sampling
repetitions = 1;                                     % Number of Times to Play Sound Each Time (Keep at 1)
startCue = 0;                                        % Start Immeadiately
pahandle = PsychPortAudio('Open', [], 1, 1, fs, nrchannels); % Open Psych-Audio Port
PsychPortAudio('Volume', pahandle, 1);
waitforDeviceStart = 1;                 % Integrate Device and Sound Driver