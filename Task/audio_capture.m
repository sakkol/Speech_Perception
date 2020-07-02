function [recordedaudio, freq] = audio_capture(pahandle)

recordedaudio = [];

% Wait for release of all keys on keyboard:
%KbReleaseWait;

% Get what freq'uency we are actually using:
s = PsychPortAudio('GetStatus', pahandle);
freq = s.SampleRate;

% Preallocate an internal audio recording  buffer with a capacity of 10 seconds:
PsychPortAudio('GetAudioData', pahandle, 10);

% Start audio capture immediately and wait for the capture to start.
% We set the number of 'repetitions' to zero,
% i.e. record until recording is manually stopped.
PsychPortAudio('Start', pahandle, 0, 0, 1);

fprintf('Audio capture started, press any key for about 0.25 second to quit.\n');


% We retrieve status once to get access to SampleRate:
s = PsychPortAudio('GetStatus', pahandle);

% Stay in a little loop until keypress:
while ~KbCheck && ((length(recordedaudio) / s.SampleRate) < maxsecs)
    % Wait 0.25 seconds...
    WaitSecs(0.25);

    % Query current capture status and print it to the Matlab window:
    s = PsychPortAudio('GetStatus', pahandle);

    % Retrieve pending audio data from the drivers internal ringbuffer:
    audiodata = PsychPortAudio('GetAudioData', pahandle);

    % And attach it to our full sound vector:
    recordedaudio = [recordedaudio audiodata];
end

% Stop capture:
PsychPortAudio('Stop', pahandle);

% Perform a last fetch operation to get all remaining data from the capture engine:
audiodata = PsychPortAudio('GetAudioData', pahandle);

% Attach it to our full sound vector:
recordedaudio = [recordedaudio audiodata];

