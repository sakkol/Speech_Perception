function [stimulus,envelope]=stim_creatorv2(cfg)
% In order to be able to make it very malleable, there is going to be some
% subfields in cfg. Here are the possible subfields:
% cfg.SNR = eg. -3.5, 2. Signal-to-noise ratio, that is going to be applied
%           in only speech+noise part. (default=-3)
% cfg.prespeech.part1.noise
%                    .length
%              .part2.signal
%                    .noise
%                    .length
%                   each part field can include .length (in s), .signal and/or .noise
%                   (combination of {'pink','white','silence'} or specific 
%                   file(s) that is/are going to be used). So that there
%                   can be max 2 noises combined, one as signal other as
%                   noise, using the SNR above. If there's length field
%                   with signal, that part will be in that length.
% cfg.speech.file : speech file to be used. (also signal)
%           .noise: same options as in prespeech (length and file).
% cfg.postspeech  : same options as in prespeech (parts, length and file).
% cfg.LvsR        : 'L','R' (defaulft: error).
% cfg.stim_save_filename    
%    .plot_save_filename
% 	 .envelope_save_filename
%                 : if there is this field, stimulus/plot will be saved with
%                 this name. Plot will include 1.normalized power in lower
%                 frequency range (1-8Hz), 2.peaks in this normalized power
%                 plot, 3.envelope of speech signal embedded in whole
%                 stimulus.
% cfg.delay       : when the speech envelope should start. If number (in
%                 seconds), stimulation will be starting with this delay
%                 relative to matrix sentence onset. If 'attention',
%                 stimulation will start with Attention sentence onset. In
%                 this condition, there should be 2 other inputs:
%                 'frequency' in which the wave will be generated; 'phase',
%                 to generate 'in'-phase or 'off'-phase.
%
% OUTPUT:
% stimulus: one channel is for presenting to patient and the other channel
% goes into stimulator, composed of only clean speech signal.
% envelope: envelope of signal (speech) in the stimuli. If there is noise
% before and/or after signal, those parts are going to stay zero.
% 
% Example:
% cfg=[];
% cfg.SNR = -4;
% cfg.prespeech.part1.length= 0.5;
% cfg.prespeech.part1.noise = 'pink';
% cfg.prespeech.part2.noise = 'pink';
% main_stim_loc = '/home/sakkol/Documents/Speech_Perception_stim/4th_Generation/';
% cfg.prespeech.part2.signal = find_sentence('Pre-stim-Attention-comma',main_stim_loc,0.9);
% cfg.speech.noise = 'pink';
% cfg.speech.file = find_sentence(curr_sentence,main_stim_loc,speech_rate);
% cfg.postspeech.part1.length=2;
% cfg.postspeech.part1.noise = 'pink';
% cfg.LvsR = 'L';
% cfg.delay = 0.05;
% cfg.stim_save_filename = 'HBML/Desktop/stimuli001.wav';
% cfg.envelope_save_filename = 'HBML/Desktop/stimuli001.jpeg';
% cfg.plot_save_filename = 'HBML/Desktop/stimuli001.jpeg';
% [stimulus,envelope]=stim_creatorv2(cfg);

if ~isfield(cfg,'SNR')
    cfg.SNR = -3;
end
if ~isfield(cfg,'prespeech')
    cfg.prespeech.part1.noise = 'silence';
    cfg.prespeech.part1.length = 0;
end
if ~isfield(cfg,'postspeech')
    cfg.postspeech.part1.noise = 'silence';
    cfg.postspeech.part1.length = 0;
end
if ~isfield(cfg,'LvsR')
    cfg.LvsR = 'both';
end

% Set sample rate same as the speech signal
sound_info = audioinfo(cfg.speech.file);
SampleRate = sound_info.SampleRate;

%% Create speech part first
% Rational to start with speech is to use the same modulation factor across noise for each part

curr_part = cfg.speech;

% process speech signal
[speech_audio, ~] = audioread(curr_part.file);
if size(speech_audio,2)==2 % it may have two channels, get only one for now
    speech_audio(:,2) = [];
end

% remove zeros in the end, so that it will not affect power calculations
speech_audio = speech_audio(1:find(speech_audio,1,'last'));

curr_part.length = length(speech_audio)/SampleRate;


% generate noise
if strcmp(curr_part.noise,'silence')
    curr_noise = zeros(curr_part.length*SampleRate,1);
elseif any(strcmp(curr_part.noise,{'pink','white','brown','blue','purple'}))
    cn = dsp.ColoredNoise('Color',curr_part.noise,'SamplesPerFrame',SampleRate,'NumChannels',1);
    curr_noise = cn();
elseif exist(curr_part.noise,'file')
    [curr_noise, noise_SampleRate] = audioread(curr_part.noise);
    curr_noise = resample(curr_noise,SampleRate,noise_SampleRate);
else
    error('Unrecognized noise in speech!')
end
if size(curr_noise,2)==2 % it may have two channels, get only one for now
    curr_noise(:,2) = [];
end

% chop/extend the noise to get desired length (same with )
if length(curr_noise) > length(speech_audio)
    curr_noise = curr_noise(1:length(speech_audio),:);
elseif length(curr_noise) <  length(speech_audio)
    long_curr_noise = [curr_noise;curr_noise;curr_noise;curr_noise;curr_noise;curr_noise;curr_noise];
    curr_noise = long_curr_noise(1:length(speech_audio),:);
    clear long_curr_noise
end


% Combine signal and noise according to given SNR (based on whole
% segment power

Npts = length(speech_audio); % Number of input time samples

Signal_Power = sum(abs(speech_audio).*abs(speech_audio))/Npts;
Noise_Power = sum(abs(curr_noise).*abs(curr_noise))/Npts;
% Initial_SNR = 10*(log10(Signal_Power/Noise_Power))

K = (Signal_Power/Noise_Power)*10^(-cfg.SNR/10);  % Scale factor (this is going to be same across parts)

modulated_curr_noise = sqrt(K)*curr_noise; % Change Noise level
% New_Noise_Power = sum(abs(modulated_curr_noise).*abs(modulated_curr_noise))/Npts
% New_SNR = 10*(log10(Signal_Power/New_Noise_Power))

% Combine speech and noise
curr_stim = speech_audio + modulated_curr_noise;


% LvsR option
if strcmp(cfg.LvsR,'L') 
    speech_stim = [curr_stim,modulated_curr_noise];
elseif strcmp(cfg.LvsR,'R') 
    speech_stim = [modulated_curr_noise,curr_stim];
end

%% Second to-do: prespeech part
pre_fields = fieldnames(cfg.prespeech);
pre_fields=sortrows(pre_fields);    % sort the fields in the order.
prespeech_stim = [];
for pre_i = 1:length(fieldnames(cfg.prespeech))
    curr_part = cfg.prespeech.(pre_fields{pre_i});
    
    % generate signal (if present); if not, need a specified length to
    % produce silent period
    if isfield(curr_part,'signal')
        [curr_signal, sign_SampRate] = audioread(curr_part.signal);
        curr_signal=resample(curr_signal,SampleRate,sign_SampRate); % resample so that it will match with signal sample rate
        if size(curr_signal,2)==2 % it may have one channel
            curr_signal(:,2) = [];
        end
        
        % if wanting to chop/extend the signal to get desired length
        if isfield(curr_part,'length')
            if length(curr_signal) > curr_part.length*SampleRate
                curr_signal = curr_signal(1:curr_part.length*SampleRate,:);
            elseif length(curr_signal) < curr_part.length*SampleRate
                long_curr_signal = [curr_signal;curr_signal;curr_signal;curr_signal;curr_signal;curr_signal;curr_signal];
                curr_signal = long_curr_signal(1:curr_part.length*SampleRate,:);
                clear long_curr_signal
            end
        else
            curr_part.length = length(curr_signal)/SampleRate;
        end
    else
        curr_signal=zeros(curr_part.length*SampleRate,1);
    end
    
    % generate noise
    if strcmp(curr_part.noise,'silence')
        curr_noise = zeros(curr_part.length*SampleRate,1);
    elseif any(strcmp(curr_part.noise,{'pink','white','brown','blue','purple'}))
        cn = dsp.ColoredNoise('Color',curr_part.noise,'SamplesPerFrame',SampleRate,'NumChannels',1);
        curr_noise = cn();
    elseif exist(curr_part.noise,'file')
        [curr_noise, ~] = audioread(curr_part.noise);
        if size(curr_noise,2)==2 % it may have two channel, make it 1
            curr_noise(:,2) = [];
        end
    else
        error('Unrecognized noise in prespeech!')
    end
    
    % chop/extend the noise to get desired length
    if length(curr_noise) > curr_part.length*SampleRate
        curr_noise = curr_noise(1:curr_part.length*SampleRate,:);
    elseif length(curr_noise) < curr_part.length*SampleRate
        long_curr_noise = [curr_noise;curr_noise;curr_noise;curr_noise;curr_noise;curr_noise;curr_noise];
        curr_noise = long_curr_noise(1:curr_part.length*SampleRate,:);
        clear long_curr_noise
    end
    
    
    % Modulate noise with specific K-modulation factor
    curr_noise = sqrt(K)*curr_noise; % Change Noise level
    %     Npts = length(curr_noise); % Number of input time samples
    %
    % Signal_Power = sum(abs(curr_signal).*abs(curr_signal))/Npts;
    % Noise_Power = sum(abs(curr_noise).*abs(curr_noise))/Npts;
    % Initial_SNR = 10*(log10(Signal_Power/Noise_Power))
    %New_Noise_Power = sum(abs(curr_noise).*abs(curr_noise))/Npts
    %New_SNR = 10*(log10(Signal_Power/New_Noise_Power))
    
    curr_stim = curr_signal + curr_noise;
    
    % Collect each step
    % LvsR option
    if strcmp(cfg.LvsR,'L')
        prespeech_stim = [prespeech_stim;[curr_stim,curr_noise]];
    elseif strcmp(cfg.LvsR,'R')
        prespeech_stim = [prespeech_stim;[curr_noise,curr_stim]];
    end
    
end


%% Third: Post speech part generation
post_fields = fieldnames(cfg.postspeech);
post_fields=sortrows(post_fields);    % sort the fields in the order.
postspeech_stim = [];
for post_i = 1:length(fieldnames(cfg.postspeech))
    curr_part = cfg.postspeech.(post_fields{post_i});
    
    % generate noise
    if strcmp(curr_part.noise,'silence')
        curr_noise = zeros(curr_part.length*SampleRate,2);
    elseif any(strcmp(curr_part.noise,{'pink','white','brown','blue','purple'}))
        cn = dsp.ColoredNoise('Color',curr_part.noise,'SamplesPerFrame',SampleRate,'NumChannels',1);
        curr_noise = cn();
    elseif exist(curr_part.noise,'file')
        [curr_noise, ~] = audioread(curr_part.noise);
        if size(curr_noise,2)==2 % it may have two channel, make it 1
            curr_noise(:,2) = [];
        end
    else
        error('Unrecognized noise in prespeech!')
    end
    
    % chop/extend the noise to get desired length
    if length(curr_noise) > curr_part.length*SampleRate
        curr_noise = curr_noise(1:curr_part.length*SampleRate,:);
    elseif length(curr_noise) < curr_part.length*SampleRate
        long_curr_noise = [curr_noise;curr_noise;curr_noise;curr_noise;curr_noise;curr_noise;curr_noise];
        curr_noise = long_curr_noise(1:curr_part.length*SampleRate,:);
        clear long_curr_noise
    end
    
    % generate signal (if present)
    if isfield(curr_part,'signal')
        [curr_signal, sign_SampRate] = audioread(curr_part.signal);
        curr_signal=resample(curr_signal,SampleRate,sign_SampRate); % resample so that it will match with signal sample rate
        if size(curr_signal,2)==2 % it may have one channel
            curr_signal(:,2) = [];
        end
        % chop/extend the signal to get desired length
        if length(curr_signal) > curr_part.length*SampleRate
            curr_signal = curr_signal(1:curr_part.length*SampleRate,:);
        elseif length(curr_signal) < curr_part.length*SampleRate
            long_curr_signal = [curr_signal;curr_signal;curr_signal;curr_signal;curr_signal;curr_signal;curr_signal];
            curr_signal = long_curr_signal(1:curr_part.length*SampleRate,:);
            clear long_curr_signal
        end
    else
        curr_signal=zeros(curr_part.length*SampleRate,1);
    end
    
    
%     % Modulate noise with specific K-modulation factor
    curr_noise = sqrt(K)*curr_noise; % Change Noise level
    
    curr_stim = curr_signal + curr_noise;
    curr_stim=[curr_stim,curr_stim]; % make it same for each ear
    
    % Collect each step
    postspeech_stim=[postspeech_stim;curr_stim];
    
end

%% Combine all parts of stimuli
stimulus = [prespeech_stim;speech_stim;postspeech_stim];

% Save (optional)
if isfield(cfg,'stim_save_filename')
    stimulus(stimulus<-1)=-1; % in order to be able to save without clipping
    stimulus(stimulus>1)=1; 
    fprintf('\tStimulus is saved to: %s\n',cfg.stim_save_filename)
    information = ['SNR is ' num2str(cfg.SNR) 'dB; '...
        'Prespeech has ' num2str(length(pre_fields)) 'part(s); '...
        'Postspeech has ' num2str(length(post_fields)) 'part(s); '...
        'Audio is given to ' cfg.LvsR ' ear(s). Other is for stimulator'];
    audiowrite(cfg.stim_save_filename,stimulus,SampleRate,'Comment',information);
end

% Add delay and create second channel that goes to stimulator
if isfield(cfg,'delay')
    if isnumeric(cfg.delay)
        stim_channel = [zeros(1,(length(prespeech_stim)+(cfg.delay*SampleRate))),speech_audio',zeros(1,(length(postspeech_stim)-(cfg.delay*SampleRate)))];
    elseif strcmpi(cfg.delay, 'attention')
        attention_onset = cfg.prespeech.(pre_fields{1}).length;
        
        if ~isfield(cfg,'frequency') || ~isfield(cfg,'phase'),error('Need to have frequency and phase inside, in order to generate sine wave'),end
        % Create sine wave
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Need to create here %%%%%%%%%%%%%%%%
        
        
        % stim_channel = zeros(1,length(attention_onset*SampleRate)), 
                %           wave + zeros
    else
        error('Unrecognized delay input')
    end
else
    stim_channel = [zeros(1,length(prespeech_stim)),speech_audio',zeros(1,length(postspeech_stim))];
end

% LvsR option: one side is stimulus, the other is only clean speech (will go into stimulator)
if strcmp(cfg.LvsR,'L')
    stimulus = [stimulus(:,1),stim_channel'];
    presented_stim=stimulus(:,1);
elseif strcmp(cfg.LvsR,'R')
    stimulus = [stim_channel',stimulus(:,2)];
    presented_stim=stimulus(:,2);
end


%% Create envelope of speech (Based on 2018_Kosem scripts)

%compute envelope of speech only
% 1. band pass filtering between 100 Hz and 4000 Hz
filter_range = [100 4000];
[b,a] = butter(2, filter_range/(SampleRate/2), 'bandpass');
yfilt = filter(b,a,speech_audio);

% 2. Compute Hilbert transform and take the absolute value to get the
% amplitude of the envelope
amp = abs(hilbert(yfilt));

% embed envelope inside whole stim length 
envelope = [zeros(1,length(prespeech_stim)),amp',zeros(1,length(postspeech_stim))];

% compute envelope's power spectrum

L = length(amp);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(amp,NFFT)/L;
freqf = SampleRate/2*linspace(0,1,NFFT/2+1);
pow = abs(Y(1:NFFT/2+1)).^2;

% Optional save for envelope
if isfield(cfg,'envelope_save_filename')
    save(cfg.envelope_save_filename,'envelope')
end

%% Make plot (optional)

if isfield(cfg,'plot_save_filename')
h1=figure('Position', [50 50  1000 1200]);
% Plot findpeaks in first row
subplot(4,1,1)
[pks,locs] = findpeaks(pow(freqf>0.7&freqf<20)/max(pow((freqf>0.7&freqf<20))),freqf(freqf>0.7&freqf<20));
findpeaks(pow/max(pow((freqf>0.7))),freqf)
text(locs+.02,pks,num2str(locs',3));
xlim([1 8])
ylim([0,1])
title(['Peaks in current sentence; speech is filtered in ' num2str(filter_range) 'Hz'])
set(gca, 'fontname','arial'); 

% plot envelope in 2nd row
subplot(4,1,2)
time_axis = 1/SampleRate:1/SampleRate:length(envelope)/1/SampleRate;
plot(time_axis,envelope)
xlim([0 length(envelope)/1/SampleRate])
xlabel('Time (s)')
ylabel('Amplitude')
title('Envelope of the speech in whole stimulus')

% plot what patient hears in 3rd row
subplot(4,1,3)
time_axis = 1/SampleRate:1/SampleRate:length(envelope)/1/SampleRate;
plot(time_axis,presented_stim)
xlim([0 length(envelope)/1/SampleRate])
xlabel('Time (s)')
ylabel('Amplitude')
title('Presented stimulus: what goes in patients ear.')

% plot what goes into stimulator in 4th row
subplot(4,1,4)
plot(time_axis,stim_channel)
xlim([0 length(envelope)/1/SampleRate])
xlabel('Time (s)')
ylabel('Amplitude')
title('Stimulator channel: what goes into stimulator.')


[~,name,~] = fileparts(cfg.speech.file);
name = replace(name,'_',' ');
h1.NextPlot = 'add';
a = axes;
%// Set the title and get the handle to it
ht = title({name;''});
%// Turn the visibility of the axes off
a.Visible = 'off';
%// Turn the visibility of the title on
ht.Visible = 'on';

fprintf('\tPlot is saved to :%s\n',cfg.plot_save_filename)
print('-r300','-djpeg',cfg.plot_save_filename)
close all

end


end