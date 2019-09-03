function [stim,envelope]=stimuli_creator(cfg)
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
%                   noise, using the SNR above.
% cfg.speech.file : speech file to be used. (also signal)
%           .noise: same options as in prespeech (length and file).
% cfg.postspeech  : same options as in prespeech (length and file).
% cfg.LvsR        : 'L','R' or 'both' (defaulft: 'both'). If not both,
%                   speech signal is going to be given only from one ear.
% cfg.save.dir    
%         .name  
%                 : if there is this field, speech will be saved to dir
%                   with name.
%
% OUTPUT:
% stim: generated stimuli.
% envelope: envelope of signal (speech) in the stimuli. If there is noise
% before and/or after signal, those parts are going to stay zero.
% 
% Example:
% cfg=[];
% cfg.SNR = -4;
% cfg.prespeech.part1.length=0.5;
% cfg.prespeech.part1.noise = 'silence';
% cfg.prespeech.part2.length = 2;
% cfg.prespeech.part2.noise = 'pink';
% cfg.prespeech.part2.signal = 'babble.wav';
% cfg.speech.file = 'example_speech.wav';
% cfg.speech.noise = 'pink';
% cfg.LvsR = 'L';
% cfg.save.dir = 'HBML/Desktop';
% cfg.save.dir = 'stimuli001.wav';
% [stim,envelope]=stimuli_creator(cfg)

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

%% Start with prespeech part
pre_fields = fieldnames(cfg.prespeech);
pre_fields=sortrows(pre_fields);    % sort the fields in the order.
prespeech_stim = [];
for pre_i = 1:length(fieldnames(cfg.prespeech))
    curr_part = cfg.prespeech.(pre_fields{pre_i});
    
    % generate noise
    if strcmp(curr_part.noise,'silence')
        curr_noise = zeros(curr_part.length*SampleRate,2);
    elseif any(strcmp(curr_part.noise,{'pink','white','brown','blue','purple'}))
        cn = dsp.ColoredNoise('Color',curr_part.noise,'SamplesPerFrame',SampleRate,'NumChannels',2);
        curr_noise = cn();
    elseif exist(curr_part.noise,'file')
        [curr_noise, ~] = audioread(curr_part.noise);
        if size(curr_noise,2)==1 % it may have one channel
            curr_noise(:,2) = curr_noise(:,1);
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
        if size(curr_signal,2)==1 % it may have one channel
            curr_signal(:,2) = curr_signal(:,1);
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
        curr_signal=zeros(curr_part.length*SampleRate,2);
    end
    
    % Combine signal and noise according to given SNR (based on whole
    % segment power (take mean from ears)
    
    Npts = length(curr_signal); % Number of input time samples
    Noise = randn(1,Npts); % Generate initial noise; mean zero, variance one
    
    Signal_Power = sum(abs(curr_signal).*abs(curr_signal))/Npts;
    Noise_Power = sum(abs(curr_noise).*abs(curr_noise))/Npts;
    %Initial_SNR = 10*(log10(Signal_Power/Noise_Power));
    
    K = (Signal_Power/Noise_Power)*10^(-cfg.SNR/10);  % Scale factor
    
    curr_noise = sqrt(K)*curr_noise; % Change Noise level
    %New_Noise_Power = sum(abs(New_Noise).*abs(New_Noise))/Npts
    %New_SNR = 10*(log10(Signal_Power/New_Noise_Power))
    
    curr_stim = curr_signal + curr_noise;
    
    % Collect each step
    prespeech_stim=[prespeech_stim;curr_stim];
    
end





end