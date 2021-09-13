function [amp_env,deriv_amp_env,spectrogram,CF] = get_speech_2features(speech,Fs,scale)
% intermediary function for get_speech_peaks.m and also can be used to
% directly get the amplitude envelope and the first derivative of the
% amplitude envelope (no rectification).
% Serdar Akkol, June 2020, HBML.

if ~exist('scale', 'var') || isempty(scale) || strcmpi(scale,'cochlear')
    CF = 440 * 2 .^ ((-31:97)/24);
elseif strcmpi(scale,'bark')
    load('bark_cutoff_freqs.mat','bark_cutoff_freqs') % I tried Bark scale but results weren't clear
    CF = bark_cutoff_freqs;
end

% Remove frequency bands that exceeds half of sampling rate
CF(CF>Fs/2) = [];

% Calculate narrow-band filtered speech envelope
for bi = 1:[length(CF)-1]
    filter_range = [CF(bi) CF(bi+1)];
    [b,a] = butter(2, filter_range/(Fs/2), 'bandpass');
    yfilt = filter(b,a,speech);
    spectrogram(bi,:) = abs(hilbert(yfilt));
end

% Average narrow band filtered signals to get broadband envelope (with minor smoothing)
amp_env = smooth(mean(spectrogram,1),0.01*Fs);

% Calculate the first derivative of envelope (which means changes in rate
% of envelope)  (with minor smoothing)
deriv_amp_env = [smooth(diff(amp_env),0.01*Fs);0];

end
