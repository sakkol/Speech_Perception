function [peakRate, peakEnv, amp_envel, deriv_amp_env] = get_speech_peaks(speech,Fs)


% CF = 440 * 2 .^ ((-31:97)/24);
load('bark_cutoff_freqs.mat')
CF = bark_cutoff_freqs;




for bi = 1:[length(CF)-1]
    filter_range = [CF(bi) CF(bi+1)];
    [b,a] = butter(2, filter_range/(Fs/2), 'bandpass');
    yfilt = filter(b,a,speech);
    hilby(bi,:) = abs(hilbert(yfilt));
end



amp_envel = smooth(mean(hilby,1),1500);


deriv_amp_env = [smooth(diff(amp_envel),1500);0];



% amp_envel = smooth(mean(hilby,1),500);
% amp_envel = sum(hilby,1);
figure
subplot(311)
imagesc(hilby)
axis xy
subplot(312)
[~,peakEnv]=findpeaks(amp_envel,'MinPeakDistance',1000);

subplot(313)
[~,peakRate]=findpeaks(deriv_amp_env,'MinPeakDistance',1000);

% remove last peaks that are close to end (probably error)
peakEnv( peakEnv > length(deriv_amp_env)-1000) = [];
peakRate( peakRate > length(deriv_amp_env)-1000) = [];





