%% Load table
sentence_table = readtable('/home/sakkol/Documents/TASKS/PREPARATIONS/IsochronousMatrixTask/vol1/English_single_words.xlsx');

% create many many sentences and write them to individual files
save_dir = '/home/sakkol/Documents/TASKS/PREPARATIONS/IsochronousMatrixTask/vol1/word-texts/English';

for y=1:6
    for x=1:5
        curr_sent = sentence_table{y,x}{1};
        fid = fopen([save_dir filesep curr_sent '.txt'],'w+');
        fwrite(fid,curr_sent);
        fclose(fid);
    end
end

%% Several words need extra care
[soundtocare,SampleRate] = audioread('/home/sakkol/Documents/TASKS/PREPARATIONS/IsochronousMatrixTask/vol1/word_sounds/English/Jane, got, Sue, to, buy, Sean, two, cheap, cups-M.wav');

sound_dir = '/home/sakkol/Documents/TASKS/PREPARATIONS/IsochronousMatrixTask/vol1/word_sounds/English';

audiowrite(fullfile(sound_dir,'Jane-M.wav'),soundtocare(923:12921),SampleRate);
audiowrite(fullfile(sound_dir,'Sue-M.wav'),soundtocare(30700:42475),SampleRate);
audiowrite(fullfile(sound_dir,'Sean-M.wav'),soundtocare(73031:86086),SampleRate);
audiowrite(fullfile(sound_dir,'cheap-M.wav'),soundtocare(103587:113300),SampleRate);
audiowrite(fullfile(sound_dir,'cups-M.wav'),soundtocare(115500:128074),SampleRate);

% for now and kept
[soundtocare,SampleRate] = audioread('/home/sakkol/Documents/TASKS/PREPARATIONS/IsochronousMatrixTask/vol1/word_sounds/English/I, kept, now-M.wav');

sound_dir = '/home/sakkol/Documents/TASKS/PREPARATIONS/IsochronousMatrixTask/vol1/word_sounds/English';

audiowrite(fullfile(sound_dir,'kept-M.wav'),soundtocare(8851:19190),SampleRate);
audiowrite(fullfile(sound_dir,'now-M.wav'),soundtocare(23323:36000),SampleRate);

%% After generating sounds through GoogleCloud_TTS_iso.py, change the sound volume

sound_dir = '/home/sakkol/Documents/TASKS/PREPARATIONS/IsochronousMatrixTask/vol1/word_sounds/English';
save_dir = fullfile(sound_dir,'word_spectrograms');
soundss = dir([sound_dir '/*.wav']);

WordsInfo=table;
for s = 1:length(soundss)
    filename = fullfile(soundss(s).folder,soundss(s).name);
    [curr_sound,Fs] = audioread(filename);
    
    % normalize to its max
    new_sound = curr_sound/max(abs(curr_sound));
    verylows = new_sound<0.005 & new_sound>-0.005;
%     new_sound = new_sound(find(new_sound,1,'first'):find(new_sound,1,'last'));
    new_sound = new_sound(find(~verylows,1,'first'):find(~verylows,1,'last'));
    audiowrite(fullfile(soundss(s).folder,soundss(s).name),new_sound,Fs);
    
    % save the spectrogram plot
    WordsInfo.StimName{s} = erase(soundss(s).name,'-M.wav');
    WordsInfo.Stim{s} = new_sound;
    WordsInfo.StimLength(s) = length(new_sound)/Fs;
    WordsInfo.SampleRate(s) = Fs;
    
    [peakRate, peakEnv, amp_envel, deriv_amp_env] = get_speech_peaks(resample(new_sound,24000,Fs),24000,1);
    WordsInfo.StimLength(s) = length(new_sound)/Fs;
    WordsInfo.peak_info{s} = table;
    WordsInfo.peak_info{s}.peakEnv{1} = peakEnv/Fs;
    WordsInfo.peak_info{s}.peakRate{1} = peakRate/Fs;
    
    % add title
    sgtitle(['Word is: ' erase(soundss(s).name,'-M.wav')], 'FontSize',16,'FontWeight','bold')
    print(fullfile(save_dir,[erase(soundss(s).name,'-M.wav') '.jpg']),'-djpeg','-r300')
    
%     obj = audioplayer(curr_sound,Fs);
%     playblocking(obj);
    close all
    
end


% Estimating the K value, to multiply with noise for SNR
Signal_Power=[];Noise_Power=[];K=[];
for s=[22,9,31,35]
    % generate noise
    cn = dsp.ColoredNoise('Color','pink','SamplesPerFrame',44100,'NumChannels',1);
    curr_noise = cn();clear cn;curr_noise=[curr_noise;curr_noise;curr_noise;curr_noise;curr_noise];
    if size(curr_noise,2)==2 % it may have two channels, get only one for now
        curr_noise(:,2) = [];
    end
    
    % Combine signal and noise according to given SNR (based on whole
    % segment power
    speech_audio=WordsInfo.Stim{s};
    Npts = length(speech_audio); % Number of input time samples
    
    Signal_Power(s,1) = sum(abs(speech_audio).*abs(speech_audio))/Npts;
    Noise_Power(s,1) = sum(abs(curr_noise(1:Npts)).*abs(curr_noise(1:Npts)))/Npts;
    
%     K(s,1) = (Signal_Power(s,1)/Noise_Power(s,1))*10^(-cfg.SNR/10);  % Scale factor (this is going to be same across parts)
%     WordsPeakInfo.sqrtScaleFactor = sqrt(K);
end

Signal_Power=Signal_Power(Signal_Power~=0);
Noise_Power=Noise_Power(Noise_Power~=0);

avg_sig_pow = mean(Signal_Power);
avg_noise_pow = mean(Noise_Power);

words_table = readtable('/home/sakkol/Documents/TASKS/PREPARATIONS/IsochronousMatrixTask/vol1/English_single_words.xlsx');

% save several things here
save(fullfile(save_dir,'EnglishWordsInfo.mat'),'WordsInfo','avg_sig_pow','avg_noise_pow','words_table')
