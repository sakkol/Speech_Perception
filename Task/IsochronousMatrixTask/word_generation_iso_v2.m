%% Load table
sentence_table = readtable('/home/sakkol/Documents/TASKS/PREPARATIONS/IsochronousMatrixTask/vol2/English_single_words.xlsx');

% create many many sentences and write them to individual files
save_dir = '/home/sakkol/Documents/TASKS/PREPARATIONS/IsochronousMatrixTask/vol2/word-texts/English';

% create text files for python
for y=1:7
    curr_sent = strjoin(table2cell(sentence_table(y,:)),', ');
    fid = fopen([save_dir filesep curr_sent '.txt'],'w+');
    fwrite(fid,curr_sent);
    fclose(fid);
end

% once they are created
sound_dir = '/home/sakkol/Documents/TASKS/PREPARATIONS/IsochronousMatrixTask/vol2/word_sounds/English';

for y=1:7
    curr_sent = strjoin(table2cell(sentence_table(y,:)),', ');
    
    [soundtocare,SampleRate] = audioread(fullfile(sound_dir,[curr_sent '-M.wav']));
    
    thr_ampl = 0.0005;
    refract_tpts = floor(0.33*SampleRate); % 1 second only
    digital_trig_chan=analog2digital_trig(soundtocare',thr_ampl,refract_tpts,0);
    trial_onsets_tpts=find(digital_trig_chan==1);trial_onsets_tpts=trial_onsets_tpts-2000;trial_onsets_tpts(trial_onsets_tpts<0)=1;
    
    digital_trig_chan=analog2digital_trig(soundtocare(end:-1:1)',thr_ampl,refract_tpts,0);
    trial_ends_tpts=find(digital_trig_chan==1);trial_ends_tpts = length(soundtocare) - trial_ends_tpts;
    trial_ends_tpts = trial_ends_tpts + 2000;trial_ends_tpts = trial_ends_tpts(end:-1:1);
    
    figure('Units','normalized','Position', [0 0  1 .5]);
    plot(1/SampleRate:1/SampleRate:(length(soundtocare)/SampleRate),soundtocare); hold on
    ylims=ylim;
    for w=1:5%length(trial_onsets_tpts)
        plot([trial_onsets_tpts(w)/SampleRate trial_onsets_tpts(w)/SampleRate],ylim)
        plot([trial_ends_tpts(w)/SampleRate trial_ends_tpts(w)/SampleRate],ylim)
        text(trial_onsets_tpts(w)/SampleRate,4*ylims(2)/5,sentence_table{y,w}{1})
    end
    ylim(ylims)
    title([num2str(length(trial_onsets_tpts)) ' onsets found: ' curr_sent]);
    xlabel('Time (s)')
    xlim([1/SampleRate length(soundtocare)/SampleRate])
    print(fullfile(sound_dir,[curr_sent '-M.jpg']),'-djpeg')
    
%     pause(2)
    
    for w = 1:5
        curr_word = sentence_table{y,w}{1};
        for_crr = soundtocare(trial_onsets_tpts(w):trial_ends_tpts(w));
        for_crr = for_crr(find(for_crr,1,'first'):find(for_crr,1,'last'));
        verylows = for_crr<0.0014 & for_crr>-0.0014;
        for_crr = for_crr(find(~verylows,1,'first'):find(~verylows,1,'last'));
        for_crr = for_crr/max(abs(for_crr));
        audiowrite(fullfile(sound_dir,[curr_word '-M.wav']),for_crr,SampleRate);
    end
end





%% After generating sounds through GoogleCloud_TTS_iso.py, change the sound volume

sound_dir = '/home/sakkol/Documents/TASKS/PREPARATIONS/IsochronousMatrixTask/vol2/word_sounds/English';
save_dir = fullfile(sound_dir,'word_spectrograms');mkdir(save_dir)
soundss = dir([sound_dir '/*.wav']);

WordsInfo=table;
for s = 1:length(soundss)
    filename = fullfile(soundss(s).folder,soundss(s).name);
    [new_sound,Fs] = audioread(filename);
    
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
for s=[25,10,39,35]  % now, catch, these, words
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

words_table = readtable('/home/sakkol/Documents/TASKS/PREPARATIONS/IsochronousMatrixTask/vol2/English_single_words.xlsx');

% save several things here
save(fullfile(save_dir,'EnglishWordsInfo.mat'),'WordsInfo','avg_sig_pow','avg_noise_pow','words_table')

%% Words in adaptation and threshold
def=load('default_conditions.mat');

% prepare adaptation and threshold word list
[events_table] = events_wrapper([],'English','both',[]) % choose with selections_GUI
cfgs_for_adapt_English = events_table.cfgs;

adaptation_word_list_English = cell(10,20);
for t = 1:10
    for p=1:5
        for w=1:4
            adaptation_word_list_English{t,(4*(p-1))+w} = cfgs_for_adapt_English{t}.(['part' num2str(p+1)]).(['word' num2str(w)]){:};
        end
    end
end

[events_table] = events_wrapper([],'English','both',[]) % choose with selections_GUI
cfgs_for_thresh_English = events_table.cfgs;

for t = 1:30
    if mod(t,3)==0
        cfgs_for_thresh_English{t}.SNR = -6;
    elseif mod(t,3)==1
        cfgs_for_thresh_English{t}.SNR = -2;
    elseif mod(t,3)==2
        cfgs_for_thresh_English{t}.SNR = 2;
    end
end
cfgs_for_thresh_English = cfgs_for_thresh_English(randperm(size(cfgs_for_thresh_English,1)));
threshold_word_list_English = cell(30,5);
for t = 1:30
    threshold_word_list_English{t,1} = cfgs_for_thresh_English{t}.SNR;
    for p=1
        for w=1:4
            
            threshold_word_list_English{t,(4*(p-1))+w+1} = cfgs_for_thresh_English{t}.(['part' num2str(p+2)]).(['word' num2str(w)]){:};
        end
    end
end

% now get the defaults again
iso_mat_24 = def.iso_mat_24;
FreeRecall_3sent = def.FreeRecall_3sent; % changing to 1.5Hz
def.FreeRecall_3sent.Frequency(1) = 1.5;
def.FreeRecall_3sent.Frequency(2) = 1.5;
adapt_thresh = def.adapt_thresh;
clear_table = def.clear_table;
passive_5sent = def.passive_5sent;

% Create all stimuli to select from for Free Recall
events_table1 = events_wrapper([],'English','both','FreeRecall_3sent'); % only sentence, 105 trials
events_table2 = events_wrapper([],'English','both','FreeRecall_3sent'); % only scrambled, 105 trials
words_list = cell(210,12);
for t=1:105
    words_list(t,:) = [events_table1.trial_words{t}];
end
for t=1:105
    words_list(t+105,:) = [events_table2.trial_words{t}];
end
events_table = [events_table1;events_table2];
cfgs_for_FreeRecall_English = table2cell(events_table(:,2:4));

save('default_conditions.mat','iso_mat_24','passive_5sent','FreeRecall_3sent','clear_table','adapt_thresh',...
    'cfgs_for_adapt_English','cfgs_for_thresh_English','cfgs_for_FreeRecall_English','words_list','adaptation_word_list_English','threshold_word_list_English')


%% Create example stimuli
sound_dir = '/home/sakkol/Documents/TASKS/PREPARATIONS/IsochronousMatrixTask/vol2/word_sounds/English';
save_dir = fullfile(sound_dir,'example_stim');mkdir(save_dir)

events_table = events_wrapper([],'English','both','FreeRecall_3sent');

for t=1:30
    audiowrite(fullfile(save_dir,[num2str(t) '-' events_table.cond_info{t}.Number_of_sentences{1} '-' events_table.cond_info{t}.Sentence_vs_Scrambled{1} '-.wav']),...
        events_table.trials{t},44100)
end
    

