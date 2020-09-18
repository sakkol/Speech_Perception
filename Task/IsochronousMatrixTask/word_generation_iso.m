%% Load table
sentence_table = readtable('/home/sakkol/Documents/TASKS/PREPARATIONS/IsochronousMatrixTask/English/vol1/English_single_words.xlsx');

%% create many many sentences and write them to individual files
save_dir = '/home/sakkol/Documents/TASKS/PREPARATIONS/IsochronousMatrixTask/English/vol1/word-texts';

for y=1:6
    for x=1:5
        curr_sent = sentence_table{y,x}{1};
        fid = fopen([save_dir filesep curr_sent '.txt'],'w+');
        fwrite(fid,curr_sent);
        fclose(fid);
    end
end


%% After generating sounds through GoogleCloud_TTS_iso.py, change the sound volume

sound_dir = '/home/sakkol/Documents/TASKS/PREPARATIONS/IsochronousMatrixTask/English/vol1/word_sounds';
save_dir = fullfile(sound_dir,'word_spectrograms');
soundss = dir([sound_dir '/*.wav']);

WordsPeakInfo=table;
for s = 1:length(soundss)
    filename = fullfile(soundss(s).folder,soundss(s).name);
    [curr_sound,Fs] = audioread(filename);
    
    % normalize to its max
    new_sound = curr_sound/max(curr_sound);
    new_sound = new_sound(find(new_sound,1,'first'):find(new_sound,1,'last'));
    audiowrite(fullfile(soundss(s).folder,soundss(s).name),new_sound,Fs);
    
    % save the spectrogram plot
    WordsPeakInfo.StimName{s} = soundss(s).name;
    
    [peakRate, peakEnv, amp_envel, deriv_amp_env] = get_speech_peaks(resample(new_sound,24000,40000),24000,1);
    
    WordsPeakInfo.peak_info{s} = table;
    WordsPeakInfo.peak_info{s}.peakEnv{1} = peakEnv/40000;
    WordsPeakInfo.peak_info{s}.peakRate{1} = peakRate/40000;
    
    % add title
    sgtitle(['Word is: ' erase(soundss(s).name,'-M.wav')], 'FontSize',16,'FontWeight','bold')
    print(fullfile(save_dir,[erase(soundss(s).name,'-M.wav') '.jpg']),'-djpeg','-r300')
    close all
    
end

save(fullfile(save_dir,'WordsPeakInfo.mat'),'WordsPeakInfo')