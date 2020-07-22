%% Adding peak info to English all info
main_stim_loc = ...
    '/home/sakkol/Documents/BACKUPS/Speech_Perception_BACKUP(2019.09.11)/SpeechPerception_Task(backup-11.09.2019)/SpeechPerception_Task/main_stim_loc';
load('English_all_info_old.mat')

savedir = fullfile(main_stim_loc,'peakInfo_plots');
ifPlot = 0;

fprintf('Iteration = 000');
for i=1:size(all_info_table,1)
    fprintf('\b\b\b%3d',i)
    
%     all_info_table.all_info{i}.Properties.VariableNames{'onset'} = 'syllable_onset';
%     all_info_table.all_info{i}.Properties.VariableNames{'offset'} = 'syllable_offset';
    
    
    if strcmp(all_info_table.rate{i},'Fast')
        speech_rate = 1.3;
    elseif strcmp(all_info_table.rate{i},'Slow')
        speech_rate = 0.9;
    end
    
    filename = find_sentence(all_info_table.sentence{i},main_stim_loc,speech_rate);
    
    [speech,Fs] = audioread(filename);
    
    [peakRate, peakEnv, amp_envel, deriv_amp_env] = get_speech_peaks(speech,Fs,ifPlot);
    
    all_info_table.peak_info{i} = table;
    all_info_table.peak_info{i}.peakEnv{1} = peakEnv/Fs;
    all_info_table.peak_info{i}.peakRate{1} = peakRate/Fs;
    
    if ifPlot
        % add some more info, word boundaries
        subplot(211)
        hold on
        y=ylim;
        onsoff = [all_info_table.word_info{i}.onset;all_info_table.word_info{i}.offset(end)];
        z1=plot([onsoff onsoff],y,'w','LineWidth',2);
        splsent = strsplit(all_info_table.sentence{i},' ');
        for s = 1:length(splsent)
            text((onsoff(s)+onsoff(s+1))/2,y(2)-700,splsent{s},'HorizontalAlignment','center','FontSize',16,'Color','w','FontWeight','bold')
        end
        set(gca,'XGrid','off');
        
        subplot(212)
        set(gca,'XGrid','off');
        sgtitle([all_info_table.sentence{i} ' (white lines: word boundaries)'], 'FontSize',16,'FontWeight','bold')
        
        print(fullfile(savedir,['Rate_' num2str(speech_rate)],[strjoin(splsent,'_') '.jpg']),'-djpeg','-r300')
        close all
    end
end
fprintf('\nSaving...\n')

save(fullfile('/home/sakkol/Documents/Codes_git/Speech_Perception/preproc/English_all_info.mat'),'all_info_table')

%% Creating peak info for Spanish version
main_stim_loc = '/home/sakkol/Documents/TASKS/Matrix_Speech_Task/Spanish_main_stim_loc';

savedir = '/home/sakkol/Documents/TASKS/Matrix_Speech_Task_additionals/Spanish_peakInfo_plots';
ifPlot = 1;

speech_rates = [0.9,1.3];
sent2use = readtable(fullfile(main_stim_loc,'Sentence_to_use.xlsx'));
allsents = sent2use.Sentence;
peakInfotable=struct;
for s = 1:2
    speech_rate = speech_rates(s);
    mkdir(fullfile(savedir,['Rate_' num2str(speech_rate)]))
    
    fprintf('Iteration = 000');
    for i=1:size(allsents,1)
        fprintf('\b\b\b%3d',i)
        
        peakInfotable.lang{i} = 'Spanish';
        if speech_rate == 1.3
            peakInfotable.rate{i}='Fast';
        elseif speech_rate == 0.9
            peakInfotable.rate{i}='Slow';
        end
        
        
        filename = find_sentence(allsents{i},main_stim_loc,speech_rate);
        [speech,Fs] = audioread(filename);
        
        [peakRate, peakEnv, amp_envel, deriv_amp_env] = get_speech_peaks(speech,Fs,ifPlot);
        
        peakInfotable.peak_info{i} = table;
        peakInfotable.peak_info{i}.peakEnv{1} = peakEnv/Fs;
        peakInfotable.peak_info{i}.peakRate{1} = peakRate/Fs;
        
        if ifPlot
            % add some more info, word boundaries
            sgtitle(['Sentence is: ' allsents{i}], 'FontSize',16,'FontWeight','bold')
            print(fullfile(savedir,['Rate_' num2str(speech_rate)],[replace(allsents{i},' ','_') '.jpg']),'-djpeg','-r300')
            close all
        end
    end
end
fprintf('\nSaving...\n')
peakInfotable = struct2table(peakInfotable);
save(fullfile('/home/sakkol/Documents/Codes_git/Speech_Perception/preproc/Spanish_all_info.mat'),'peakInfotable')