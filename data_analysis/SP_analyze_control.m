%% TO-DO:
% Add Rayleigh stats, because it is important only for control
control_wlt = addRayleigh_to_ft(control_wlt); % ?????????

%% Prepare
data_root = '/media/sakkol/HDD1/HBML/';
project_name = 'Speech_Perception';
sbj_ID = 'NS148';
Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'

% Select blocks to import
blocklistsheet = [char(Sbj_Metadata.project_root) filesep char(Sbj_Metadata.project_name) '_BlockInfo.xlsx']; % "F:\HBML\PROJECTS_DATA\CL_Train\CL_Train_BlockLists.xlsx";
blocklistall = readtable(blocklistsheet);
control_blocks = blocklistall.BlockList(strcmpi(blocklistall.sbj_ID,Sbj_Metadata.sbj_ID) & contains(blocklistall.conditions_code,'1'));

clear blocklistall blocklistsheet

%% bring in these blocks and combine only the control events

for b = 1:length(control_blocks)
    curr_block = control_blocks{b};
    
    % Load iEEG
    load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_wlt.mat']))
    events = epoched_wlt.events;
    
    % Select only control events
    control_idx = events.Cond_code == 1;
    
    % Select ecog data
    cfg = [];
    cfg.trials = control_idx;
    [epoched_wlt.wlt] = ft_selectdata(cfg, epoched_wlt.wlt);
    
    % Select events
    events = events(control_idx,:);
    
    % Append to overall list
    if b == 1
        control_wlt = epoched_wlt.wlt;
        control_events = events;
    else
        cfg = [];
        cfg.parameter  = 'fourierspctrm';
        control_wlt = ft_appendfreq(cfg, control_wlt, epoched_wlt.wlt);
        
        control_events = [control_events;events];
    end
    
end

%% Separate fourier spectrums of correct responses and others
corr_rspn_fft_word = {0};
no_rspn_fft_word = {0};
wrng_rspn_fft_word = {0};

corr_rspn_fft_syll = {0};
no_rspn_fft_syll = {0};
wrng_rspn_fft_syll = {0};

corr_rspn_fft_pho = {0};
no_rspn_fft_pho = {0};
wrng_rspn_fft_pho = {0};

% Which freq bands
freq_bands = {[1 4], [4 8],[8 12]};


for f = 1:length(freq_bands)
    cl_freqs = [nearest(control_wlt.freq, freq_bands{f}(1)), ...
        nearest(control_wlt.freq,freq_bands{f}(2))];
    
    wcorr_ind = 1;
    wno_ind = 1;
    wwrng_ind = 1;
    scorr_ind = 1;
    sno_ind = 1;
    swrng_ind = 1;
    pcorr_ind = 1;
    pno_ind = 1;
    pwrng_ind = 1;
    
    for t = 1:size(control_events,1)
        
        for w = 1:5
            % Separate words
            % Check if correct
            % Collect fft in a cell structure
            if strcmp(control_events.word_info{t}.response{w},'1')
                % find closest timepoint
                cl_times = [nearest(control_wlt.time, control_events.word_info{t}.onset(w))...
                    nearest(control_wlt.time,control_events.word_info{t}.offset(w))];
                corr_rspn_fft_word{wcorr_ind,f} = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)));
                wcorr_ind = wcorr_ind+1;
            elseif strcmp(control_events.word_info{t}.response{w},'0')
                % find closest timepoint
                cl_times = [nearest(control_wlt.time, control_events.word_info{t}.onset(w))...
                    nearest(control_wlt.time,control_events.word_info{t}.offset(w))];
                no_rspn_fft_word{wno_ind,f} = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)));
                wno_ind = wno_ind+1;
                
            else % wrong responses
                % find closest timepoint
                cl_times = [nearest(control_wlt.time, control_events.word_info{t}.onset(w))...
                    nearest(control_wlt.time,control_events.word_info{t}.offset(w))];
                wrng_rspn_fft_word{wwrng_ind,f} = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)));
                wwrng_ind = wwrng_ind+1;
                
            end
            
        end
        
        for s = 1:size(control_events.syllable_info{t},1)
            % Separate syllables like different trials
            % Check if correct
            % Collect fft in a cell structure
            if strcmp(control_events.word_info{t}.response{control_events.syllable_info{t}.word_id(s)},'1')
                % find closest timepoint
                cl_times = [nearest(control_wlt.time, control_events.syllable_info{t}.onset(s))...
                    nearest(control_wlt.time,control_events.syllable_info{t}.offset(s))];
                corr_rspn_fft_syll{scorr_ind,f} = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)));
                scorr_ind = scorr_ind+1;
            elseif strcmp(control_events.word_info{t}.response{control_events.syllable_info{t}.word_id(s)},'0')
                % find closest timepoint
                cl_times = [nearest(control_wlt.time, control_events.syllable_info{t}.onset(s))...
                    nearest(control_wlt.time,control_events.syllable_info{t}.offset(s))];
                no_rspn_fft_syll{sno_ind,f} = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)));
                sno_ind = sno_ind+1;
                
            else % wrong responses
                % find closest timepoint
                cl_times = [nearest(control_wlt.time, control_events.syllable_info{t}.onset(s))...
                    nearest(control_wlt.time,control_events.syllable_info{t}.offset(s))];
                wrng_rspn_fft_syll{swrng_ind,f} = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)));
                swrng_ind = swrng_ind+1;
                
            end
        end
        
        for p = 1:size(control_events.all_info{t},1)
            % Separate phonemes like different trials
            % Check if correct
            % Collect fft in a cell structure
            if strcmp(control_events.word_info{t}.response{control_events.all_info{t}.word_id(p)},'1')
                % find closest timepoint
                cl_times = [nearest(control_wlt.time, control_events.syllable_info{t}.onset(s))...
                    nearest(control_wlt.time,control_events.syllable_info{t}.offset(s))];
                corr_rspn_fft_pho{pcorr_ind,f} = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)));
                pcorr_ind = pcorr_ind+1;
            elseif strcmp(control_events.word_info{t}.response{control_events.all_info{t}.word_id(p)},'0')
                % find closest timepoint
                cl_times = [nearest(control_wlt.time, control_events.syllable_info{t}.onset(s))...
                    nearest(control_wlt.time,control_events.syllable_info{t}.offset(s))];
                no_rspn_fft_pho{pno_ind,f} = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)));
                pno_ind = pno_ind+1;
                
            else % wrong responses
                % find closest timepoint
                cl_times = [nearest(control_wlt.time, control_events.all_info{t}.onset(s))...
                    nearest(control_wlt.time,control_events.syllable_info{t}.offset(s))];
                wrng_rspn_fft_pho{pwrng_ind,f} = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)));
                pwrng_ind = pwrng_ind+1;
                
            end
        end
    end
end

fft_of_words = [];

fft_of_words.freq_bands = freq_bands;
fft_of_words.freq_band_dtls{1} = control_wlt.freq(nearest(control_wlt.freq, freq_bands{1}(1)):...
    nearest(control_wlt.freq,freq_bands{1}(2)));
fft_of_words.freq_band_dtls{2} = control_wlt.freq(nearest(control_wlt.freq, freq_bands{2}(1)):...
    nearest(control_wlt.freq,freq_bands{2}(2)));
fft_of_words.freq_band_dtls{3} = control_wlt.freq(nearest(control_wlt.freq, freq_bands{3}(1)):...
    nearest(control_wlt.freq,freq_bands{3}(2)));

fft_of_words.corr_rspn_fft_word = corr_rspn_fft_word;
fft_of_words.no_rspn_fft_word = no_rspn_fft_word;
fft_of_words.wrng_rspn_fft_word = wrng_rspn_fft_word;

fft_of_words.corr_rspn_fft_syll = corr_rspn_fft_syll;
fft_of_words.no_rspn_fft_syll = no_rspn_fft_syll;
fft_of_words.wrng_rspn_fft_syll = wrng_rspn_fft_syll;

fft_of_words.corr_rspn_fft_pho = corr_rspn_fft_pho;
fft_of_words.no_rspn_fft_pho = no_rspn_fft_pho;
fft_of_words.wrng_rspn_fft_pho = wrng_rspn_fft_pho;

if ~exist(fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_')),'dir'),mkdir(fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'))),end
save(fullfile(Sbj_Metadata.results,strjoin(control_blocks,'_'), [strjoin(control_blocks,'_') '_ctrl_word_fft.mat']),'fft_of_words','-v7.3');

%% Compare angles of whole segment
comp_angle_res{size(corr_rspn_fft_word{1},1),size(corr_rspn_fft_word,2)} = [];

for el = 81:96 %1:size(corr_rspn_fft_word{1},1) % loop electrode
    clear noncorr_angles corr_angles
    %% First work on word level
    figure('Units','normalized','Position', [0 0  .5 1]);
    for f = 1:size(corr_rspn_fft_word,2) % loop frequency band
        % average frequency and time dimension per trial
        for i=1:size(corr_rspn_fft_word,1)
            corr_angles(i) = squeeze(mean(angle(corr_rspn_fft_word{i,f}(el,:,:)),[2 3]));
        end
        
        for i=1:size(no_rspn_fft_word,1)
            noncorr_angles(i) = squeeze(mean(angle(no_rspn_fft_word{i,f}(el,:,:)),[2 3]));
        end
        for i=1:size(wrng_rspn_fft_word,1)
            noncorr_angles(end+1) = squeeze(mean(angle(wrng_rspn_fft_word{i,f}(el,:,:)),[2 3]));
        end
        
        % Test for difference in angles
        [pval, table] = circ_wwtest(corr_angles', noncorr_angles');
        comp_angle_res{el,f}.wwtest_word_pval = pval;
        comp_angle_res{el,f}.wwtest_word_table = table;
        
        % plot rose plot
        subplot(size(corr_rspn_fft_word,2),2,2*f-1)
        circ_plot(corr_angles,'hist',[],90,true,true,'linewidth',2,'color','r');
        hold on
        ylabel([num2str(freq_bands{f}(1)),'-',num2str(freq_bands{f}(2)), ' Hz range'],'FontSize',15,'FontWeight','bold')
        %         text(0,0.5,[num2str(freq_bands{f}(1)),'-',num2str(freq_bands{f}(2)), ' Hz range'],'FontSize',15,'FontWeight','bold','Units','normalized','HorizontalAlignment', 'left', ...
        %             'VerticalAlignment', 'middle','BackgroundColor', 'white','Rotation',90)
        
        if pval > 0.05
            title(['Correct vs Non-correct (wrong+unheard) word report circular stats p-value: ' num2str(pval)],'HorizontalAlignment','left')
        else
            title(['Correct vs Non-correct (wrong+unheard) word report circular stats p-value: ' num2str(pval)],'Color','r','HorizontalAlignment','left')
        end
        subplot(size(corr_rspn_fft_word,2),2,2*f)
        circ_plot(noncorr_angles,'hist',[],90,true,true,'linewidth',2,'color','r');
        
    end
    
    sgtitle(['Electrode label: ' control_wlt.label{el}])
    print(fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),[control_wlt.label{el} , 'angle_comp_words.jpeg']),'-djpeg','-r300')
    close all
    
    clear noncorr_angles corr_angles
    %% Second work on Syllables
    figure('Units','normalized','Position', [0 0  .5 1]);
    for f = 1:size(corr_rspn_fft_syll,2) % loop frequency band
        % average frequency and time dimension per trial
        for i=1:size(corr_rspn_fft_syll,1)
            corr_angles(i) = squeeze(mean(angle(corr_rspn_fft_syll{i,f}(el,:,:)),[2 3]));
        end
        
        for i=1:size(no_rspn_fft_syll,1)
            noncorr_angles(i) = squeeze(mean(angle(no_rspn_fft_syll{i,f}(el,:,:)),[2 3]));
        end
        for i=1:size(wrng_rspn_fft_syll,1)
            
            noncorr_angles(end+1) = squeeze(mean(angle(wrng_rspn_fft_syll{i,f}(el,:,:)),[2 3]));
        end
        
        % Test for difference in angles
        [pval, table] = circ_wwtest(corr_angles', noncorr_angles');
        
        comp_angle_res{el,f}.wwtest_syll_pval = pval;
        comp_angle_res{el,f}.wwtest_syll_table = table;
        
        % plot rose plot
        subplot(size(corr_rspn_fft_syll,2),2,2*f-1)
        circ_plot(corr_angles,'hist',[],90,true,true,'linewidth',2,'color','r');
        hold on
        ylabel([num2str(freq_bands{f}(1)),'-',num2str(freq_bands{f}(2)), ' Hz range'],'FontSize',15,'FontWeight','bold')
        %         text(0,0.5,[num2str(freq_bands{f}(1)),'-',num2str(freq_bands{f}(2)), ' Hz range'],'FontSize',15,'FontWeight','bold','Units','normalized','HorizontalAlignment', 'left', ...
        %             'VerticalAlignment', 'middle','BackgroundColor', 'white','Rotation',90)
        
        if pval > 0.05
            title(['Correct vs Non-correct (wrong+unheard) syllable report circular stats p-value: ' num2str(pval)],'HorizontalAlignment','left')
        else
            title(['Correct vs Non-correct (wrong+unheard) syllable report circular stats p-value: ' num2str(pval)],'Color','r','HorizontalAlignment','left')
        end
        subplot(size(corr_rspn_fft_syll,2),2,2*f)
        circ_plot(noncorr_angles,'hist',[],90,true,true,'linewidth',2,'color','r');
        
    end
    
    sgtitle(['Electrode label: ' control_wlt.label{el}])
    print(fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),[control_wlt.label{el} , 'angle_comp_syll.jpeg']),'-djpeg','-r300')
    close all
    clear noncorr_angles corr_angles
    
    %% Third: work on phonemes
    figure('Units','normalized','Position', [0 0  .5 1]);
    for f = 1:size(corr_rspn_fft_pho,2) % loop frequency band
        % average frequency and time dimension per trial
        for i=1:size(corr_rspn_fft_pho,1)
            corr_angles(i) = squeeze(mean(angle(corr_rspn_fft_pho{i,f}(el,:,:)),[2 3]));
        end
        
        for i=1:size(no_rspn_fft_pho,1)
            noncorr_angles(i) = squeeze(mean(angle(no_rspn_fft_pho{i,f}(el,:,:)),[2 3]));
        end
        for i=1:size(wrng_rspn_fft_pho,1)
            
            noncorr_angles(end+1) = squeeze(mean(angle(wrng_rspn_fft_pho{i,f}(el,:,:)),[2 3]));
        end
        
        % Test for difference in angles
        [pval, table] = circ_wwtest(corr_angles', noncorr_angles');
        
        comp_angle_res{el,f}.wwtest_pho_pval = pval;
        comp_angle_res{el,f}.wwtest_pho_table = table;
        
        % plot rose plot
        subplot(size(corr_rspn_fft_pho,2),2,2*f-1)
        circ_plot(corr_angles,'hist',[],90,true,true,'linewidth',2,'color','r');
        hold on
        ylabel([num2str(freq_bands{f}(1)),'-',num2str(freq_bands{f}(2)), ' Hz range'],'FontSize',15,'FontWeight','bold')
        %         text(0,0.5,[num2str(freq_bands{f}(1)),'-',num2str(freq_bands{f}(2)), ' Hz range'],'FontSize',15,'FontWeight','bold','Units','normalized','HorizontalAlignment', 'left', ...
        %             'VerticalAlignment', 'middle','BackgroundColor', 'white','Rotation',90)
        
        if pval > 0.05
            title(['Correct vs Non-correct (wrong+unheard) phoneme report circular stats p-value: ' num2str(pval)],'HorizontalAlignment','left')
        else
            title(['Correct vs Non-correct (wrong+unheard) phoneme report circular stats p-value: ' num2str(pval)],'Color','r','HorizontalAlignment','left')
        end
        subplot(size(corr_rspn_fft_pho,2),2,2*f)
        circ_plot(noncorr_angles,'hist',[],90,true,true,'linewidth',2,'color','r');
        
    end
    
    sgtitle(['Electrode label: ' control_wlt.label{el}])
    print(fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),[control_wlt.label{el} , 'angle_comp_phonemes.jpeg']),'-djpeg','-r300')
    close all
    
end

save(fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),[control_wlt.label{el} , 'angle_comp_stats.mat']),'comp_angle_res')


%% Compare angles of 20 ms response
clear comp_angle_res
comp_angle_res{size(corr_rspn_fft_word{1},1),size(corr_rspn_fft_word,2)} = [];

for el = 81:96 %1:size(corr_rspn_fft_word{1},1) % loop electrode
    clear noncorr_angles corr_angles
    %% First work on word level
    figure('Units','normalized','Position', [0 0  .5 1]);
    for f = 1:size(corr_rspn_fft_word,2) % loop frequency band
        % average frequency and time dimension per trial
        for i=1:size(corr_rspn_fft_word,1)
            corr_angles(i) = squeeze(mean(angle(corr_rspn_fft_word{i,f}(el,:,3)),[2 3]));
        end
        
        for i=1:size(no_rspn_fft_word,1)
            noncorr_angles(i) = squeeze(mean(angle(no_rspn_fft_word{i,f}(el,:,3)),[2 3]));
        end
        for i=1:size(wrng_rspn_fft_word,1)
            noncorr_angles(end+1) = squeeze(mean(angle(wrng_rspn_fft_word{i,f}(el,:,3)),[2 3]));
        end
        
        % Test for difference in angles
        [pval, table] = circ_wwtest(corr_angles', noncorr_angles');
        comp_angle_res{el,f}.wwtest_word_pval = pval;
        comp_angle_res{el,f}.wwtest_word_table = table;
        
        % plot rose plot
        subplot(size(corr_rspn_fft_word,2),2,2*f-1)
        circ_plot(corr_angles,'hist',[],90,true,true,'linewidth',2,'color','r');
        hold on
        ylabel([num2str(freq_bands{f}(1)),'-',num2str(freq_bands{f}(2)), ' Hz range'],'FontSize',15,'FontWeight','bold')
        %         text(0,0.5,[num2str(freq_bands{f}(1)),'-',num2str(freq_bands{f}(2)), ' Hz range'],'FontSize',15,'FontWeight','bold','Units','normalized','HorizontalAlignment', 'left', ...
        %             'VerticalAlignment', 'middle','BackgroundColor', 'white','Rotation',90)
        
        if pval > 0.05
            title(['Correct vs Non-correct (wrong+unheard) word report (20ms) circular stats p-value: ' num2str(pval)],'HorizontalAlignment','left')
        else
            title(['Correct vs Non-correct (wrong+unheard) word report (20ms) circular stats p-value: ' num2str(pval)],'Color','r','HorizontalAlignment','left')
        end
        subplot(size(corr_rspn_fft_word,2),2,2*f)
        circ_plot(noncorr_angles,'hist',[],90,true,true,'linewidth',2,'color','r');
        
    end
    
    sgtitle(['Electrode label: ' control_wlt.label{el}])
    print(fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),[control_wlt.label{el} , 'angle_comp_words_20ms.jpeg']),'-djpeg','-r300')
    close all
    
    clear noncorr_angles corr_angles
    %% Second work on Syllables
    figure('Units','normalized','Position', [0 0  .5 1]);
    for f = 1:size(corr_rspn_fft_syll,2) % loop frequency band
        % average frequency and time dimension per trial
        for i=1:size(corr_rspn_fft_syll,1)
            corr_angles(i) = squeeze(mean(angle(corr_rspn_fft_syll{i,f}(el,:,3)),[2 3]));
        end
        
        for i=1:size(no_rspn_fft_syll,1)
            noncorr_angles(i) = squeeze(mean(angle(no_rspn_fft_syll{i,f}(el,:,3)),[2 3]));
        end
        for i=1:size(wrng_rspn_fft_syll,1)
            
            noncorr_angles(end+1) = squeeze(mean(angle(wrng_rspn_fft_syll{i,f}(el,:,3)),[2 3]));
        end
        
        % Test for difference in angles
        [pval, table] = circ_wwtest(corr_angles', noncorr_angles');
        
        comp_angle_res{el,f}.wwtest_syll_pval = pval;
        comp_angle_res{el,f}.wwtest_syll_table = table;
        
        % plot rose plot
        subplot(size(corr_rspn_fft_syll,2),2,2*f-1)
        circ_plot(corr_angles,'hist',[],90,true,true,'linewidth',2,'color','r');
        hold on
        ylabel([num2str(freq_bands{f}(1)),'-',num2str(freq_bands{f}(2)), ' Hz range'],'FontSize',15,'FontWeight','bold')
        %         text(0,0.5,[num2str(freq_bands{f}(1)),'-',num2str(freq_bands{f}(2)), ' Hz range'],'FontSize',15,'FontWeight','bold','Units','normalized','HorizontalAlignment', 'left', ...
        %             'VerticalAlignment', 'middle','BackgroundColor', 'white','Rotation',90)
        
        if pval > 0.05
            title(['Correct vs Non-correct (wrong+unheard) syllable (20ms) report circular stats p-value: ' num2str(pval)],'HorizontalAlignment','left')
        else
            title(['Correct vs Non-correct (wrong+unheard) syllable (20ms) report circular stats p-value: ' num2str(pval)],'Color','r','HorizontalAlignment','left')
        end
        subplot(size(corr_rspn_fft_syll,2),2,2*f)
        circ_plot(noncorr_angles,'hist',[],90,true,true,'linewidth',2,'color','r');
        
    end
    
    sgtitle(['Electrode label: ' control_wlt.label{el}])
    print(fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),[control_wlt.label{el} , 'angle_comp_syll_20ms.jpeg']),'-djpeg','-r300')
    close all
    clear noncorr_angles corr_angles
    
    %% Third: work on phonemes
    figure('Units','normalized','Position', [0 0  .5 1]);
    for f = 1:size(corr_rspn_fft_pho,2) % loop frequency band
        % average frequency and time dimension per trial
        for i=1:size(corr_rspn_fft_pho,1)
            corr_angles(i) = squeeze(mean(angle(corr_rspn_fft_pho{i,f}(el,:,3)),[2 3]));
        end
        
        for i=1:size(no_rspn_fft_pho,1)
            noncorr_angles(i) = squeeze(mean(angle(no_rspn_fft_pho{i,f}(el,:,3)),[2 3]));
        end
        for i=1:size(wrng_rspn_fft_pho,1)
            
            noncorr_angles(end+1) = squeeze(mean(angle(wrng_rspn_fft_pho{i,f}(el,:,3)),[2 3]));
        end
        
        % Test for difference in angles
        [pval, table] = circ_wwtest(corr_angles', noncorr_angles');
        
        comp_angle_res{el,f}.wwtest_pho_pval = pval;
        comp_angle_res{el,f}.wwtest_pho_table = table;
        
        % plot rose plot
        subplot(size(corr_rspn_fft_pho,2),2,2*f-1)
        circ_plot(corr_angles,'hist',[],90,true,true,'linewidth',2,'color','r');
        hold on
        ylabel([num2str(freq_bands{f}(1)),'-',num2str(freq_bands{f}(2)), ' Hz range'],'FontSize',15,'FontWeight','bold')
        %         text(0,0.5,[num2str(freq_bands{f}(1)),'-',num2str(freq_bands{f}(2)), ' Hz range'],'FontSize',15,'FontWeight','bold','Units','normalized','HorizontalAlignment', 'left', ...
        %             'VerticalAlignment', 'middle','BackgroundColor', 'white','Rotation',90)
        
        if pval > 0.05
            title(['Correct vs Non-correct (wrong+unheard) phoneme (20ms) report circular stats p-value: ' num2str(pval)],'HorizontalAlignment','left')
        else
            title(['Correct vs Non-correct (wrong+unheard) phoneme (20ms) report circular stats p-value: ' num2str(pval)],'Color','r','HorizontalAlignment','left')
        end
        subplot(size(corr_rspn_fft_pho,2),2,2*f)
        circ_plot(noncorr_angles,'hist',[],90,true,true,'linewidth',2,'color','r');
        
    end
    
    sgtitle(['Electrode label: ' control_wlt.label{el}])
    print(fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),[control_wlt.label{el} , 'angle_comp_phonemes_20ms.jpeg']),'-djpeg','-r300')
    close all
    
end

save(fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),[control_wlt.label{el} , 'angle_comp_stats_20ms.mat']),'comp_angle_res')

%% Compare angles of 0 ms response
clear comp_angle_res
comp_angle_res{size(corr_rspn_fft_word{1},1),size(corr_rspn_fft_word,2)} = [];

for el = 81:96 %1:size(corr_rspn_fft_word{1},1) % loop electrode
    clear noncorr_angles corr_angles
    %% First work on word level
    figure('Units','normalized','Position', [0 0  .5 1]);
    for f = 1:size(corr_rspn_fft_word,2) % loop frequency band
        % average frequency and time dimension per trial
        for i=1:size(corr_rspn_fft_word,1)
            corr_angles(i) = squeeze(mean(angle(corr_rspn_fft_word{i,f}(el,:,1)),[2 3]));
        end
        
        for i=1:size(no_rspn_fft_word,1)
            noncorr_angles(i) = squeeze(mean(angle(no_rspn_fft_word{i,f}(el,:,1)),[2 3]));
        end
        for i=1:size(wrng_rspn_fft_word,1)
            noncorr_angles(end+1) = squeeze(mean(angle(wrng_rspn_fft_word{i,f}(el,:,1)),[2 3]));
        end
        
        % Test for difference in angles
        [pval, table] = circ_wwtest(corr_angles', noncorr_angles');
        comp_angle_res{el,f}.wwtest_word_pval = pval;
        comp_angle_res{el,f}.wwtest_word_table = table;
        
        % plot rose plot
        subplot(size(corr_rspn_fft_word,2),2,2*f-1)
        circ_plot(corr_angles,'hist',[],90,true,true,'linewidth',2,'color','r');
        hold on
        ylabel([num2str(freq_bands{f}(1)),'-',num2str(freq_bands{f}(2)), ' Hz range'],'FontSize',15,'FontWeight','bold')
        %         text(0,0.5,[num2str(freq_bands{f}(1)),'-',num2str(freq_bands{f}(2)), ' Hz range'],'FontSize',15,'FontWeight','bold','Units','normalized','HorizontalAlignment', 'left', ...
        %             'VerticalAlignment', 'middle','BackgroundColor', 'white','Rotation',90)
        
        if pval > 0.05
            title(['Correct vs Non-correct (wrong+unheard) word report (0ms) circular stats p-value: ' num2str(pval)],'HorizontalAlignment','left')
        else
            title(['Correct vs Non-correct (wrong+unheard) word report (0ms) circular stats p-value: ' num2str(pval)],'Color','r','HorizontalAlignment','left')
        end
        subplot(size(corr_rspn_fft_word,2),2,2*f)
        circ_plot(noncorr_angles,'hist',[],90,true,true,'linewidth',2,'color','r');
        
    end
    
    sgtitle(['Electrode label: ' control_wlt.label{el}])
    print(fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),[control_wlt.label{el} , 'angle_comp_words_0ms.jpeg']),'-djpeg','-r300')
    close all
    
    clear noncorr_angles corr_angles
    %% Second work on Syllables
    figure('Units','normalized','Position', [0 0  .5 1]);
    for f = 1:size(corr_rspn_fft_syll,2) % loop frequency band
        % average frequency and time dimension per trial
        for i=1:size(corr_rspn_fft_syll,1)
            corr_angles(i) = squeeze(mean(angle(corr_rspn_fft_syll{i,f}(el,:,1)),[2 3]));
        end
        
        for i=1:size(no_rspn_fft_syll,1)
            noncorr_angles(i) = squeeze(mean(angle(no_rspn_fft_syll{i,f}(el,:,1)),[2 3]));
        end
        for i=1:size(wrng_rspn_fft_syll,1)
            
            noncorr_angles(end+1) = squeeze(mean(angle(wrng_rspn_fft_syll{i,f}(el,:,1)),[2 3]));
        end
        
        % Test for difference in angles
        [pval, table] = circ_wwtest(corr_angles', noncorr_angles');
        
        comp_angle_res{el,f}.wwtest_syll_pval = pval;
        comp_angle_res{el,f}.wwtest_syll_table = table;
        
        % plot rose plot
        subplot(size(corr_rspn_fft_syll,2),2,2*f-1)
        circ_plot(corr_angles,'hist',[],90,true,true,'linewidth',2,'color','r');
        hold on
        ylabel([num2str(freq_bands{f}(1)),'-',num2str(freq_bands{f}(2)), ' Hz range'],'FontSize',15,'FontWeight','bold')
        %         text(0,0.5,[num2str(freq_bands{f}(1)),'-',num2str(freq_bands{f}(2)), ' Hz range'],'FontSize',15,'FontWeight','bold','Units','normalized','HorizontalAlignment', 'left', ...
        %             'VerticalAlignment', 'middle','BackgroundColor', 'white','Rotation',90)
        
        if pval > 0.05
            title(['Correct vs Non-correct (wrong+unheard) syllable (0ms) report circular stats p-value: ' num2str(pval)],'HorizontalAlignment','left')
        else
            title(['Correct vs Non-correct (wrong+unheard) syllable (0ms) report circular stats p-value: ' num2str(pval)],'Color','r','HorizontalAlignment','left')
        end
        subplot(size(corr_rspn_fft_syll,2),2,2*f)
        circ_plot(noncorr_angles,'hist',[],90,true,true,'linewidth',2,'color','r');
        
    end
    
    sgtitle(['Electrode label: ' control_wlt.label{el}])
    print(fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),[control_wlt.label{el} , 'angle_comp_syll_0ms.jpeg']),'-djpeg','-r300')
    close all
    clear noncorr_angles corr_angles
    
    %% Third: work on phonemes
    figure('Units','normalized','Position', [0 0  .5 1]);
    for f = 1:size(corr_rspn_fft_pho,2) % loop frequency band
        % average frequency and time dimension per trial
        for i=1:size(corr_rspn_fft_pho,1)
            corr_angles(i) = squeeze(mean(angle(corr_rspn_fft_pho{i,f}(el,:,1)),[2 3]));
        end
        
        for i=1:size(no_rspn_fft_pho,1)
            noncorr_angles(i) = squeeze(mean(angle(no_rspn_fft_pho{i,f}(el,:,1)),[2 3]));
        end
        for i=1:size(wrng_rspn_fft_pho,1)
            
            noncorr_angles(end+1) = squeeze(mean(angle(wrng_rspn_fft_pho{i,f}(el,:,1)),[2 3]));
        end
        
        % Test for difference in angles
        [pval, table] = circ_wwtest(corr_angles', noncorr_angles');
        
        comp_angle_res{el,f}.wwtest_pho_pval = pval;
        comp_angle_res{el,f}.wwtest_pho_table = table;
        
        % plot rose plot
        subplot(size(corr_rspn_fft_pho,2),2,2*f-1)
        circ_plot(corr_angles,'hist',[],90,true,true,'linewidth',2,'color','r');
        hold on
        ylabel([num2str(freq_bands{f}(1)),'-',num2str(freq_bands{f}(2)), ' Hz range'],'FontSize',15,'FontWeight','bold')
        %         text(0,0.5,[num2str(freq_bands{f}(1)),'-',num2str(freq_bands{f}(2)), ' Hz range'],'FontSize',15,'FontWeight','bold','Units','normalized','HorizontalAlignment', 'left', ...
        %             'VerticalAlignment', 'middle','BackgroundColor', 'white','Rotation',90)
        
        if pval > 0.05
            title(['Correct vs Non-correct (wrong+unheard) phoneme (0ms) report circular stats p-value: ' num2str(pval)],'HorizontalAlignment','left')
        else
            title(['Correct vs Non-correct (wrong+unheard) phoneme (0ms) report circular stats p-value: ' num2str(pval)],'Color','r','HorizontalAlignment','left')
        end
        subplot(size(corr_rspn_fft_pho,2),2,2*f)
        circ_plot(noncorr_angles,'hist',[],90,true,true,'linewidth',2,'color','r');
        
    end
    
    sgtitle(['Electrode label: ' control_wlt.label{el}])
    print(fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),[control_wlt.label{el} , 'angle_comp_phonemes_0ms.jpeg']),'-djpeg','-r300')
    close all
    
end

save(fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),[control_wlt.label{el} , 'angle_comp_stats_0ms.mat']),'comp_angle_res')
