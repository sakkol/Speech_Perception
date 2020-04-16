function SP_word_evoked(Sbj_Metadata,control_blocks)
% This is to get the mean of evoked responses (ERP and HFA) word heard
% vs others.

%% Select blocks to import
vars=who;
if ~ismember(vars,'control_blocks')
    control_blocks = select_cont_blocks(Sbj_Metadata);
end
clear vars

%% bring in these blocks and combine only the control events
fprintf('These blocks are going to be used: %s\n',strjoin(control_blocks,', '))
for b = 1:length(control_blocks)
    curr_block = control_blocks{b};
    fprintf('...loading:%s\n',curr_block)
    
    % Load iEEG
    load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_wlt.mat']))
    load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']))
    events = info.events;
    
    % Select only control events
    control_idx = events.Cond_code == 1;
    
    % Select ecog data
    cfg = [];
    cfg.trials = control_idx;
    [epoched_wlt.wlt] = ft_selectdata(cfg, epoched_wlt.wlt);
    curr_ERP = ft_selectdata(cfg, epoched_wlt);
    
    % from fourierspectrum to powerspectrum
    cfg = [];
    cfg.output='abs';
    cfg.keeptrials = 'yes';
    epoched_wlt.wlt=ft_freqdescriptives(cfg,epoched_wlt.wlt);
    
    % Select events
    events = events(control_idx,:);
    
    % Append to overall list
    if b == 1
        control_wlt = epoched_wlt.wlt;
        control_events = events;
        control_ERP = curr_ERP;
    else
        cfg = [];
        cfg.parameter  = 'powspctrm';
        control_wlt = ft_appendfreq(cfg, control_wlt, epoched_wlt.wlt);
        cfg = [];
        control_ERP = ft_appendtimelock(cfg,control_ERP,curr_ERP);
        control_events = [control_events;events];
    end
    clear epoched_wlt events info curr_ERP control_idx
end

%%  Baseline correct time-freq data
cfg              = [];
cfg.baseline     = [-3.45 -3.05]; % seconds (prespeech part is 3.5455 seconds) [0.5sec + 3.0455sec]
cfg.baselinetype = 'db';
cfg.parameter    = 'powspctrm';
[control_wlt]         = ft_freqbaseline(cfg, control_wlt);

cfg              = [];
cfg.baseline     = [-3.45 -3.05];
cfg.channel      = 'all';
cfg.parameter    = 'trial';
control_ERP      = ft_timelockbaseline(cfg, control_ERP);

%% Separate conditions

fprintf('There are total of %d events\n',size(control_events,1))

wcorr_ind = 1;
wno_ind = 1;
wwrng_ind = 1;
scorr_ind = 1;
sno_ind = 1;
swrng_ind = 1;
pcorr_ind = 1;
pno_ind = 1;
pwrng_ind = 1;
fprintf('\nCurrent event:')
cl_freqs = [nearest(control_wlt.freq, 70), ...
    nearest(control_wlt.freq,150)];
for t = 1:size(control_events,1)
    
    fprintf('-%d',t)
    for w = 1:5
        % Separate words
        % Check if correct
        % Collect fouri in a cell structure
        if strcmp(control_events.word_info{t}.response{w},'1')
            % find closest timepoint
            cl_times = [nearest(control_wlt.time, control_events.word_info{t}.onset(w))...
                nearest(control_wlt.time,control_events.word_info{t}.offset(w))];
            corr_rspn_powspec_word{wcorr_ind,1} = squeeze(mean(control_wlt.powspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)),3));
            cl_times = [nearest(control_ERP.time, control_events.word_info{t}.onset(w))...
                nearest(control_ERP.time,control_events.word_info{t}.offset(w))];
            corr_rspn_ERP_word{wcorr_ind,1} = squeeze(control_ERP.trial(t,:,cl_times(1):cl_times(2)));
            wcorr_ind = wcorr_ind+1;
        elseif strcmp(control_events.word_info{t}.response{w},'0')
            % find closest timepoint
            cl_times = [nearest(control_wlt.time, control_events.word_info{t}.onset(w))...
                nearest(control_wlt.time,control_events.word_info{t}.offset(w))];
            no_rspn_powspec_word{wno_ind,1} =  squeeze(mean(control_wlt.powspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)),3));
            cl_times = [nearest(control_ERP.time, control_events.word_info{t}.onset(w))...
                nearest(control_ERP.time,control_events.word_info{t}.offset(w))];
            no_rspn_ERP_word{wno_ind,1} = squeeze(control_ERP.trial(t,:,cl_times(1):cl_times(2)));
            wno_ind = wno_ind+1;
            
        else % wrong responses
            % find closest timepoint
            cl_times = [nearest(control_wlt.time, control_events.word_info{t}.onset(w))...
                nearest(control_wlt.time,control_events.word_info{t}.offset(w))];
            wrng_rspn_powspec_word{wwrng_ind,1} =  squeeze(mean(control_wlt.powspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)),3));
            cl_times = [nearest(control_ERP.time, control_events.word_info{t}.onset(w))...
                nearest(control_ERP.time,control_events.word_info{t}.offset(w))];
            wrng_rspn_ERP_word{wwrng_ind,1} = squeeze(control_ERP.trial(t,:,cl_times(1):cl_times(2)));
            wwrng_ind = wwrng_ind+1;
            
        end
        
    end
    
    for s = 1:size(control_events.syllable_info{t},1)
        % Separate syllables like different trials
        % Check if correct
        % Collect fouri in a cell structure
        if strcmp(control_events.word_info{t}.response{control_events.syllable_info{t}.word_id(s)},'1')
            % find closest timepoint
            cl_times = [nearest(control_wlt.time, control_events.syllable_info{t}.onset(s))...
                nearest(control_wlt.time,control_events.syllable_info{t}.offset(s))];
            corr_rspn_powspec_syll{scorr_ind,1} = squeeze(mean(control_wlt.powspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)),3));
            cl_times = [nearest(control_ERP.time, control_events.syllable_info{t}.onset(s))...
                nearest(control_ERP.time,control_events.syllable_info{t}.offset(s))];
            corr_rspn_ERP_syll{scorr_ind,1} = squeeze(control_ERP.trial(t,:,cl_times(1):cl_times(2)));
            scorr_ind = scorr_ind+1;
        elseif strcmp(control_events.word_info{t}.response{control_events.syllable_info{t}.word_id(s)},'0')
            % find closest timepoint
            cl_times = [nearest(control_wlt.time, control_events.syllable_info{t}.onset(s))...
                nearest(control_wlt.time,control_events.syllable_info{t}.offset(s))];
            no_rspn_powspec_syll{sno_ind,1} = squeeze(mean(control_wlt.powspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)),3));
            cl_times = [nearest(control_ERP.time, control_events.syllable_info{t}.onset(s))...
                nearest(control_ERP.time,control_events.syllable_info{t}.offset(s))];
            no_rspn_ERP_syll{sno_ind,1} = squeeze(control_ERP.trial(t,:,cl_times(1):cl_times(2)));
            sno_ind = sno_ind+1;
            
        else % wrong responses
            % find closest timepoint
            cl_times = [nearest(control_wlt.time, control_events.syllable_info{t}.onset(s))...
                nearest(control_wlt.time,control_events.syllable_info{t}.offset(s))];
            wrng_rspn_powspec_syll{swrng_ind,1} = squeeze(mean(control_wlt.powspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)),3));
            cl_times = [nearest(control_ERP.time, control_events.syllable_info{t}.onset(s))...
                nearest(control_ERP.time,control_events.syllable_info{t}.offset(s))];
            wrng_rspn_ERP_syll{swrng_ind,1} = squeeze(control_ERP.trial(t,:,cl_times(1):cl_times(2)));
            swrng_ind = swrng_ind+1;
            
        end
    end
    
    for p = 1:size(control_events.all_info{t},1)
        % Separate phonemes like different trials
        % Check if correct
        % Collect fouri in a cell structure
        if strcmp(control_events.word_info{t}.response{control_events.all_info{t}.word_id(p)},'1')
            % find closest timepoint
            cl_times = [nearest(control_wlt.time, control_events.all_info{t}.onset(p))...
                nearest(control_wlt.time,control_events.all_info{t}.offset(p))];
            corr_rspn_powspec_pho{pcorr_ind,1} = squeeze(mean(control_wlt.powspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)),3));
            cl_times = [nearest(control_ERP.time, control_events.all_info{t}.onset(p))...
                nearest(control_ERP.time,control_events.all_info{t}.offset(p))];
            corr_rspn_ERP_pho{pcorr_ind,1} = squeeze(control_ERP.trial(t,:,cl_times(1):cl_times(2)));
            pcorr_ind = pcorr_ind+1;
        elseif strcmp(control_events.word_info{t}.response{control_events.all_info{t}.word_id(p)},'0')
            % find closest timepoint
            cl_times = [nearest(control_wlt.time, control_events.all_info{t}.onset(p))...
                nearest(control_wlt.time,control_events.all_info{t}.offset(p))];
            no_rspn_powspec_pho{pno_ind,1} = squeeze(mean(control_wlt.powspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)),3));
            cl_times = [nearest(control_ERP.time, control_events.all_info{t}.onset(p))...
                nearest(control_ERP.time,control_events.all_info{t}.offset(p))];
            no_rspn_ERP_pho{pno_ind,1} = squeeze(control_ERP.trial(t,:,cl_times(1):cl_times(2)));
            pno_ind = pno_ind+1;
            
        else % wrong responses
            % find closest timepoint
            cl_times = [nearest(control_wlt.time, control_events.all_info{t}.onset(p))...
                nearest(control_wlt.time,control_events.all_info{t}.offset(p))];
            wrng_rspn_powspec_pho{pwrng_ind,1} = squeeze(mean(control_wlt.powspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)),3));
            cl_times = [nearest(control_ERP.time, control_events.all_info{t}.onset(p))...
                nearest(control_ERP.time,control_events.all_info{t}.offset(p))];
            wrng_rspn_ERP_pho{pwrng_ind,1} = squeeze(control_ERP.trial(t,:,cl_times(1):cl_times(2)));
            pwrng_ind = pwrng_ind+1;
            
        end
    end
    
end
fprintf('\n')


%% Average word evoked activity in different conditions overlaying with single trials
save_folder = fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),'word_evoked');
if ~exist(save_folder,'dir'),mkdir(save_folder),end

for el = 1:size(corr_rspn_ERP_word{1},1)
    % 3x2 plot, one side HFA, one side ERP, 3 conditions (corr, no, wrng)
    figure('Units','normalized','Position', [0 0  .6 1]);
    
    % plot single trials and average of corr - ERP
    subplot(321)
    for c = 1:size(corr_rspn_ERP_word,1)
        plot(corr_rspn_ERP_word{c}(el,:),'Color',[0.9412,0.9412,0.9412])
        hold on
        tr_lenght(c,1) = size(corr_rspn_ERP_word{c},2);
    end
    for_avg=[];
    for c = 1:size(corr_rspn_ERP_word,1)
        for_avg(c,:) = [corr_rspn_ERP_word{c}(el,:),NaN(1,max(tr_lenght)-tr_lenght(c))];
    end
    shadedErrorBar(1:size(for_avg,2),nanmean(for_avg,1),nanstd(for_avg),'lineprops','r')
    title([num2str(size(corr_rspn_ERP_word,1)) ' correct words ERP (mean+SD)'])
    xlim([1 max(tr_lenght)])
    xticks([200 400 600])
    xticklabels({'0.2','0.4','0.6'})
    
    % plot single trials and average of no resp - ERP
    subplot(323)
    tr_lenght=[];
    for c = 1:size(no_rspn_ERP_word,1)
        plot(no_rspn_ERP_word{c}(el,:),'Color',[0.9412,0.9412,0.9412])
        hold on
        tr_lenght(c,1) = size(no_rspn_ERP_word{c},2);
    end
    for_avg=[];
    for c = 1:size(no_rspn_ERP_word,1)
        for_avg(c,:) = [no_rspn_ERP_word{c}(el,:),NaN(1,max(tr_lenght)-tr_lenght(c))];
    end
    shadedErrorBar(1:size(for_avg,2),nanmean(for_avg,1),nanstd(for_avg),'lineprops','r')
    title([num2str(size(no_rspn_ERP_word,1)) ' no response words ERP (mean+SD)'])
    xlim([1 max(tr_lenght)])
    xticks([200 400 600])
    xticklabels({'0.2','0.4','0.6'})
    
    % plot single trials and average of wrong resp - ERP
    subplot(325)
    tr_lenght=[];
    for c = 1:size(wrng_rspn_ERP_word,1)
        plot(wrng_rspn_ERP_word{c}(el,:),'Color',[0.9412,0.9412,0.9412])
        hold on
        tr_lenght(c,1) = size(wrng_rspn_ERP_word{c},2);
    end
    for_avg=[];
    for c = 1:size(wrng_rspn_ERP_word,1)
        for_avg(c,:) = [wrng_rspn_ERP_word{c}(el,:),NaN(1,max(tr_lenght)-tr_lenght(c))];
    end
    shadedErrorBar(1:size(for_avg,2),nanmean(for_avg,1),nanstd(for_avg),'lineprops','r')
    title([num2str(size(wrng_rspn_ERP_word,1)) ' wrong response words ERP (mean+SD)'])
    xlim([1 max(tr_lenght)])
    xticks([200 400 600])
    xticklabels({'0.2','0.4','0.6'})
    xlabel('Time (s)')
    
    % plot single trials and average of corr - HFA
    subplot(322)
    for c = 1:size(corr_rspn_powspec_word,1)
        plot(corr_rspn_powspec_word{c}(el,:),'Color',[0.9412,0.9412,0.9412])
        hold on
        tr_lenght(c,1) = size(corr_rspn_powspec_word{c},2);
    end
    for_avg=[];
    for c = 1:size(corr_rspn_powspec_word,1)
        for_avg(c,:) = [corr_rspn_powspec_word{c}(el,:),NaN(1,max(tr_lenght)-tr_lenght(c))];
    end
    shadedErrorBar(1:size(for_avg,2),nanmean(for_avg,1),nanstd(for_avg),'lineprops','r')
    title([num2str(size(corr_rspn_powspec_word,1)) ' correct words HFA (mean+SD)'])
    xlim([1 max(tr_lenght)])
    xticks([20 40 60])
    xticklabels({'0.2','0.4','0.6'})
    
    % plot single trials and average of no resp - HFA
    subplot(324)
    for c = 1:size(no_rspn_powspec_word,1)
        plot(no_rspn_powspec_word{c}(el,:),'Color',[0.9412,0.9412,0.9412])
        hold on
        tr_lenght(c,1) = size(no_rspn_powspec_word{c},2);
    end
    for_avg=[];
    for c = 1:size(no_rspn_powspec_word,1)
        for_avg(c,:) = [no_rspn_powspec_word{c}(el,:),NaN(1,max(tr_lenght)-tr_lenght(c))];
    end
    shadedErrorBar(1:size(for_avg,2),nanmean(for_avg,1),nanstd(for_avg),'lineprops','r')
    title([num2str(size(no_rspn_powspec_word,1)) ' no response words HFA (mean+SD)'])
    xlim([1 max(tr_lenght)])
    xticks([20 40 60])
    xticklabels({'0.2','0.4','0.6'})
    
    % plot single trials and average of wrong resp - HFA
    subplot(326)
    for c = 1:size(wrng_rspn_powspec_word,1)
        plot(wrng_rspn_powspec_word{c}(el,:),'Color',[0.9412,0.9412,0.9412])
        hold on
        tr_lenght(c,1) = size(wrng_rspn_powspec_word{c},2);
    end
    for_avg=[];
    for c = 1:size(wrng_rspn_powspec_word,1)
        for_avg(c,:) = [wrng_rspn_powspec_word{c}(el,:),NaN(1,max(tr_lenght)-tr_lenght(c))];
    end
    shadedErrorBar(1:size(for_avg,2),nanmean(for_avg,1),nanstd(for_avg),'lineprops','r')
    title([num2str(size(wrng_rspn_powspec_word,1)) ' wrong response words HFA (mean+SD)'])
    xlim([1 max(tr_lenght)])
    xticks([20 40 60])
    xticklabels({'0.2','0.4','0.6'})
    xlabel('Time (s)')
    
    sgtitle({[control_ERP.label{el}, ' - word onset locked ERP and HFA'];...
        'Baseline corrected to prespeech only noise part ([-3.45 -3.05]sec of sentence onset)';...
        ['from ' num2str(length(control_blocks)) ' blocks']})
    
    % Save the figure
    fprintf('\t-Saving electrode #%d-%s, out of %d\n',el,control_ERP.label{el},size(control_ERP.label,1))
    print(fullfile(save_folder,[control_wlt.label{el} , '_words_only.jpg']),'-djpeg','-r300')
    close all
    
end


end
