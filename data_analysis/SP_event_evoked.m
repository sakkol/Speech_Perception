function SP_event_evoked(Sbj_Metadata,control_blocks,runagain)
% This is to get the mean of evoked responses (ERP and HFA) word heard
% vs others.

%% Select blocks to import
vars=who;
if ~ismember(vars,'control_blocks')
    control_blocks = select_cont_blocks(Sbj_Metadata);
elseif isempty(control_blocks)
    control_blocks = select_cont_blocks(Sbj_Metadata);
end
if ~ismember(vars,'runagain')
    runagain = 1;
end
clear vars

%% bring in these blocks and combine only the control events
fprintf('These blocks are going to be used: %s\n',strjoin(control_blocks,', '))

save_folder = fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'));
if ~exist(save_folder,'dir'),mkdir(save_folder),end
% check if this has already been run
if ~exist(fullfile(save_folder,[strjoin(control_blocks,'_') '_control_wltERP.mat']),'file') || runagain
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
        [epoched_wlt] = ft_selectdata(cfg, epoched_wlt);
        curr_ERP = ft_selectdata(cfg, epoched_data);
        
        % from fourierspectrum to powerspectrum
        cfg = [];
        cfg.output='abs';
        cfg.keeptrials = 'yes';
        epoched_wlt=ft_freqdescriptives(cfg,epoched_wlt);
        
        % Select events
        events = events(control_idx,:);
        
        % Append to overall list
        if b == 1
            control_wlt = epoched_wlt;
            control_events = events;
            control_ERP = curr_ERP;
        else
            cfg = [];
            cfg.parameter  = 'powspctrm';
            control_wlt = ft_appendfreq(cfg, control_wlt, epoched_wlt);
            cfg = [];
            control_ERP = ft_appendtimelock(cfg,control_ERP,curr_ERP);
            control_events = [control_events;events];
        end
        clear epoched_wlt events info curr_ERP control_idx epoched_data
    end
    
    fprintf('Saving to:\n\t->%s\n',fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),[strjoin(control_blocks,'_') '_control_wltERP.mat']))
    save(fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),[strjoin(control_blocks,'_') '_control_wltERP.mat']),'control_events','control_ERP','control_wlt','-v7.3')
else
    fprintf('Loading from:\n\t->%s\n',fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),[strjoin(control_blocks,'_') '_control_wltERP.mat']))
    load(fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),[strjoin(control_blocks,'_') '_control_wltERP.mat']),'control_events','control_ERP','control_wlt')
    
end

%% Baseline correct TF data but not ERPs
%  Baseline correct time-freq data
cfg              = [];
cfg.baseline     = [-3.45 -3.05]; % seconds (prespeech part is 3.5455 seconds) [0.5sec + 3.0455sec]
cfg.baselinetype = 'db';
cfg.parameter    = 'powspctrm';
[control_wlt]    = ft_freqbaseline(cfg, control_wlt);


% cfg              = [];
% cfg.baseline     = [-3.45 -3.05];
% cfg.baselinetype = 'db';
% cfg.channel      = 'all';
% cfg.parameter    = 'trial';
% control_ERP      = ft_timelockbaseline(cfg, control_ERP);

%% Separate conditions

fprintf('There are total of %d events\n',size(control_events,1))

corr_rspn_ERP_peakEnv = [];
no_rspn_ERP_peakEnv = [];
wrng_rspn_ERP_peakEnv = [];
corr_rspn_spect_peakEnv = [];
no_rspn_spect_peakEnv = [];
wrng_rspn_spect_peakEnv = [];

corr_rspn_ERP_peakRate = [];
no_rspn_ERP_peakRate = [];
wrng_rspn_ERP_peakRate = [];
corr_rspn_spect_peakRate = [];
no_rspn_spect_peakRate = [];
wrng_rspn_spect_peakRate = [];

wcorr_ind = 1;
wno_ind = 1;
wwrng_ind = 1;
scorr_ind = 1;
sno_ind = 1;
swrng_ind = 1;
fprintf('\nCurrent event:')


for t = 1:size(control_events,1)
    
    fprintf('-%d',t)
    % if only 1 correct answer, skip this trial
    if sum(strcmp(control_events.word_info{t}.response,'1')) < 2
        fprintf('Skipped')
        continue
    end
    
    % Separate peakEnv events
    % Collect fouri in a cell structure
    for pE = 1:length(control_events.peak_info{t}.peakEnv{1})
        
        % first check if the points are too close to each other
        if pE == length(control_events.peak_info{t}.peakEnv{1})
            % do nothing
        elseif control_events.peak_info{t}.peakEnv{1}(pE+1) - control_events.peak_info{t}.peakEnv{1}(pE) < 0.25
            continue
        end
        
        % check what is the response when that word was heard
        wword_resp = control_events.word_info{t}.response{...
            control_events.peak_info{t}.peakEnv{1}(pE)>=control_events.word_info{t}.onset & ...
            control_events.peak_info{t}.peakEnv{1}(pE)<control_events.word_info{t}.offset};
        
        % find closest timepoint for ERP
        cl_times = [nearest(control_ERP.time, control_events.peak_info{t}.peakEnv{1}(pE)-0.05)...
            nearest(control_ERP.time,control_events.peak_info{t}.peakEnv{1}(pE)+.5)];
        if strcmp(wword_resp,'1')
            corr_rspn_ERP_peakEnv(wcorr_ind,:,:) = squeeze(control_ERP.trial(t,:,cl_times(1):cl_times(2)));
        elseif strcmp(wword_resp,'0')
            no_rspn_ERP_peakEnv(wno_ind,:,:) = squeeze(control_ERP.trial(t,:,cl_times(1):cl_times(2)));
        else % wrong responses
            wrng_rspn_ERP_peakEnv(wwrng_ind,:,:) = squeeze(control_ERP.trial(t,:,cl_times(1):cl_times(2)));
        end
        
        % find closest timepoint for HFA and spectrogram
        cl_times = [nearest(control_wlt.time, control_events.peak_info{t}.peakEnv{1}(pE)-0.05)...
            nearest(control_wlt.time,control_events.peak_info{t}.peakEnv{1}(pE)+.5)];
        if strcmp(wword_resp,'1')
            corr_rspn_spect_peakEnv(wcorr_ind,:,:,:) = squeeze(control_wlt.powspctrm(t,:,:,cl_times(1):cl_times(2)));
            wcorr_ind = wcorr_ind+1;
        elseif strcmp(wword_resp,'0')
            no_rspn_spect_peakEnv(wno_ind,:,:,:) = squeeze(control_wlt.powspctrm(t,:,:,cl_times(1):cl_times(2)));
            wno_ind = wno_ind+1;
        else % wrong responses
            wrng_rspn_spect_peakEnv(wwrng_ind,:,:,:) = squeeze(control_wlt.powspctrm(t,:,:,cl_times(1):cl_times(2)));
            wwrng_ind = wwrng_ind+1;
        end
        
    end
    
    % Separate peakRate like different trials
    % Collect fouri in a cell structure
    for pR = 1:length(control_events.peak_info{t}.peakRate{1})
        
        % first check if the points are too close to each other
        if pR == length(control_events.peak_info{t}.peakRate{1})
            % do nothing
        elseif control_events.peak_info{t}.peakRate{1}(pR+1) - control_events.peak_info{t}.peakRate{1}(pR) < 0.25
            continue
        end
        
        % check what is the response when that word was heard
        wword_resp = control_events.word_info{t}.response{...
            control_events.peak_info{t}.peakRate{1}(pR)>=control_events.word_info{t}.onset & ...
            control_events.peak_info{t}.peakRate{1}(pR)<control_events.word_info{t}.offset};
        
        % find closest timepoint for ERP
        cl_times = [nearest(control_ERP.time, control_events.peak_info{t}.peakRate{1}(pR)-0.05)...
            nearest(control_ERP.time, control_events.peak_info{t}.peakRate{1}(pR)+.5)];
        if strcmp(wword_resp,'1')
            corr_rspn_ERP_peakRate(scorr_ind,:,:) = squeeze(control_ERP.trial(t,:,cl_times(1):cl_times(2)));
        elseif strcmp(wword_resp,'0')
            no_rspn_ERP_peakRate(sno_ind,:,:) = squeeze(control_ERP.trial(t,:,cl_times(1):cl_times(2)));
        else % wrong responses
            wrng_rspn_ERP_peakRate(swrng_ind,:,:) = squeeze(control_ERP.trial(t,:,cl_times(1):cl_times(2)));
        end
        
        
        % find closest timepoint for HFA and spectrogram
        cl_times = [nearest(control_wlt.time, control_events.peak_info{t}.peakRate{1}(pR)-0.05)...
            nearest(control_wlt.time, control_events.peak_info{t}.peakRate{1}(pR)+.5)];
        if strcmp(wword_resp,'1')
            corr_rspn_spect_peakRate(scorr_ind,:,:,:) = squeeze(control_wlt.powspctrm(t,:,:,cl_times(1):cl_times(2)));
            scorr_ind = scorr_ind+1;
        elseif strcmp(wword_resp,'0')
            no_rspn_spect_peakRate(sno_ind,:,:,:) = squeeze(control_wlt.powspctrm(t,:,:,cl_times(1):cl_times(2)));
            sno_ind = sno_ind+1;
        else % wrong responses
            wrng_rspn_spect_peakRate(swrng_ind,:,:,:) = squeeze(control_wlt.powspctrm(t,:,:,cl_times(1):cl_times(2)));
            swrng_ind = swrng_ind+1;
        end
    end
    
    
end
fprintf('\n')



%% Calculate ITPC
save_dir = fullfile(Sbj_Metadata.results, [strjoin(control_blocks,'_') '_v3']);
load(fullfile(save_dir, [strjoin(control_blocks,'_') '_ctrl_word_fouri.mat']),'fouri_of_words');
bwr = load('bwr_cmap.mat');

itpc = [];
for pp = 1:2 % first peakRate second peakEnv
    for cond = 1:2 % first correct second no response
        % get data
        if cond == 1 && pp == 2
            curr_fouri_all = fouri_of_words.corr_rspn_fouri_peakEnv{1};
        elseif cond == 2 && pp ==2
            curr_fouri_all = fouri_of_words.no_rspn_fouri_peakEnv{1};
        elseif cond == 1 && pp ==1
            curr_fouri_all = fouri_of_words.corr_rspn_fouri_peakRate{1};
        elseif cond == 2 && pp ==1
            curr_fouri_all = fouri_of_words.no_rspn_fouri_peakRate{1};
        end
        
        % compute inter-trial phase coherence (itpc) for each conditions
        tmp      = curr_fouri_all./abs(curr_fouri_all);    % divide by amplitude
        tmp      = sum(tmp,1);                            % sum angles across trials
        tmp      = abs(tmp)/size(curr_fouri_all,1);       % take the absolute value and normalize
        
        if isempty(tmp)
            sztmp = size(itpc);
            tmp = zeros(sztmp(3:5));
        end
        itpc(pp,cond,:,:,:) = squeeze(tmp);    % remove the first singleton dimension
        
    end
end

% for ITPC plot
load(fullfile(Sbj_Metadata.iEEG_data,Sbj_Metadata.BlockLists{1},[Sbj_Metadata.BlockLists{1} '_info.mat']),'info');
load(fullfile(Sbj_Metadata.sbjDir,[Sbj_Metadata.sbj_ID, '_channel_OI.mat']),'channel_OI')
freq_ITPC = fouri_of_words.freq_band_dtls{1};
time_ITPC = linspace(-0.05,0.5,size(fouri_of_words.corr_rspn_fouri_peakEnv{1},4));


% calculate stats
for i=1:2
    
%     if i==1
%         to_freq1 = fouri_of_words.corr_rspn_fouri_peakRate{1};
%         to_freq2 = fouri_of_words.no_rspn_fouri_peakRate{1};
%     else
%         to_freq1 = fouri_of_words.corr_rspn_fouri_peakEnv{1};
%         to_freq2 = fouri_of_words.no_rspn_fouri_peakEnv{1};
%     end
%     freq1 = [];
%     freq1.label = info.channelinfo.Label;
%     freq1.time = fouri_of_words.time_dtls;
%     freq1.freq = fouri_of_words.freq_band_dtls{1};
%     freq1.dimord = 'rpt_chan_freq_time';
%     freq2=freq1;
%     freq1.fourierspectrum = to_freq1;
%     freq2.fourierspectrum = to_freq2;
%     
%     cfg=[];
%     cfg.channel = channel_OI;
%     cfg.latency = [0 0.5];
%     cfg.statistic = 'ft_statfun_diff_itc';
% %     cfg.statistic       = 'indepsamplesT';
%     cfg.method='montecarlo';
%     cfg.correctm = 'bonferroni';
% %     cfg.correctm = 'cluster';
% %     cfg.clusterthreshold = 'nonparametric_individual';
%     cfg.parameter = 'fourierspectrum';
%     cfg.numrandomization = 1000;
%     cfg.design = [ones(1,size(to_freq1,1)),2*ones(1,size(to_freq2,1))];
%     [ITPC_stats(i)] = ft_freqstatistics(cfg, freq1,freq2);

    if i==1
        to_freq = cat(1,fouri_of_words.corr_rspn_fouri_peakRate{1},fouri_of_words.no_rspn_fouri_peakRate{1});
        size11 = size(fouri_of_words.corr_rspn_fouri_peakRate{1},1);
        size22 = size(fouri_of_words.no_rspn_fouri_peakRate{1},1);
    else
        to_freq = cat(1,fouri_of_words.corr_rspn_fouri_peakEnv{1},fouri_of_words.no_rspn_fouri_peakEnv{1});
        size11 = size(fouri_of_words.corr_rspn_fouri_peakEnv{1},1);
        size22 = size(fouri_of_words.no_rspn_fouri_peakEnv{1},1);
    end
    freq = [];
    freq.label = info.channelinfo.Label;
    freq.time = fouri_of_words.time_dtls;
    freq.freq = fouri_of_words.freq_band_dtls{1};
    freq.dimord = 'rpt_chan_freq_time';
    freq.fourierspectrum = to_freq;
    
    cfg=[];
    cfg.channel = channel_OI;
    cfg.latency = [0 0.5];
    cfg.statistic = 'ft_statfun_diff_itc';
%     cfg.statistic       = 'indepsamplesT';
    cfg.method='montecarlo';
    cfg.alpha       = 0.05;
    cfg.tail        = 0; % two-sided test
    cfg.correcttail = 'alpha';

    
    cfg.correctm = 'bonferroni';
%     cfg.correctm = 'cluster';
%     cfg.clusterthreshold = 'nonparametric_individual';
    cfg.parameter = 'fourierspectrum';
    cfg.numrandomization = 1000;
    cfg.design = [ones(1,size11),2*ones(1,size22)];
    [ITPC_stats(i)] = ft_freqstatistics(cfg, freq);


end

clear fouri_of_words pp cond tmp curr_fouri_all to_freq1 to_freq2

%%

% corr_rspn_ERP_peakEnv = [];
% no_rspn_ERP_peakEnv = [];
% wrng_rspn_ERP_peakEnv = [];
% corr_rspn_spect_peakEnv = [];
% no_rspn_spect_peakEnv = [];
% wrng_rspn_spect_peakEnv = [];
%
% corr_rspn_ERP_peakRate = [];
% no_rspn_ERP_peakRate = [];
% wrng_rspn_ERP_peakRate = [];
% corr_rspn_spect_peakRate = [];
% no_rspn_spect_peakRate = [];
% wrng_rspn_spect_peakRate = [];


time_ERP = linspace(-0.05,0.5,size(corr_rspn_ERP_peakEnv,3));
freq_spec = control_wlt.freq;
time_spec = linspace(-0.05,0.5,size(corr_rspn_spect_peakEnv,4));

cl_freqs=[];
for ii = 70:150
    cl_freqs(end+1) = nearest(control_wlt.freq, ii);
end
cl_freqs = unique(cl_freqs);

for el = 1:size(corr_rspn_ERP_peakRate,2)
    
    if ~any(ismember(channel_OI,control_ERP.label{el}))
        continue
    end
    
    figure('Units','normalized','Position', [0 0  1 1]);
    
    %% plot ERPs
    % plot single trials and average of corr - ERP - pR
    subplot(4,4,1)
    for_avg=squeeze(corr_rspn_ERP_peakRate(:,el,:));
    plot(time_ERP,for_avg,'Color',[0.9412,0.9412,0.9412])
    hold on
    shadedErrorBar(time_ERP,mean(for_avg,1),stderr(for_avg),'lineprops','r')
    title([num2str(size(corr_rspn_ERP_peakRate,1)) ' correct response - peakRate locked - ERP'])
    xlim([time_ERP(1) time_ERP(end)])
    set(gca, 'FontSize',13,'FontWeight','bold');
    plot([0 0], ylim,'k')
    
    % plot single trials and average of no - ERP - pR
    subplot(4,4,5)
    for_avg=squeeze(no_rspn_ERP_peakRate(:,el,:));
    plot(time_ERP,for_avg,'Color',[0.9412,0.9412,0.9412])
    hold on
    shadedErrorBar(time_ERP,mean(for_avg,1),stderr(for_avg),'lineprops','r')
    title([num2str(size(no_rspn_ERP_peakRate,1)) ' no response - peakRate locked - ERP'])
    xlim([time_ERP(1) time_ERP(end)])
    set(gca, 'FontSize',13,'FontWeight','bold');
    plot([0 0], ylim,'k')
    
    % plot single trials and average of corr - ERP - pE
    subplot(4,4,9)
    for_avg=squeeze(corr_rspn_ERP_peakEnv(:,el,:));
    plot(time_ERP,for_avg,'Color',[0.9412,0.9412,0.9412])
    hold on
    shadedErrorBar(time_ERP,mean(for_avg,1),stderr(for_avg),'lineprops','r')
    title([num2str(size(corr_rspn_ERP_peakEnv,1)) ' correct response - peakEnv locked - ERP'])
    xlim([time_ERP(1) time_ERP(end)])
    set(gca, 'FontSize',13,'FontWeight','bold');
    plot([0 0], ylim,'k')
    
    % plot single trials and average of no - ERP - pE
    subplot(4,4,13)
    for_avg=squeeze(no_rspn_ERP_peakEnv(:,el,:));
    plot(time_ERP,for_avg,'Color',[0.9412,0.9412,0.9412])
    hold on
    shadedErrorBar(time_ERP,mean(for_avg,1),stderr(for_avg),'lineprops','r')
    title([num2str(size(no_rspn_ERP_peakEnv,1)) ' no response - peakEnv locked - ERP'])
    xlim([time_ERP(1) time_ERP(end)])
    set(gca, 'FontSize',13,'FontWeight','bold');
    plot([0 0], ylim,'k')
    
    
    %% plot spectrograms
    % spectrogram of corr - pR
    subplot(4,4,2)
    avg_spec = squeeze(mean(corr_rspn_spect_peakRate(:,el,:,:),1));
    h = pcolor(time_spec,freq_spec,avg_spec);
    h.EdgeColor = 'none';
%     imagesc(time_spec,freq_spec,avg_spec)
    axis xy
    set(gca,'YScale','log')
    set(gca, 'FontSize',13,'FontWeight','bold');
    caxis([-7 7]);hold on
    plot([0 0], ylim,'k')
    title('Spectrograms')
    % spectrogram of no - pR
    subplot(4,4,6)
    avg_spec = squeeze(mean(no_rspn_spect_peakRate(:,el,:,:),1));
%     imagesc(time_spec,freq_spec,avg_spec)
    h = pcolor(time_spec,freq_spec,avg_spec);
    h.EdgeColor = 'none';
    axis xy
    set(gca,'YScale','log')
    set(gca, 'FontSize',13,'FontWeight','bold');
    caxis([-7 7]);hold on
    plot([0 0], ylim,'k')
    % spectrogram of corr - pE
    subplot(4,4,10)
    avg_spec = squeeze(mean(corr_rspn_spect_peakEnv(:,el,:,:),1));
%     imagesc(time_spec,freq_spec,avg_spec)
    h = pcolor(time_spec,freq_spec,avg_spec);
    h.EdgeColor = 'none';
    set(gca,'YScale','log')
    set(gca, 'FontSize',13,'FontWeight','bold');
    axis xy;hold on
    caxis([-7 7])
    plot([0 0], ylim,'k')
    % spectrogram of no - pE
    subplot(4,4,14)
    avg_spec = squeeze(mean(no_rspn_spect_peakEnv(:,el,:,:),1));
%     imagesc(time_spec,freq_spec,avg_spec)
    h = pcolor(time_spec,freq_spec,avg_spec);
    h.EdgeColor = 'none';
    set(gca,'YScale','log')
    set(gca, 'FontSize',13,'FontWeight','bold');
    axis xy;hold on
    caxis([-7 7])
    plot([0 0], ylim,'k')
    
    
    
    %% plot HFA
    % cl_freqs: HFA frequency range
    % HFA of corr - pR
    subplot(4,4,3)
    avg_HFA = squeeze(mean(corr_rspn_spect_peakRate(:,el,cl_freqs,:),1));
    plot(time_spec,avg_HFA,'Color',[0.9412,0.9412,0.9412])
    hold on
    shadedErrorBar(time_spec,mean(avg_HFA,1),stderr(avg_HFA),'lineprops','r')
    title([num2str(size(corr_rspn_spect_peakRate,1)) ' correct response - peakRate locked - HFA'])
    xlim([time_spec(1) time_spec(end)])
    ylim([-8 8])
    plot(xlim,[0 0],'k')
    set(gca, 'FontSize',13,'FontWeight','bold');
    ylabel('dB')
    plot([0 0], ylim,'k')
    % HFA of no - pR
    subplot(4,4,7)
    avg_HFA = squeeze(mean(no_rspn_spect_peakRate(:,el,cl_freqs,:),1));
    plot(time_spec,avg_HFA,'Color',[0.9412,0.9412,0.9412])
    hold on
    shadedErrorBar(time_spec,mean(avg_HFA,1),stderr(avg_HFA),'lineprops','r')
    title([num2str(size(no_rspn_spect_peakRate,1)) ' no response - peakRate locked - HFA'])
    xlim([time_spec(1) time_spec(end)])
    ylim([-8 8])
    plot(xlim,[0 0],'k')
    set(gca, 'FontSize',13,'FontWeight','bold');
    ylabel('dB')
    plot([0 0], ylim,'k')
    % HFA of corr - pE
    subplot(4,4,11)
    avg_HFA = squeeze(mean(corr_rspn_spect_peakEnv(:,el,cl_freqs,:),1));
    plot(time_spec,avg_HFA,'Color',[0.9412,0.9412,0.9412])
    hold on
    shadedErrorBar(time_spec,mean(avg_HFA,1),stderr(avg_HFA),'lineprops','r')
    title([num2str(size(corr_rspn_spect_peakEnv,1)) ' correct response - peakEnv locked - HFA'])
    xlim([time_spec(1) time_spec(end)])
    ylim([-8 8])
    plot(xlim,[0 0],'k')
    set(gca, 'FontSize',13,'FontWeight','bold');
    ylabel('dB')
    plot([0 0], ylim,'k')
    % HFA of no - pE
    subplot(4,4,15)
    avg_HFA = squeeze(mean(no_rspn_spect_peakEnv(:,el,cl_freqs,:),1));
    plot(time_spec,avg_HFA,'Color',[0.9412,0.9412,0.9412])
    hold on
    shadedErrorBar(time_spec,mean(avg_HFA,1),stderr(avg_HFA),'lineprops','r')
    title([num2str(size(no_rspn_spect_peakEnv,1)) ' no response - peakEnv locked - HFA'])
    xlim([time_spec(1) time_spec(end)])
    ylim([-8 8])
    plot(xlim,[0 0],'k')
    set(gca, 'FontSize',13,'FontWeight','bold');
    ylabel('dB')
    plot([0 0], ylim,'k')
    
    
    %% plot ITPC
    for pp = 1:2 % first peakRate second peakEnv
        % cond = 1:2 % first correct second no response
        
        subplot(4, 4, [subplotno(4,2*pp+1-2,4) subplotno(4,2*pp+2-2,4)])
        
        imagesc(time_ITPC, freq_ITPC, [squeeze(itpc(pp,1,el,:,:))-squeeze(itpc(pp,2,el,:,:))]);
        axis xy
        set(gca, 'FontSize',13,'FontWeight','bold');
        caxis([-0.5 0.5])
        
        if pp == 1
            ylabel({'peakRate locked';'Frequency (Hz)'});
        elseif pp == 2
            ylabel({'peakEnv locked';'Frequency (Hz)'});
        end
        title('ITPC [Correct-No] responses');
        
        hold on
        signplot = [zeros(size(ITPC_stats(pp).mask,2),find(time_ITPC==0)-1),double(squeeze(ITPC_stats(pp).mask(strcmp(ITPC_stats(pp).label,control_ERP.label{el}),:,:)))];
        contour(time_ITPC, freq_ITPC,signplot,1,'LineColor','k','LineWidth',3)
        plot([0 0], ylim,'k')
        
        
        
    end
    
    % Create and delete new axes to plot colorbar of spectrogram
    ax = axes;
    colormap(bwr.rgb_vals);
    cmaph = colorbar(ax);
    cmaph.Ticks = linspace(0,1,5);
    cmaph.TickLabels = num2cell(linspace(-7,7,5));
    cmaph.FontSize = 13;cmaph.FontWeight='bold';
    cmaph.LineWidth = 1;
    colorTitleHandle = get(cmaph,'Title');
    set(colorTitleHandle ,'String','dB-power','FontSize',13,'FontWeight','bold','Position',[149.5800 -30 0]);
    
    a=get(cmaph); %gets properties of colorbar
    a = a.Position; %gets the positon and size of the color bar
    set(cmaph,'Location','southoutside') % to change orientation
    set(cmaph,'Position',[a(1)/4+0.11 0.04 0.16 0.02]) % To change size    
    ax.Visible = 'off';
    
    % Create and delete new axes to plot colorbar of ITPC
    cmaph2 = colorbar(ax);
    cmaph2.Ticks = linspace(0,1,6);
    cmaph2.TickLabels = num2cell(linspace(-0.5,0.5,6));
    cmaph2.FontSize = 13;cmaph2.FontWeight='bold';
    cmaph2.LineWidth = 1;
    colorTitleHandle = get(cmaph2,'Title');
    set(colorTitleHandle ,'String','ITPC','FontSize',13,'FontWeight','bold','Position',[149.5800 -30 0]);
    
    a2=get(cmaph2); %gets properties of colorbar
    a2 = a2.Position; %gets the positon and size of the color bar
    set(cmaph2,'Location','southoutside') % to change orientation
    set(cmaph2,'Position',[3*a(1)/4+0.08 0.04 0.16 0.02]) % To change size    
    ax.Visible = 'off';
    text(-0.07,0.6,'peakRate locked events','Units','normalized','Rotation',90,'FontSize',18,'FontWeight','bold')
    text(-0.07,0.1,'peakEnv locked events','Units','normalized','Rotation',90,'FontSize',18,'FontWeight','bold')
    
    sgtitle(['Elec: ' info.channelinfo.Label{el} ' - ERP-spectrogram-HFA-ITPC'], 'FontSize',15,'FontWeight','bold')
    
    % Save the figure
    fprintf('\t-Saving electrode #%d-%s, out of %d\n',el,control_ERP.label{el},size(control_ERP.label,1))
    print(fullfile(save_folder,[control_wlt.label{el} , '_pEvents.jpg']),'-djpeg','-r300')
    close all
    
end


end
