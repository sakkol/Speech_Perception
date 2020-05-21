function SP_control_spectrogram(Sbj_Metadata,control_blocks,runagain,onlySuccess)
% This function streamlines: loading wavelet output from each block,
% combines only control trials (if not done or if need to be redone),
% baseline correct time-freq map, then plots spectrogram, HFA and ERP in
% 1x3 figure.

%% Select blocks to import
vars=who;
if ~ismember(vars,'control_blocks')
    control_blocks = select_cont_blocks(Sbj_Metadata);
elseif isempty(control_blocks)
    control_blocks = select_cont_blocks(Sbj_Metadata);
end
if ~ismember(vars,'runagain')
    runagain = 0;
end
if ~ismember(vars,'onlySuccess')
    onlySuccess = 0;
    save_folder = fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),'control_spectrogram');
elseif onlySuccess == 1
    save_folder = fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),'control_spectrogram_onlySuccess');
elseif onlySuccess == 0
    save_folder = fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),'control_spectrogram');
end
if ~exist(save_folder,'dir'),mkdir(save_folder),end
clear vars

%% bring in these blocks and combine only the control events
fprintf('These blocks are going to be used: %s\n',strjoin(control_blocks,', '))
% load info.mat file
load(fullfile(Sbj_Metadata.iEEG_data,Sbj_Metadata.BlockLists{1},[Sbj_Metadata.BlockLists{1} '_info.mat']),'info');

% check if this has already been run or need to be re-run
if ~exist(fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),[strjoin(control_blocks,'_') '_control_wltERP.mat']),'file') || runagain
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

if onlySuccess
    % select the trials that the accuracy is higher than 3 words
    succ_idx = control_events.Acc_word_count > 3;
    fprintf('Only %d successful trials are to be used.\n',sum(succ_idx))
    cfg = [];
    cfg.trials = succ_idx;
    [control_wlt] = ft_selectdata(cfg, control_wlt);
    control_ERP = ft_selectdata(cfg, control_ERP);
end

%% Baseline correct TF data but not ERPs
%  Baseline correct time-freq data
baseline = [-3.4 -3.1];
cfg              = [];
cfg.baseline     =  baseline;% seconds (prespeech part is 3.5455 seconds) [0.5sec + 3.0455sec]
cfg.baselinetype = 'db';
cfg.parameter    = 'powspctrm';
[control_wlt]    = ft_freqbaseline(cfg, control_wlt);

cfg               = [];
cfg.frequency     = [70 150];
cfg.avgoverfreq   = 'yes';
freq_sel = ft_selectdata(cfg, control_wlt);

% cfg              = [];
% cfg.baseline     = [-3.45 -3.05];
% cfg.baselinetype = 'db';
% cfg.channel      = 'all';
% cfg.parameter    = 'trial';
% control_ERP      = ft_timelockbaseline(cfg, control_ERP);

%% Prepare for plot
% load channels of interest
load(fullfile(Sbj_Metadata.sbjDir,[Sbj_Metadata.sbj_ID, '_channel_OI.mat']),'channel_OI')

freq_spec = control_wlt.freq;
time_spec = control_wlt.time;

time_ERP = control_ERP.time;
% cl_freqs=[];
% for ii = 70:150
%     cl_freqs(end+1) = nearest(control_wlt.freq, ii);
% end
% cl_freqs = unique(cl_freqs);
bwr = load('bwr_cmap.mat');

%% Plot

for el = 1:size(control_ERP.trial,2)
    
    if ~any(ismember(channel_OI,control_ERP.label{el}))
        continue
    end
    
    figure('Units','normalized','Position', [0 0  .4 1]);
    
    % plot ERPs
    subplot(3,1,1)
    for_avg=squeeze(control_ERP.trial(:,el,:));
    plot(time_ERP,for_avg,'Color',[0.9412,0.9412,0.9412])
    hold on
    shadedErrorBar(time_ERP,mean(for_avg,1),stderr(for_avg),'lineprops','r')
    title(['Target sentence onset locked ERP'])
    set(gca, 'FontSize',13,'FontWeight','bold');
    xlim([time_spec(1) time_spec(end)])
    hold on
    ylims = ylim;
    plot([-3.54 -3.54], ylim,'k') % trial onset
    text(-3.54,ylims(2)-(ylims(2)-ylims(1))/10,'Trial onset')
    
    plot([-3.045 -3.045], ylim,'k') % 3.045
    text(-3.045,ylims(2)-2*(ylims(2)-ylims(1))/10,'Attention sent. onset')
    
    plot([0 0], ylim,'k')
    text(0,ylims(2)-(ylims(2)-ylims(1))/10,'Target sent. onset')
    
    plot([2.2 2.2], ylim,'k') % sentence offset
    text(2.2,ylims(2)-(ylims(2)-ylims(1))/10,'Target sent. ends')
    
    plot([4.2 4.2], ylim,'k') % trial ends
    text(4.2,ylims(2)-(ylims(2)-ylims(1))/10,'Trial ends')
    ylabel('Voltage (uV)')
    ylim(ylims)
    
    % spectrogram
    subplot(3,1,2)
    avg_spec = squeeze(nanmean(control_wlt.powspctrm(:,el,:,:),1));
    h = pcolor(time_spec,freq_spec,avg_spec);
    h.EdgeColor = 'none';
    %     imagesc(time_spec,freq_spec,avg_spec)
    axis xy
    set(gca,'YScale','log')
    yticks([1,4,8,12,20,50,70,100,150,200])
    set(gca, 'FontSize',13,'FontWeight','bold');
    caxis([-7 7]);
    hold on
    ylims = ylim;
    plot([-3.54 -3.54], ylim,'k') % trial onset
    text(-3.54,ylims(2)-(ylims(2)-ylims(1))/10,'Trial onset')
    
    plot([-3.045 -3.045], ylim,'k') % 3.045
    text(-3.045,ylims(2)-3*(ylims(2)-ylims(1))/10,'Attention sent. onset')
    
    plot([0 0], ylim,'k')
    text(0,ylims(2)-(ylims(2)-ylims(1))/10,'Target sent. onset')
    
    plot([2.2 2.2], ylim,'k') % sentence offset
    text(2.2,ylims(2)-(ylims(2)-ylims(1))/10,'Target sent. ends')
    
    plot([4.2 4.2], ylim,'k') % trial ends
    text(4.2,ylims(2)-(ylims(2)-ylims(1))/10,'Trial ends')
    xlim([time_spec(1) time_spec(end)])
    ylabel('log frequency (Hz)')
    colorbar
    
    % plot HFA
    % cl_freqs: HFA frequency range
    subplot(3,1,3)
    avg_HFA = squeeze(freq_sel.powspctrm(:,el,:,:));
    plot(time_spec,avg_HFA,'Color',[0.9412,0.9412,0.9412])
    hold on
    shadedErrorBar(time_spec,mean(avg_HFA,1),stderr(avg_HFA),'lineprops','r')
    title(['Target sentence onset locked HFA'])
    xlim([time_spec(1) time_spec(end)])
    ylim([-8 8])
    plot(xlim,[0 0],'k')
    set(gca, 'FontSize',13,'FontWeight','bold');
    ylabel('HFA (dB)')
    ylims = ylim;
    plot([-3.54 -3.54], ylim,'k') % trial onset
    text(-3.54,ylims(2)-(ylims(2)-ylims(1))/10,'Trial onset')
    
    plot([-3.045 -3.045], ylim,'k') % 3.045
    text(-3.045,ylims(2)-2*(ylims(2)-ylims(1))/10,'Attention sent. onset')
    
    plot([0 0], ylim,'k')
    text(0,ylims(2)-(ylims(2)-ylims(1))/10,'Target sent. onset')
    
    plot([2.2 2.2], ylim,'k') % sentence offset
    text(2.2,ylims(2)-(ylims(2)-ylims(1))/10,'Target sent. ends')
    
    plot([4.2 4.2], ylim,'k') % trial ends
    text(4.2,ylims(2)-(ylims(2)-ylims(1))/10,'Trial ends')
    xlim([time_spec(1) time_spec(end)])
    
    colormap(bwr.rgb_vals);
    
    sgtitle(['Elec: ' info.channelinfo.Label{el} ' - ERP-spectrogram-HFA; ' num2str(size(control_wlt.powspctrm,1)) ' trials - baseline correction in '...
        num2str(baseline(1)) ' - ' num2str(baseline(2))], 'FontSize',15,'FontWeight','bold')
    
    % Save the figure
    fprintf('\t-Saving electrode #%d-%s, out of %d\n',el,control_ERP.label{el},size(channel_OI,1))
    print(fullfile(save_folder,[control_wlt.label{el} , '_spectrograms.jpg']),'-djpeg','-r300')
    close all
end