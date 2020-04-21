%% First attempt

freq1 = [];
freq1.label = info.channelinfo.Label;
freq1.time = fouri_of_words.time_dtls;
freq1.freq = fouri_of_words.freq_band_dtls{1};
freq1.dimord = 'rpt_chan_freq_time';
freq2=freq1;
freq1.fourierspectrum = fouri_of_words.corr_rspn_fouri_peakEnv{1};
freq2.fourierspectrum = fouri_of_words.no_rspn_fouri_peakEnv{1};

cfg=[];
cfg.channel = ft_channelselection('RTs*', freq1);
cfg.latency = [0.05 0.5];
cfg.statistic = 'ft_statfun_diff_itc';
cfg.statistic       = 'indepsamplesT';
cfg.method='montecarlo';
cfg.correctm = 'bonferroni';
cfg.parameter = 'fourierspectrum';
cfg.numrandomization = 1000;
cfg.design = [ones(1,size(fouri_of_words.corr_rspn_fouri_peakEnv{1},1)),2*ones(1,size(fouri_of_words.no_rspn_fouri_peakEnv{1},1))];
[stat] = ft_freqstatistics(cfg, freq1,freq2);

%% Second attempt

i = [];
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


%% quick peak at which have significance
for el=113:128
    for pp=1:2
        fprintf('El:%d, pp:%d - any ITPC:%d\n',el,pp,sum(sum(double(squeeze(ITPC_stats(pp).mask(strcmp(ITPC_stats(pp).label,control_ERP.label{el}),:,:))))))
        
    end
end

%% Third attempt

% Select blocks to import
control_blocks = select_cont_blocks(Sbj_Metadata);
save_dir = fullfile(Sbj_Metadata.results, [strjoin(control_blocks,'_') '_v3']);
fprintf('Loading ''fouri_of_words'' from:\n-->%s\n',fullfile(save_dir, [strjoin(control_blocks,'_') '_ctrl_word_fouri.mat']))
load(fullfile(save_dir, [strjoin(control_blocks,'_') '_ctrl_word_fouri.mat']));
load(fullfile(Sbj_Metadata.iEEG_data,Sbj_Metadata.BlockLists{1},[Sbj_Metadata.BlockLists{1} '_info.mat']),'info');

freq = fouri_of_words.freq_band_dtls{1};
time = linspace(-0.05,0.5,size(fouri_of_words.corr_rspn_fouri_peakEnv{1},4));

itpc = [];
for pp = 1:2
    for cond = 1:3
        % get data
        if cond == 1 && pp ==1
            curr_fouri_all = fouri_of_words.corr_rspn_fouri_peakEnv{1};
        elseif cond == 2 && pp ==1
            curr_fouri_all = fouri_of_words.no_rspn_fouri_peakEnv{1};
        elseif cond == 3 && pp ==1
            curr_fouri_all = fouri_of_words.wrng_rspn_fouri_peakEnv{1};
        elseif cond == 1 && pp ==2
            curr_fouri_all = fouri_of_words.corr_rspn_fouri_peakRate{1};
        elseif cond == 2 && pp ==2
            curr_fouri_all = fouri_of_words.no_rspn_fouri_peakRate{1};
        elseif cond == 3 && pp ==2
            curr_fouri_all = fouri_of_words.wrng_rspn_fouri_peakRate{1};
        end
        
        % compute inter-trial phase coherence (itpc) for each conditions
        tmp      = curr_fouri_all./abs(curr_fouri_all);    % divide by amplitude
        tmp      = sum(tmp,1);                            % sum angles across trials
        tmp      = abs(tmp)/size(curr_fouri_all,1);       % take the absolute value and normalize
        
        if isempty(tmp)
            sztmp = size(itpc);
            tmp = zeros(sztmp(3:5));
        end
        itpc(pp,cond,:,:,:) = squeeze(tmp);                % remove the first singleton dimension
        
        % z-score ITPC??????? or n*ITPC^2????
        
        
    end
end



i = [];
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
