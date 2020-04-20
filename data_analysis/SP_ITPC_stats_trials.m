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