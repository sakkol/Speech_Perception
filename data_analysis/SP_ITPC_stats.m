

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