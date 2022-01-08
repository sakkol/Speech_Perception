%% simulate something
fsample = 1000;
nsample = 100*fsample;

% start with a continuous representation of 100 seconds of data
data.label = {'eeg'};
data.time = {((1:nsample)/fsample)};
data.trial{1,1} = .1*randn(size(data.time{1}));

% add a trigger every 1/2.4s
i = 1/2.4;
smalljump = 4*i:i:data.time{1}(end-1500);
midjumps = smalljump(1:2:end);
bigjumps = midjumps(1:2:end);
trials = bigjumps(1:5:end);trials(end)=[];
jumpp = .5+(-.5*linspace(-1, 1, 100).^2);

for s=1:length(smalljump)-1
    data.trial{1}(nearest(data.time{1},smalljump(s)):nearest(data.time{1},smalljump(s))+99) = ...
        data.trial{1}(nearest(data.time{1},smalljump(s)):nearest(data.time{1},smalljump(s))+99) + jumpp;
end
for s=1:length(midjumps)
    data.trial{1}(nearest(data.time{1},midjumps(s)):nearest(data.time{1},midjumps(s))+99) = ...
        data.trial{1}(nearest(data.time{1},midjumps(s)):nearest(data.time{1},midjumps(s))+99) + jumpp;
end
for s=1:length(bigjumps)
    data.trial{1}(nearest(data.time{1},bigjumps(s)):nearest(data.time{1},bigjumps(s))+99) = ...
        data.trial{1}(nearest(data.time{1},bigjumps(s)):nearest(data.time{1},bigjumps(s))+99) + jumpp;
end

% Make trial structure
pre=0;post=8.33;
trl           = [];
trl(:,1)      = floor( trials*fsample - fsample*(pre+.1) ); % adding 2 second for lower frequencies, will cut it out later
trl(:,2)      = floor( trials*fsample + fsample*(post+.1) );
trl(:,3)      = floor( -(pre+.1)*fsample );

% Epoch
cfg           = [];
cfg.trl       = trl;
data_epo  = ft_redefinetrial(cfg,data);
% compute trials
cfg             = [];
cfg.keeptrials  = 'yes';
data_epo     = ft_timelockanalysis(cfg,data_epo);

% plot briefly
figure;plot(data_epo.time,squeeze(data_epo.trial(1,:,:)))
xlim([data_epo.time(1) data_epo.time(end)])

%% run and plot the TFs
figure('Units','normalized','Position', [0 0  1 1]);

% Subplot-1 : my original mtmconvol-hanning code
cfg                   = [];
cfg.method            = 'mtmconvol';
cfg.taper             = 'hanning';
cfg.foi               = [0.2:0.1:5,50:5:200];
cfg.toi               = -pre:0.01:post;
cfg.t_ftimwin         = [5*ones(sum(cfg.foi<10),1);.1*ones(sum(cfg.foi>10),1)]; % for low freqs, time-window is 5 sec, for high freq, it is 0.1 sec
cfg.keeptrials        = 'yes';
cfg.keeptaper         = 'yes';
cfg.output            = 'pow';
cfg.pad               = 'maxperlen';
wlt1           = ft_freqanalysis(cfg,data_epo);
subplot(241)
plot(wlt1.freq,squeeze(nanmean(wlt1.powspctrm,4)))
xlim([0 3])

% Subplot-2 : mtmconvol-hanning trial
cfg                   = [];
cfg.method            = 'mtmconvol';
cfg.taper             = 'hanning';
cfg.foi               = [0.2:0.1:5,50:5:200];
cfg.toi               = -pre:0.01:post;
cfg.t_ftimwin         = 3./cfg.foi; % different time windows for each freq: 3 /frequency
cfg.keeptrials        = 'yes';
cfg.keeptaper         = 'yes';
cfg.output            = 'pow';
cfg.pad               = 'maxperlen';
wlt2           = ft_freqanalysis(cfg,data_epo);
subplot(242)
plot(wlt2.freq,squeeze(nanmean(wlt2.powspctrm,4)))
xlim([0 3])

% Subplot-3 : mtmconvol-dpss trial1
cfg                   = [];
cfg.method            = 'mtmconvol';
cfg.taper             = 'dpss';
cfg.tapsmofrq         = 0.15;
cfg.foi               = [0.2:0.1:5];
cfg.toi               = -pre:0.01:post;
cfg.t_ftimwin         = [6*ones(sum(cfg.foi<10),1)]; % for low freqs, time-window is 6 sec
cfg.keeptrials        = 'yes';
cfg.keeptaper         = 'yes';
cfg.output            = 'pow';
cfg.pad               = 'maxperlen';
wlt3           = ft_freqanalysis(cfg,data_epo);
subplot(243)
plot(wlt3.freq,squeeze(nanmean(wlt3.powspctrm,4)))
xlim([0 3])

% Subplot-4 : mtmconvol-dpss trial2
cfg                   = [];
cfg.method            = 'mtmconvol';
cfg.taper             = 'dpss';
cfg.tapsmofrq         = 0.2;
cfg.foi               = [0.2:0.1:5];
cfg.toi               = -pre:0.01:post;
cfg.t_ftimwin         = [20./cfg.foi(cfg.foi<10)]; % different time windows for each freq: 20 /frequency
cfg.keeptrials        = 'yes';
cfg.keeptaper         = 'yes';
cfg.output            = 'pow';
cfg.pad               = 'maxperlen';
wlt4           = ft_freqanalysis(cfg,data_epo);
subplot(244)
plot(wlt4.freq,squeeze(nanmean(wlt4.powspctrm,4)))
xlim([0 3])


% Subplot-5 : my original wavelet params
cfg                   = [];
cfg.method            = 'wavelet';
cfg.foi               = 0.2:0.1:5;
cfg.width             = 4*ones(sum(cfg.foi<10),1);
cfg.toi               = -pre:0.01:post;
cfg.keeptrials        = 'yes';
cfg.keeptaper         = 'yes';
cfg.output            = 'pow';
cfg.pad               = 'nextpow2';
wlt5           = ft_freqanalysis(cfg,data_epo);
subplot(245)
plot(wlt5.freq,squeeze(nanmean(wlt5.powspctrm,4)))
xlim([0 3])

% Subplot-5 : wavelet trial
cfg                   = [];
cfg.method            = 'wavelet';
cfg.foi               = 0.2:0.1:5;
cfg.width             = 16*ones(sum(cfg.foi<10),1);
cfg.toi               = -pre:0.01:post;
cfg.keeptrials        = 'yes';
cfg.keeptaper         = 'yes';
cfg.output            = 'pow';
cfg.pad               = 'nextpow2';
wlt6           = ft_freqanalysis(cfg,data_epo);
subplot(246)
plot(wlt6.freq,squeeze(nanmean(wlt6.powspctrm,4)))
xlim([0 3])

% Subplot-7 : superlet trial
cfg                   = [];
cfg.method            = 'superlet';
cfg.foi               = 0.2:0.1:5;
cfg.toi               = -pre:0.01:post;
cfg.superlet          = [];
cfg.keeptrials        = 'yes';
cfg.keeptaper         = 'yes';
cfg.output            = 'pow';
wlt7           = ft_freqanalysis(cfg,data_epo);
subplot(247)
plot(wlt7.freq,squeeze(nanmean(wlt7.powspctrm,4)))
xlim([0 3])





