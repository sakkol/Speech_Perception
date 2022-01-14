function [epoched_data, epoched_wlt] = IL_mtmconvolTF(Sbj_Metadata,curr_block,ftrip,event_secs,pre,post,elec,foi,fourier_pow)
% This function performs the time-frequency decomposition so that there is
% no clutters in scripts. Also, during the process the epoched_data is
% produced. ftrip is the fieldtrip structure of continuous data (probably
% common average referenced); events_sec is the onsets of where to epoch
% data (=time zero); pre and post are positive values of epoching seconds.
% Example code:
% pre  = 1;
% post = 4;
% [epoched_data, epoched_wlt] = master_TF(ecog_avg.ftrip,events.Click1,pre,post);

%% Set the stage
if ~exist('event_secs','var') || isempty(event_secs)
    error('There is no events info.')
end
if ~exist('pre','var') || isempty(pre)
    pre = 2;
end
if ~exist('post','var') || isempty(post)
    post = 4;
end
if ~exist('elec','var') || isempty(elec)
    error('Provide electrode name')
end
if ~exist('foi','var') || isempty(foi)
    foi = [0.2:0.1:6,50:7:200];
end
if ~exist('fourier_pow','var') || isempty(fourier_pow)
    fourier_pow = 'fourier';
end

%% Preprocessing
% select channel
cfg           = [];
cfg.channel   = elec;
ftrip         = ft_selectdata(cfg,ftrip);
% First resample for saving space and memory, if necessary
if ftrip.fsample < 998 || ftrip.fsample > 1002
    cfg             = [];
    cfg.resamplefs  = 1000;
    ftrip  = ft_resampledata(cfg,ftrip);
end
fs = ftrip.fsample; % sampling rate

% Make trial structure
trl           = [];
trl(:,1)      = floor( event_secs*fs - fs*(pre+6) ); % adding 2 second for lower frequencies, will cut it out later
trl(:,2)      = floor( event_secs*fs + fs*(post+6) );
trl(:,3)      = floor( -(pre+6)*fs );

% Epoch
cfg           = [];
cfg.trl       = trl;
epoched_data  = ft_redefinetrial(cfg,ftrip);

% compute trials
cfg             = [];
cfg.keeptrials  = 'yes';
epoched_data     = ft_timelockanalysis(cfg,epoched_data);

% check if there is NaNs due to very early start of trial in the beginning mostly
if sum(sum(sum(isnan(epoched_data.trial))))
    warning('\n\tTHERE WERE %s TIME-POINTS THAT HAVE NANs, REPLACING WITH ZEROS\n',num2str(sum(sum(sum(isnan(epoched_data.trial))))))
    % replace NaNs with zeros
    epoched_data.trial(isnan(epoched_data.trial)) = 0;
end

%% Wavelet
cfg                   = [];
cfg.method            = 'mtmconvol';
cfg.taper             = 'hanning';
cfg.foi               = foi;
cfg.toi               = -pre:0.01:post;
cfg.t_ftimwin         = [6*ones(sum(cfg.foi<10),1);.1*ones(sum(cfg.foi>10),1)];
% cfg.t_ftimwin         = 16./cfg.foi;
cfg.keeptrials        = 'yes';
cfg.keeptaper         = 'yes';
cfg.output            = fourier_pow;
cfg.pad               = 'maxperlen';
epoched_wlt           = ft_freqanalysis(cfg,epoched_data);

% Shorten the epoched_data
cfg                = [];
cfg.latency      = [-pre post];
epoched_data     = ft_selectdata(cfg, epoched_data);

% save
save_folder = fullfile(Sbj_Metadata.iEEG_data, curr_block, 'elec_TF');
if ~exist(save_folder,'dir'),mkdir(save_folder),end
    
save(fullfile(save_folder, [elec '_' curr_block '_mtmconvol_fourier.mat']),'epoched_wlt','epoched_data','-v7.3');

end
