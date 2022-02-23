function [epoched_hfa,epoched_data,select_events] = SP_getTargetHFA(Sbj_Metadata, curr_block)
% just to grab the HFA file easily
hfa_file = fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_TargetHFA.mat']);

if ~exist(hfa_file,'file')
    fprintf('Couldn''t find the HFA file, running the analysis and saving to: %s\n', hfa_file)
    tmp=load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_info.mat']));info=tmp.info;clear tmp
    load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_avg.mat']))
    
    % only the control HFA of target sentences
    select_events = info.events(info.events.Cond_code==1,:);
    speech_onsets = info.events.speech_onsets(info.events.Cond_code==1);
    
    % run wavelet
    pre  = 4.55;
    post = 9;
    freq = 70:5:200;freq(ismember(freq,[120 180]))=[];
    fourier_pow = 'pow';
    [epoched_data, epoched_wlt] = master_waveletTF(ecog_avg.ftrip,speech_onsets,pre,post,freq,fourier_pow);
    
    % get and baseline correct HFA
    cfg=[];
    cfg.frequency=[69 201];
    cfg.avgoverfreq='no';
    epoched_hfa = ft_selectdata(cfg, epoched_wlt);
    cfg              = [];
    cfg.baseline     = [-4.5 -4.1];
    cfg.baselinetype = 'relative';
    cfg.parameter    = 'powspctrm';
    [epoched_hfa]    = ft_freqbaseline(cfg, epoched_hfa);
    cfg=[];
    cfg.frequency=[70 200];
    cfg.avgoverfreq='yes';
    cfg.latency = [-1 9];  % get to save space
    epoched_hfa = ft_selectdata(cfg, epoched_hfa);
    cfg=[];
    cfg.latency = [-1 9];  % get to save space
    epoched_data = ft_selectdata(cfg, epoched_data);
    
    save(hfa_file,'epoched_hfa','epoched_data','select_events')
else
    fprintf('Loading the HFA file from: %s\n',hfa_file)
    load(hfa_file)
end

end