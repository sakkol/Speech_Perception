function [epoched_hfa,epoched_data] = SP_getAttHFA(Sbj_Metadata, curr_block)
% just to grab the HFA file easily
hfa_file = fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_AttHFA.mat']);

if ~exist(hfa_file,'file')
    fprintf('Couldn''t find the HFA file, running the analysis and saving to: %s\n', hfa_file)
    tmp=load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_info.mat']));info=tmp.info;clear tmp
    load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_avg.mat']))
    % run wavelet
    pre  = 1;
    post = 3.5;
    freq = 70:5:200;freq(ismember(freq,[120 180]))=[];
    fourier_pow = 'pow';
    [epoched_data, epoched_wlt] = master_waveletTF(ecog_avg.ftrip,info.events.att_sent_onset,pre,post,freq,fourier_pow);
    
    % run HFA
    cfg=[];
    cfg.frequency=[69 201];
    cfg.avgoverfreq='no';
    epoched_hfa = ft_selectdata(cfg, epoched_wlt);
    cfg              = [];
    cfg.baseline     = [-.95 -.55];
    cfg.baselinetype = 'relative';
    cfg.parameter    = 'powspctrm';
    [epoched_hfa]    = ft_freqbaseline(cfg, epoched_hfa);
    cfg=[];
    cfg.frequency=[70 200];
    cfg.avgoverfreq='yes';
    epoched_hfa = ft_selectdata(cfg, epoched_hfa);
    
    save(hfa_file,'epoched_hfa','epoched_data')
else
    fprintf('Loading the HFA file from: %s\n',hfa_file)
    load(hfa_file)
end

end