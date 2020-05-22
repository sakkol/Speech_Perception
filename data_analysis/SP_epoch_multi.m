function [all_output,freq,time] = SP_epoch_multi(ft_data_TF,trial_nos,onsets,before,after)
% this is to get segments from ERP or time-frequency data where some
% segments can be in different time points and trials.

if isfield(ft_data_TF,'powspctrm')
    lchan = size(ft_data_TF.powspctrm,2);
    lsize = length(size(ft_data_TF.powspctrm));
else
    lchan = size(ft_data_TF.trial,2);
    lsize = size(ft_data_TF.trial);
end

% loop trial_nos
fprintf('Extracting data from %d trials\n',length(trial_nos))
for t = 1:length(trial_nos)
    % get -before and +after relative to onsets use ft_selectdata
    cfg=[];
    cfg.trials = trial_nos(t);
    cfg.latency = [onsets(t)+before,onsets(t)+after];
    tmp = ft_selectdata(cfg,ft_data_TF);
        
    if lsize == 4
        if t==1;all_output = zeros([length(trial_nos) lchan length(tmp.freq) length(tmp.time)]);end
        all_output(t,:,:,:) = tmp.powspctrm;
        freq = tmp.freq;
    else
        if t==1;all_output = zeros([length(trial_nos) lchan length(tmp.time)]);end
        all_output(t,:,:) = tmp.trial;
        freq = 0;
    end
end
time = linspace(before,after,length(tmp.time));
end