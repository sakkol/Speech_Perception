function [all_output,freq,time] = SP_epoch_multi(ft_data_TF,trial_nos,onsets,before,after)
% this is to get segments from ERP or time-frequency data where some
% segments can be in different time points and trials.

if isfield(ft_data_TF,'fourierspctrm')
    lchan = size(ft_data_TF.fourierspctrm,2);
    fieldOI = 'fourierspctrm';
elseif isfield(ft_data_TF,'powspctrm') && length(size(ft_data_TF.powspctrm)) == 4
    lchan = size(ft_data_TF.powspctrm,2);
    fieldOI = 'powspctrm';
elseif isfield(ft_data_TF,'powspctrm') && length(size(ft_data_TF.powspctrm)) == 3
    lchan = size(ft_data_TF.powspctrm,1);
    fieldOI = 'powspctrm';
else % probably time domain ERP
    lchan = size(ft_data_TF.trial,2);
end

% loop trial_nos
fprintf('Extracting data from %d trials\n',length(trial_nos))
for t = 1:length(trial_nos)
    % get -before and +after relative to onsets use ft_selectdata
    cfg=[];
    cfg.trials = trial_nos(t);
    cfg.latency = [onsets(t)+before,onsets(t)+after];
    tmp = ft_selectdata(cfg,ft_data_TF);
        
    if isfield(ft_data_TF,'powspctrm') || isfield(ft_data_TF,'fourierspctrm')
        if t==1;all_output = zeros([length(trial_nos) lchan length(tmp.freq) length(tmp.time)]);end
        all_output(t,:,:,:) = tmp.(fieldOI);
        freq = tmp.freq;
    else
        if t==1;all_output = zeros([length(trial_nos) lchan length(tmp.time)]);end
        all_output(t,:,:) = tmp.trial;
        freq = 0;
    end
end
time = linspace(before,after,length(tmp.time));
end