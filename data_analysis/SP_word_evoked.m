function SP_word_evoked(Sbj_Metadata,control_blocks,HFAvsERP)
% This is to get the mean of evoked responses (ERP and/or HFA) word heard
% vs others.

%% Load wlt
% Select blocks to import
vars=who;
if ~ismember(vars,'control_blocks')
    control_blocks = select_cont_blocks(Sbj_Metadata);
end
clear vars

% bring in these blocks and combine only the control events
fprintf('These blocks are going to be used: %s\n',strjoin(control_blocks,', '))
for b = 1:length(control_blocks)
    curr_block = control_blocks{b};
    fprintf('...loading:%s\n',curr_block)
    
    % Load iEEG
    load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_wlt.mat']))
    events = epoched_wlt.events;
    
    % Select only control events
    control_idx = events.Cond_code == 1;
    
    % Select ecog data
    cfg = [];
    cfg.trials = control_idx;
    [epoched_wlt.wlt] = ft_selectdata(cfg, epoched_wlt.wlt);
    
    % from fourierspectrum to powerspectrum
    cfg = [];
    cfg.output='abs';
    cfg.keeptrials = 'yes';
    epoched_wlt.wlt=ft_freqdescriptives(cfg,epoched_wlt.wlt);

    % Select events
    events = events(control_idx,:);
    
    % Append to overall list
    if b == 1
        control_wlt = epoched_wlt.wlt;
        control_events = events;
    else
        cfg = [];
        cfg.parameter  = 'powspctrm';
        control_wlt = ft_appendfreq(cfg, control_wlt, epoched_wlt.wlt);
        control_events = [control_events;events];
    end
    clear epoched_wlt events
end

%%  Baseline correct time-freq data
cfg              = [];
cfg.baseline     = [-3.45 -3.05]; % seconds (prespeech part is 3.5455 seconds) [0.5sec + 3.0455sec]
cfg.baselinetype = 'db';
cfg.parameter    = 'powspctrm';
[wlt_bc]         = ft_freqbaseline(cfg, control_wlt);

%% Separate conditions

fprintf('There are total of %d events\n',size(control_events,1))


for c = 1:length(HFAvsERP)
    fprintf('\nCurrent analysis:')
    
    
    fprintf('-%d',t)
    for w = 1:5
        % Separate words
        % Check if correct
        % Collect fouri in a cell structure
        if strcmp(control_events.word_info{t}.response{w},'1')
            % find closest timepoint
            cl_times = [nearest(control_wlt.time, control_events.word_info{t}.onset(w))...
                nearest(control_wlt.time,control_events.word_info{t}.offset(w))];
            corr_rspn_fouri_word{wcorr_ind,f} = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)));
            wcorr_ind = wcorr_ind+1;
        elseif strcmp(control_events.word_info{t}.response{w},'0')
            % find closest timepoint
            cl_times = [nearest(control_wlt.time, control_events.word_info{t}.onset(w))...
                nearest(control_wlt.time,control_events.word_info{t}.offset(w))];
            no_rspn_fouri_word{wno_ind,f} = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)));
            wno_ind = wno_ind+1;
            
        else % wrong responses
            % find closest timepoint
            cl_times = [nearest(control_wlt.time, control_events.word_info{t}.onset(w))...
                nearest(control_wlt.time,control_events.word_info{t}.offset(w))];
            wrng_rspn_fouri_word{wwrng_ind,f} = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)));
            wwrng_ind = wwrng_ind+1;
            
        end
        
    end
    
    for s = 1:size(control_events.syllable_info{t},1)
        % Separate syllables like different trials
        % Check if correct
        % Collect fouri in a cell structure
        if strcmp(control_events.word_info{t}.response{control_events.syllable_info{t}.word_id(s)},'1')
            % find closest timepoint
            cl_times = [nearest(control_wlt.time, control_events.syllable_info{t}.onset(s))...
                nearest(control_wlt.time,control_events.syllable_info{t}.offset(s))];
            corr_rspn_fouri_syll{scorr_ind,f} = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)));
            scorr_ind = scorr_ind+1;
        elseif strcmp(control_events.word_info{t}.response{control_events.syllable_info{t}.word_id(s)},'0')
            % find closest timepoint
            cl_times = [nearest(control_wlt.time, control_events.syllable_info{t}.onset(s))...
                nearest(control_wlt.time,control_events.syllable_info{t}.offset(s))];
            no_rspn_fouri_syll{sno_ind,f} = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)));
            sno_ind = sno_ind+1;
            
        else % wrong responses
            % find closest timepoint
            cl_times = [nearest(control_wlt.time, control_events.syllable_info{t}.onset(s))...
                nearest(control_wlt.time,control_events.syllable_info{t}.offset(s))];
            wrng_rspn_fouri_syll{swrng_ind,f} = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)));
            swrng_ind = swrng_ind+1;
            
        end
    end
    
    for p = 1:size(control_events.all_info{t},1)
        % Separate phonemes like different trials
        % Check if correct
        % Collect fouri in a cell structure
        if strcmp(control_events.word_info{t}.response{control_events.all_info{t}.word_id(p)},'1')
            % find closest timepoint
            cl_times = [nearest(control_wlt.time, control_events.syllable_info{t}.onset(s))...
                nearest(control_wlt.time,control_events.syllable_info{t}.offset(s))];
            corr_rspn_fouri_pho{pcorr_ind,f} = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)));
            pcorr_ind = pcorr_ind+1;
        elseif strcmp(control_events.word_info{t}.response{control_events.all_info{t}.word_id(p)},'0')
            % find closest timepoint
            cl_times = [nearest(control_wlt.time, control_events.syllable_info{t}.onset(s))...
                nearest(control_wlt.time,control_events.syllable_info{t}.offset(s))];
            no_rspn_fouri_pho{pno_ind,f} = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)));
            pno_ind = pno_ind+1;
            
        else % wrong responses
            % find closest timepoint
            cl_times = [nearest(control_wlt.time, control_events.all_info{t}.onset(s))...
                nearest(control_wlt.time,control_events.syllable_info{t}.offset(s))];
            wrng_rspn_fouri_pho{pwrng_ind,f} = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)));
            pwrng_ind = pwrng_ind+1;
            
        end
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
end


end