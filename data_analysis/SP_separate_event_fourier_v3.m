function SP_separate_event_fourier_v3(Sbj_Metadata,freq_bands)
% This separate automatically the fourier spectrum of each event based on
% different frequency bands (can be defined) and patients response
% (correct, wrong and no response).

% Which freq bands
vars = who;
if ~ismember(vars,'freq_bands')
    freq_bands = {[1 4],[4 8],[8 12],[70 150]};
end
clear vars

% Select blocks to import
control_blocks = select_cont_blocks(Sbj_Metadata);

%% bring in these blocks and combine only the control events
fprintf('For subject %s, these blocks are going to be used:%s\n',Sbj_Metadata.sbj_ID,strjoin(control_blocks,', '))
for b = 1:length(control_blocks)
    curr_block = control_blocks{b};
    fprintf('...loading:%s\n',curr_block)
    % Load iEEG
    load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_wlt.mat']))
    load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_info.mat']))
    events = info.events; clear info
    
    % Select only control events
    control_idx = events.Cond_code == 1;
    
    % Select ecog data
    cfg = [];
    cfg.trials = control_idx;
    [epoched_wlt.wlt] = ft_selectdata(cfg, epoched_wlt.wlt);
    
    % Select events
    events = events(control_idx,:);
    
    % Append to overall list
    if b == 1
        control_wlt = epoched_wlt.wlt;
        control_events = events;
    else
        cfg = [];
        cfg.parameter  = 'fourierspctrm';
        control_wlt = ft_appendfreq(cfg, control_wlt, epoched_wlt.wlt);
        
        control_events = [control_events;events];
    end
    
end

%% Separate fourier spectrums of correct responses and others
fouri_of_words = [];
fouri_of_words.freq_bands = freq_bands;

fouri_of_words.corr_rspn_fouri_peakEnv = {[]};
fouri_of_words.no_rspn_fouri_peakEnv = {[]};
fouri_of_words.wrng_rspn_fouri_peakEnv = {[]};

fouri_of_words.corr_rspn_fouri_peakRate = {[]};
fouri_of_words.no_rspn_fouri_peakRate = {[]};
fouri_of_words.wrng_rspn_fouri_peakRate = {[]};


fprintf('There are total of %d events\n',size(control_events,1))
for f = 1:length(freq_bands)
    fprintf('Separating data for frequency range: %d-%dHz\n',freq_bands{f}(1),freq_bands{f}(2))
    
    for ii = freq_bands{1}(1):freq_bands{1}(2)
        cl_freqs(ii) = nearest(control_wlt.freq, ii);
    end
    
    fouri_of_words.freq_band_dtls{f} = control_wlt.freq(cl_freqs);
    
    wcorr_ind = 1;
    wno_ind = 1;
    wwrng_ind = 1;
    scorr_ind = 1;
    sno_ind = 1;
    swrng_ind = 1;

    fprintf('Current event:')
    for t = 1:size(control_events,1)
        fprintf('-%d',t)
        
        % Separate peakEnv events
        % Collect fouri in a cell structure
        for pE = 1:length(control_events.peak_info{t}.peakEnv{1})
            
            % first check if the points are too close to each other
            if pE == length(control_events.peak_info{t}.peakEnv{1})
                % do nothing
            elseif control_events.peak_info{t}.peakEnv{1}(pE+1) - control_events.peak_info{t}.peakEnv{1}(pE) < 0.25
                continue
            end
            
            % check what is the response when that word was heard
            wword_resp = control_events.word_info{t}.response{...
                control_events.peak_info{t}.peakEnv{1}(pE)>=control_events.word_info{t}.onset & ...
                control_events.peak_info{t}.peakEnv{1}(pE)<control_events.word_info{t}.offset};
            
            % find closest timepoint
            cl_times = [nearest(control_wlt.time, control_events.peak_info{t}.peakEnv{1}(pE)-0.05)...
                        nearest(control_wlt.time,control_events.peak_info{t}.peakEnv{1}(pE)+.5)];
            if strcmp(wword_resp,'1')
                fouri_of_words.corr_rspn_fouri_peakEnv{1,f}(wcorr_ind,:,:,:) = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs,cl_times(1):cl_times(2)));
                wcorr_ind = wcorr_ind+1;
            elseif strcmp(wword_resp,'0')
                fouri_of_words.no_rspn_fouri_peakEnv{1,f}(wno_ind,:,:,:) = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs,cl_times(1):cl_times(2)));
                wno_ind = wno_ind+1;
                
            else % wrong responses
                fouri_of_words.wrng_rspn_fouri_peakEnv{1,f}(wwrng_ind,:,:,:) = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs,cl_times(1):cl_times(2)));
                wwrng_ind = wwrng_ind+1;
                
            end
            
        end
        
        % Separate peakRate like different trials
        % Collect fouri in a cell structure
        for pR = 1:length(control_events.peak_info{t}.peakRate{1})
            
            % first check if the points are too close to each other
            if pR == length(control_events.peak_info{t}.peakRate{1})
                % do nothing
            elseif control_events.peak_info{t}.peakRate{1}(pR+1) - control_events.peak_info{t}.peakRate{1}(pR) < 0.25
                continue
            end
            
            % check what is the response when that word was heard
            wword_resp = control_events.word_info{t}.response{...
                control_events.peak_info{t}.peakRate{1}(pR)>=control_events.word_info{t}.onset & ...
                control_events.peak_info{t}.peakRate{1}(pR)<control_events.word_info{t}.offset};
            
            % find closest timepoint
            cl_times = [nearest(control_wlt.time, control_events.peak_info{t}.peakRate{1}(pR)-0.05)...
                        nearest(control_wlt.time, control_events.peak_info{t}.peakRate{1}(pR)+.5)];
                    
            if strcmp(wword_resp,'1')
                fouri_of_words.corr_rspn_fouri_peakRate{1,f}(scorr_ind,:,:,:) = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs,cl_times(1):cl_times(2)));
                scorr_ind = scorr_ind+1;
            elseif strcmp(wword_resp,'0')
                fouri_of_words.no_rspn_fouri_peakRate{1,f}(sno_ind,:,:,:) = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs,cl_times(1):cl_times(2)));
                sno_ind = sno_ind+1;
            else % wrong responses
                fouri_of_words.wrng_rspn_fouri_peakRate{1,f}(swrng_ind,:,:,:) = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs,cl_times(1):cl_times(2)));
                swrng_ind = swrng_ind+1;
                
            end
        end

    end
    fprintf('\n\n')
end

fouri_of_words.time_dtls = linspace(-0.05,0.5,size(fouri_of_words.corr_rspn_fouri_peakEnv{1},4));

%% Save v3
save_dir = fullfile(Sbj_Metadata.results, [strjoin(control_blocks,'_') '_v3']);
if ~exist(save_dir,'dir'),mkdir(save_dir),end
fprintf('Saving ''fouri_of_words'' to:\n-->%s\n',fullfile(save_dir, [strjoin(control_blocks,'_') '_ctrl_word_fouri.mat']))
save(fullfile(save_dir, [strjoin(control_blocks,'_') '_ctrl_word_fouri.mat']),'fouri_of_words','-v7.3');

end
