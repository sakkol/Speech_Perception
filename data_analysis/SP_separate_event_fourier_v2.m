function SP_separate_event_fourier_v2(Sbj_Metadata,freq_bands)
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

fouri_of_words.corr_rspn_fouri_word = {[]};
fouri_of_words.no_rspn_fouri_word = {[]};
fouri_of_words.wrng_rspn_fouri_word = {[]};

fouri_of_words.corr_rspn_fouri_syll = {[]};
fouri_of_words.no_rspn_fouri_syll = {[]};
fouri_of_words.wrng_rspn_fouri_syll = {[]};

fouri_of_words.corr_rspn_fouri_pho = {[]};
fouri_of_words.no_rspn_fouri_pho = {[]};
fouri_of_words.wrng_rspn_fouri_pho = {[]};


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
    pcorr_ind = 1;
    pno_ind = 1;
    pwrng_ind = 1;
    fprintf('Current event:')
    for t = 1:size(control_events,1)
        fprintf('-%d',t)
        for w = 1:5
            % Separate words
            % Collect fouri in a cell structure
            
            % find closest timepoint
            cl_times = [nearest(control_wlt.time, control_events.word_info{t}.onset(w)-0.05)...
                nearest(control_wlt.time,control_events.word_info{t}.onset(w)+.7)];
            if strcmp(control_events.word_info{t}.response{w},'1')
                fouri_of_words.corr_rspn_fouri_word{1,f}(wcorr_ind,:,:,:) = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs,cl_times(1):cl_times(2)));
                wcorr_ind = wcorr_ind+1;
            elseif strcmp(control_events.word_info{t}.response{w},'0')
                fouri_of_words.no_rspn_fouri_word{1,f}(wno_ind,:,:,:) = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs,cl_times(1):cl_times(2)));
                wno_ind = wno_ind+1;
                
            else % wrong responses
                fouri_of_words.wrng_rspn_fouri_word{1,f}(wwrng_ind,:,:,:) = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs,cl_times(1):cl_times(2)));
                wwrng_ind = wwrng_ind+1;
                
            end
            
        end
        
        for s = 1:size(control_events.syllable_info{t},1)
            % Separate syllables like different trials
            % Collect fouri in a cell structure
            % find closest timepoint
            cl_times = [nearest(control_wlt.time, control_events.syllable_info{t}.onset(s)-0.05)...
                nearest(control_wlt.time,control_events.syllable_info{t}.onset(s)+.7)];
            if strcmp(control_events.word_info{t}.response{control_events.syllable_info{t}.word_id(s)},'1')
                fouri_of_words.corr_rspn_fouri_syll{1,f}(scorr_ind,:,:,:) = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs,cl_times(1):cl_times(2)));
                scorr_ind = scorr_ind+1;
            elseif strcmp(control_events.word_info{t}.response{control_events.syllable_info{t}.word_id(s)},'0')
                fouri_of_words.no_rspn_fouri_syll{1,f}(sno_ind,:,:,:) = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs,cl_times(1):cl_times(2)));
                sno_ind = sno_ind+1;
            else % wrong responses
                fouri_of_words.wrng_rspn_fouri_syll{1,f}(swrng_ind,:,:,:) = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs,cl_times(1):cl_times(2)));
                swrng_ind = swrng_ind+1;
                
            end
        end
        
        for p = 1:size(control_events.all_info{t},1)
            % Separate phonemes like different trials
            % Collect fouri in a cell structure
            % find closest timepoint
            cl_times = [nearest(control_wlt.time, control_events.all_info{t}.phoneme_onset(p)-0.05)...
                nearest(control_wlt.time,control_events.all_info{t}.phoneme_onset(p)+.7)];
            if strcmp(control_events.word_info{t}.response{control_events.all_info{t}.word_id(p)},'1')
                fouri_of_words.corr_rspn_fouri_pho{1,f}(pcorr_ind,:,:,:) = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs,cl_times(1):cl_times(2)));
                pcorr_ind = pcorr_ind+1;
            elseif strcmp(control_events.word_info{t}.response{control_events.all_info{t}.word_id(p)},'0')
                fouri_of_words.no_rspn_fouri_pho{1,f}(pno_ind,:,:,:) = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs,cl_times(1):cl_times(2)));
                pno_ind = pno_ind+1;
            else % wrong responses
                fouri_of_words.wrng_rspn_fouri_pho{1,f}(pwrng_ind,:,:,:) = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs,cl_times(1):cl_times(2)));
                pwrng_ind = pwrng_ind+1;
                
            end
        end
    end
    fprintf('\n\n')
end

fouri_of_words.time_dtls = control_wlt.time(cl_times(1):cl_times(2));
%% Save v2
save_dir = fullfile(Sbj_Metadata.results, [strjoin(control_blocks,'_') '_v2']);
if ~exist(save_dir,'dir'),mkdir(save_dir),end
fprintf('Saving ''fouri_of_words'' to:\n-->%s\n',fullfile(save_dir, [strjoin(control_blocks,'_') '_ctrl_word_fouri.mat']))
save(fullfile(save_dir, [strjoin(control_blocks,'_') '_ctrl_word_fouri.mat']),'fouri_of_words','-v7.3');

end
