% this is in case I need to bulk run the re-referencing and wavelet
% analysis from the beginning

%% Prepare
data_root = '/media/sakkol/HDD1/HBML/';
project_name = 'Speech_Perception';

AllBlockInfo = readtable(fullfile(data_root,'PROJECTS_DATA',project_name,[project_name '_BlockInfo.xlsx']));
subjects = unique(AllBlockInfo.sbj_ID);

[indx,~] = listdlg('ListString',subjects);

for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    
    for t = 1:length(Sbj_Metadata.BlockLists)
        % Get params directly from BlockList excel sheet
        curr_block = Sbj_Metadata.BlockLists{t}
        params = create_Params(Sbj_Metadata,curr_block)
        
        
        load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog.mat']),'ecog');
        
        % Notch filter, demean
        cfg=[];
        cfg.demean         = 'yes';
        cfg.bsfilter       = 'yes';
        cfg.bsfiltord      = 3;
        cfg.bsfreq         = [59 61; 119 121; 179 181];
        cont_notched = ft_preprocessing(cfg,ecog.ftrip);
        
        
        %% Re-reference
        % Average ref
        ecog.ftrip = cont_notched; % nothched or not-noteched
        plot_stuff=0;
        ignore_szr_chans=1;
        ecog_avg=ecog_avg_ref(ecog,plot_stuff,ignore_szr_chans);
        save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_avg.mat']),'ecog_avg','-v7.3');
        
        % Also BIPOLAR reference
        ecog_bp=ecog_bipolarize(ecog);
        save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_bp.mat']),'ecog_bp','-v7.3');
        
        clear ecog_bp
        
        %% Wavelet analysis
        load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']))
        events = info.events; clear info
        
        % First resample for saving space and memory
        cfg             = [];
        cfg.resamplefs  = 1000;
        ecog_avg.ftrip  = ft_resampledata(cfg,ecog_avg.ftrip);
        fs              = ecog_avg.ftrip.fsample;
        
        % Make trial structure
        % speech onset locked
        pre  = 6; % seconds (prespeech part is 3.5455 seconds) [0.5sec + 3.0455sec] - to not include any NaNs in wavelet increase by 2 seconds (from 4 to 6)
        post = 8; % seconds (longest trial is ~8.5 seconds)
        trl           = [];
        trl(:,1)      = floor( events.speech_onsets*fs - fs*pre );
        trl(:,2)      = floor( events.speech_onsets*fs + fs*post );
        trl(:,3)      = floor( -pre*fs );
        
        % Epoch
        cfg      = [];
        cfg.trl  = trl;
        epoched_data    = ft_redefinetrial(cfg,ecog_avg.ftrip); clear ecog_avg % save some memory
        % replace NaNs with zeros in first trial (if trial started fast)
        epoched_data.trial{1}(isnan(epoched_data.trial{1})) = 0;
        
        % def preproc parameters
        cfg_preproc                     = [];
        cfg_preproc.channel             = 'all';
        % cfg_preproc.padding             = 4;
        % cfg_preproc.demean              = 'yes';
        cfg_preproc.detrend             = 'yes';
        cfg_preproc.dftfilter           = 'yes';
        cfg_preproc.dftfreq             = [60 120 180];
        % cfg_preproc.baselinewindow      = [-.5 -.05];
        epoched_data                     = ft_preprocessing(cfg_preproc, epoched_data);
        
        % compute trials
        cfg             = [];
        cfg.keeptrials  = 'yes';
        epoched_data     = ft_timelockanalysis(cfg,epoched_data);
        
        if sum(sum(sum(isnan(epoched_data.trial))))
            warning('\n\n\n\t\t\tTHERE WERE %s TRIALS THAT HAVE NANs, REPLACING WITH ZEROS',sum(sum(sum(isnan(epoched_data.trial)))))
            % replace NaNs with zeros
            epoched_data.trial(isnan(epoched_data.trial)) = 0;
        end
        
        
        % Wavelet
        cfg                   = [];
        cfg.method            = 'wavelet';
        freq                  = [1 200];
        nf                    = length(freq(1):3:freq(2));
        cfg.foi               = logspace(log10(freq(1)),log10(freq(2)),nf);
        cfg.width             = 6;
        cfg.toi               = -6:0.01:8;
        cfg.keeptrials        = 'yes';
        cfg.keeptaper         = 'yes';
        cfg.output            = 'fourier';
        
        cfg.outputfile = fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_wlt_tmp.mat']); % save wavelet tmp file
        % epoched_wlt.wlt       = ft_freqanalysis(cfg,epoched_wlt);
        ft_freqanalysis(cfg,epoched_data);
        % epoched_wlt.wlt       = addRayleigh_to_ft(epoched_wlt.wlt); Rayleigh can
        % be added for control only. But because it combines all trials, it is not
        % clever to do it here.
        
        
        % read tmp file and curb 1-1 sec from ends
        cfg=[];
        cfg.latency        = [-4 6];
        cfg.inputfile = fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_wlt_tmp.mat']);
        epoched_wlt=ft_selectdata(cfg);
        
        % save wlt and data and delete tmp
        save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_wlt.mat']),'epoched_wlt','epoched_data','-v7.3');
        delete(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_wlt_tmp.mat']));
    end
end
