function SP_wordvspeak_allplot(Sbj_Metadata)

















% load channels of interest
load(fullfile(Sbj_Metadata.sbjDir,[Sbj_Metadata.sbj_ID, '_channel_OI.mat']),'channel_OI')

% select control blocks
control_blocks = select_cont_blocks(Sbj_Metadata);

% Load wavelet and ERP results
fprintf('Loading from:\n\t->%s\n',fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),[strjoin(control_blocks,'_') '_control_wltERP.mat']))
load(fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),[strjoin(control_blocks,'_') '_control_wltERP.mat']),'control_events','control_ERP','control_wlt')

%% select only channels of interest
cfg               = [];
cfg.channel       = channel_OI;
control_wlt       = ft_selectdata(cfg, control_wlt);
control_ERP       = ft_selectdata(cfg, control_ERP);

%% Baseline correct time-freq data
baseline = [-3.4 -3.1];
cfg              = [];
cfg.baseline     =  baseline;% seconds (prespeech part is 3.5455 seconds) [0.5sec + 3.0455sec]
cfg.baselinetype = 'relchange';
cfg.parameter    = 'powspctrm';
[control_wlt]    = ft_freqbaseline(cfg, control_wlt);

cfg               = [];
cfg.frequency     = [70 150];
tmp_wlt       = ft_selectdata(cfg, control_wlt);

cfg               = [];
cfg.frequency     = [0.5 12];
control_wlt       = ft_selectdata(cfg, control_wlt);

%% find word onsets, group acc to accuracy and get TF data organized
clear word_trial_nos pE_trial pR_trial pE_onset pR_onset
word_trial_nos{2}=[];word_onsets{2}=[];
for t = 1:size(control_events,1)
    for w = 1:5
        if strcmp(control_events.word_info{t}.response{w},'1')
            word_trial_nos{1}(end+1,1) = t;
            word_onsets{1}(end+1,1) = control_wlt.time(nearest(control_wlt.time, control_events.word_info{t}.onset(w)));
        elseif strcmp(control_events.word_info{t}.response{w},'0')
            word_trial_nos{2}(end+1,1) = t;
            word_onsets{2}(end+1,1) = control_wlt.time(nearest(control_wlt.time, control_events.word_info{t}.onset(w)));
        end
    end
end
accr_word_tf = SP_epoch_multi(control_wlt,word_trial_nos{1},word_onsets{1},-0.05,0.5);
wrong_word_tf = SP_epoch_multi(control_wlt,word_trial_nos{2},word_onsets{2},-0.05,0.5);
accr_word_erp = SP_epoch_multi(control_ERP,word_trial_nos{1},word_onsets{1},-0.05,0.5);
wrong_word_erp = SP_epoch_multi(control_ERP,word_trial_nos{2},word_onsets{2},-0.05,0.5);

%% find peak onsets, group acc to accuracy and get TF data organized
pE_trial{2}=[];pE_onset{2}=[];
pR_trial{2}=[];pR_onset{2}=[];

for t = 1:size(control_events,1)
    % Separate peakEnv events
    % Collect fouri in a cell structure
    for pE = 1:length(control_events.peak_info{t}.peakEnv{1})
        % first check if the points are too close to each other
        if pE == length(control_events.peak_info{t}.peakEnv{1}) % last one in the list
            % do nothing
        elseif control_events.peak_info{t}.peakEnv{1}(pE+1) - control_events.peak_info{t}.peakEnv{1}(pE) < 0.33
            continue
        end
        % check what is the response when that word was heard
        wword_resp = control_events.word_info{t}.response{...
            control_events.peak_info{t}.peakEnv{1}(pE)>=control_events.word_info{t}.onset & ...
            control_events.peak_info{t}.peakEnv{1}(pE)<control_events.word_info{t}.offset};
        if strcmp(wword_resp,'1')
            pE_trial{1}(end+1,1) = t;
            pE_onset{1}(end+1,1) = control_wlt.time(nearest(control_wlt.time, control_events.peak_info{t}.peakEnv{1}(pE)));
        elseif strcmp(wword_resp,'0')
            pE_trial{2}(end+1,1) = t;
            pE_onset{2}(end+1,1) = control_wlt.time(nearest(control_wlt.time, control_events.peak_info{t}.peakEnv{1}(pE)));
        end
    end
    % Separate peakRate like different trials
    % Collect fouri in a cell structure
    for pR = 1:length(control_events.peak_info{t}.peakRate{1})
        % first check if the points are too close to each other
        if pR == length(control_events.peak_info{t}.peakRate{1})
            % do nothing
        elseif control_events.peak_info{t}.peakRate{1}(pR+1) - control_events.peak_info{t}.peakRate{1}(pR) < 0.33
            continue
        end
        % check what is the response when that word was heard
        wword_resp = control_events.word_info{t}.response{...
            control_events.peak_info{t}.peakRate{1}(pR)>=control_events.word_info{t}.onset & ...
            control_events.peak_info{t}.peakRate{1}(pR)<control_events.word_info{t}.offset};
        if strcmp(wword_resp,'1')
            pR_trial{1}(end+1,1) = t;
            pR_onset{1}(end+1,1) = control_wlt.time(nearest(control_wlt.time, control_events.peak_info{t}.peakRate{1}(pR)));
        elseif strcmp(wword_resp,'0')
            pR_trial{2}(end+1,1) = t;
            pR_onset{2}(end+1,1) = control_wlt.time(nearest(control_wlt.time, control_events.peak_info{t}.peakRate{1}(pR)));
        end
        
    end
end
[accr_pR_tf, freq_ITPC, time_ITPC] = SP_epoch_multi(control_wlt,pR_trial{1},pR_onset{1},-0.05,0.5);
wrong_pR_tf = SP_epoch_multi(control_wlt,pR_trial{2},pR_onset{2},-0.05,0.5);
[accr_pR_erp,~,time_ERP] = SP_epoch_multi(control_ERP,pR_trial{1},pR_onset{1},-0.05,0.5);
wrong_pR_erp = SP_epoch_multi(control_ERP,pR_trial{2},pR_onset{2},-0.05,0.5);
accr_pE_tf = SP_epoch_multi(control_wlt,pE_trial{1},pE_onset{1},-0.05,0.5);
wrong_pE_tf = SP_epoch_multi(control_wlt,pE_trial{2},pE_onset{2},-0.05,0.5);
accr_pE_erp = SP_epoch_multi(control_ERP,pE_trial{1},pE_onset{1},-0.05,0.5);
wrong_pE_erp = SP_epoch_multi(control_ERP,pE_trial{2},pE_onset{2},-0.05,0.5);

%% get HFA parts
accr_word_hfa = SP_epoch_multi(tmp_wlt,word_trial_nos{1},word_onsets{1},-0.05,0.5);
wrong_word_hfa = SP_epoch_multi(tmp_wlt,word_trial_nos{2},word_onsets{2},-0.05,0.5);
accr_pR_hfa = SP_epoch_multi(tmp_wlt,pR_trial{1},pR_onset{1},-0.05,0.5);
wrong_pR_hfa = SP_epoch_multi(tmp_wlt,pR_trial{2},pR_onset{2},-0.05,0.5);
accr_pE_hfa = SP_epoch_multi(tmp_wlt,pE_trial{1},pE_onset{1},-0.05,0.5);
wrong_pE_hfa = SP_epoch_multi(tmp_wlt,pE_trial{2},pE_onset{2},-0.05,0.5);

clear tmp_wlt

%% Calculate POS statistics
% equalize number of events(trials)
if size(accr_word_tf,1) > size(wrong_word_tf,1)
    nt = randperm(size(accr_word_tf,1),size(wrong_word_tf,1));
    accr_word_tf = accr_word_tf(nt,:,:,:);
else
    nt = randperm(size(wrong_word_tf,1),size(accr_word_tf,1));
    wrong_word_tf = wrong_word_tf(nt,:,:,:);
end
data1 = permute(accr_word_tf,[2,3,4,1]);
data2 = permute(wrong_word_tf,[2,3,4,1]);
% calculate different phase opposition values, best seems to be p_zPOS
[p_circWW{1}, p_POS{1}, p_zPOS{1}] = PhaseOpposition(data1, data2, 1000, 3);

for pp=2:3 % loop peakRate and peakEnv
    % move the trials to last dimension
    if pp==2 % peakRate events
        % equalize number of events(trials)
        if size(accr_pR_tf,1) > size(wrong_pR_tf,1)
            nt = randperm(size(accr_pR_tf,1),size(wrong_pR_tf,1));
            accr_pR_tf = accr_pR_tf(nt,:,:,:);
        else
            nt = randperm(size(wrong_pR_tf,1),size(accr_pR_tf,1));
            wrong_pR_tf = wrong_pR_tf(nt,:,:,:);
        end
        data1 = permute(accr_pR_tf,[2,3,4,1]);
        data2 = permute(wrong_pR_tf,[2,3,4,1]);
    else
        % equalize number of events(trials)
        if size(accr_pE_tf,1) > size(wrong_pE_tf,1)
            nt = randperm(size(accr_pE_tf,1),size(wrong_pE_tf,1));
            accr_pE_tf = accr_pE_tf(nt,:,:,:);
        else
            nt = randperm(size(wrong_pE_tf,1),size(accr_pE_tf,1));
            wrong_pE_tf = wrong_pE_tf(nt,:,:,:);
        end
        data1 = permute(accr_pE_tf,[2,3,4,1]);
        data2 = permute(wrong_pE_tf,[2,3,4,1]);
    end
    % calculate different phase opposition values, best seems to be p_zPOS
    [p_circWW{pp}, p_POS{pp}, p_zPOS{pp}] = PhaseOpposition(data1, data2, 1000, 3);
end


%% Calculate ITPCz for each conditions
itpc = [];

% for accurate words
tmp      = accr_word_tf./abs(accr_word_tf);   % divide by amplitude
tmp      = sum(tmp,1);                            % sum angles across trials
tmp      = abs(tmp)/size(accr_word_tf,1);       % take the absolute value and normalize
itpc(1,1,:,:,:) = size(accr_word_tf,1) * squeeze(tmp).^2; % remove the first singleton dimension and n*ITPC^2
% for wrong words
tmp      = wrong_word_tf./abs(wrong_word_tf);   % divide by amplitude
tmp      = sum(tmp,1);                            % sum angles across trials
tmp      = abs(tmp)/size(wrong_word_tf,1);       % take the absolute value and normalize
itpc(1,2,:,:,:) = size(wrong_word_tf,1) * squeeze(tmp).^2; % remove the first singleton dimension and n*ITPC^2

for pp = 2:3
    for cond = 1:2
        % get data
        if cond == 1 && pp == 2
            curr_fouri_all = accr_pR_tf;
        elseif cond == 2 && pp == 2
            curr_fouri_all = wrong_pR_tf;
        elseif cond == 1 && pp == 3
            curr_fouri_all = accr_pE_tf;
        elseif cond == 2 && pp == 3
            curr_fouri_all = wrong_pE_tf;
        end
        
        % compute inter-trial phase coherence (itpc) for each conditions
        tmp      = curr_fouri_all./abs(curr_fouri_all);   % divide by amplitude
        tmp      = sum(tmp,1);                            % sum angles across trials
        tmp      = abs(tmp)/size(curr_fouri_all,1);       % take the absolute value and normalize
        itpc(pp,cond,:,:,:) = size(curr_fouri_all,1) * squeeze(tmp).^2; % remove the first singleton dimension and n*ITPC^2
    end
end

%% The greatest plot so far
save_folder = fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'), 'wordvspeak_allplot');
if ~exist(save_folder,'dir'),mkdir(save_folder),end
bwr = load('bwr_cmap.mat');

plh = 3; % plot height: word onset locked, peakRate, peakEnv
plw = 5; % plot width: ERP, HFA, ITPC-correct, ITPC-noresp, POS

for el = 1:length(control_wlt.label)
    curr_label = control_wlt.label{el};
    
    figure('Units','normalized','Position', [0 0  1 1]);
    
    % plot single trials and average of corr&no - ERP - word
    subplot(plh,plw,subplotno(plw,1,1))
    for_avg=squeeze(accr_pR_erp(:,el,:));
    hold on
    shadedErrorBar(time_ERP,mean(for_avg,1),stderr(for_avg),'lineprops','b');
    for_avg=squeeze(wrong_pR_erp(:,el,:));
    shadedErrorBar(time_ERP,mean(for_avg,1),stderr(for_avg),'lineprops','r');
    title({'ERP (mean+stderr)';[num2str(size(accr_pR_erp,1)) ' correct response - ' num2str(size(wrong_pR_erp,1)) ' no response']})
    xlim([time_ERP(1) time_ERP(end)])
    set(gca, 'FontSize',12);
    plot([0 0], ylim,'k')
    plot(xlim,[0 0],'k')
    legend({'Correct','No resp'})
    
    % plot single trials and average of corr&no - ERP - pR
    subplot(plh,plw,subplotno(plw,1,1))
    for_avg=squeeze(accr_pR_erp(:,el,:));
    hold on
    shadedErrorBar(time_ERP,mean(for_avg,1),stderr(for_avg),'lineprops','b');
    for_avg=squeeze(wrong_pR_erp(:,el,:));
    shadedErrorBar(time_ERP,mean(for_avg,1),stderr(for_avg),'lineprops','r');
    title({'ERP (mean+stderr)';[num2str(size(accr_pR_erp,1)) ' correct response - ' num2str(size(wrong_pR_erp,1)) ' no response']})
    xlim([time_ERP(1) time_ERP(end)])
    set(gca, 'FontSize',12);
    plot([0 0], ylim,'k')
    plot(xlim,[0 0],'k')
    legend({'Correct','No resp'})
    
    % plot single trials and average of corr&no - ERP - pE
    subplot(plh,plw,subplotno(plw,1,1))
    for_avg=squeeze(accr_pR_erp(:,el,:));
    hold on
    shadedErrorBar(time_ERP,mean(for_avg,1),stderr(for_avg),'lineprops','b');
    for_avg=squeeze(wrong_pR_erp(:,el,:));
    shadedErrorBar(time_ERP,mean(for_avg,1),stderr(for_avg),'lineprops','r');
    title({'ERP (mean+stderr)';[num2str(size(accr_pR_erp,1)) ' correct response - ' num2str(size(wrong_pR_erp,1)) ' no response']})
    xlim([time_ERP(1) time_ERP(end)])
    set(gca, 'FontSize',12);
    plot([0 0], ylim,'k')
    plot(xlim,[0 0],'k')
    legend({'Correct','No resp'})
    
    
    
    
    
    % Loop to create the ITPC and POs plots
    for pp = 1:3 % first word onset, second peakRate, third peakEnv
        for cond = 1:2 % loop columns: first correct, second no response
            
            subplot(plh, plw, subplotno(plw,pp,2+cond))
            surf(time_ITPC, freq_ITPC, squeeze(itpc(pp,cond,el,:,:)));set(gcf,'renderer','zbuffer');
            axis xy
            set(gca, 'FontSize',13,'FontWeight','bold');
            caxis([0 5])
            view(0,90); axis tight;
            
            hold on
            shading interp;
            %             clim_v = get(gca,'clim');  %colorbar;
            plot3([0 0],ylim,[15 15],'k');
            if cond == 1 && pp == 1
                title({'Correct response events';['Random ' num2str(size(accr_word_tf,1)) ' events']});
            elseif cond == 2 && pp == 1
                title({'Correct response events';['Random ' num2str(size(wrong_word_tf,1)) ' events']});
%             elseif cond == 1 && pp == 2
%                 title(['Random ' num2str(size(accr_pR_tf,1)) ' events during correct responses']);
%             elseif cond == 2 && pp == 2
%                 title(['Random ' num2str(size(wrong_pR_tf,1)) ' events during no-responses']);
%             elseif cond == 1 && pp == 3
%                 title(['Random ' num2str(size(accr_pE_tf,1)) ' events during correct responses']);
%             elseif cond == 2 && pp == 3
%                 title(['Random ' num2str(size(wrong_pE_tf,1)) ' events during no-responses']);
            end
            
            if cond==1 && pp == 1
                ylabel({'Word onset locked';'Frequency (Hz)'});
            end
            if cond==1 && pp == 2
                ylabel({'peakRate locked';'Frequency (Hz)'});
            end
            if cond==1 && pp == 3
                ylabel({'peakEnv locked';'Frequency (Hz)'});
            end
            if pp==3
                xlabel('Time (s)');
            end
            
        end
        
        % plot raw p-values
        subplot(plh,plw,subplotno(plw,pp,5))
        surf(time_ITPC, freq_ITPC, -log10(squeeze(p_zPOS{pp}(el,:,:))));set(gcf,'renderer','zbuffer');
        axis xy
        set(gca, 'FontSize',13,'FontWeight','bold');
        hold on
        shading interp;
        view(0,90); axis tight;
        yl=ylim;
        ylim([1 yl(2)])
        title('-log10(p-values) of z-scored POS')
        clim_p = [0 3];
        set(gca,'CLim',clim_p)  %colorbar;
        plot3([0 0],ylim,[15 15],'k');
        %         signplot = double(squeeze(p_zPOS{pp}(el,:,:)>p_thresh));
        signplot = double(squeeze(p_zPOS{pp}(el,:,:))<0.05);
        contour3(time_ITPC, freq_ITPC,15*signplot,1,'LineColor','k','LineWidth',3)
        if pp==1
            xlabel('Time (s)');
        end
        
    end
    
    
    
    
    
    
    
    
    
    
    
    
end
end
    
    
    
    
    
