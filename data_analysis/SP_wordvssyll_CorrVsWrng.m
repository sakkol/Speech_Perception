function SP_wordvssyll_CorrVsWrng(Sbj_Metadata)
% Similar to SP_wordvspeak_allplot function, but only plots word and
% syllable onset locked events and compares correct and wrong responses.
% This function loads control blocks (wavelet and ERP results), focusing on
% channels of interest, it calculates ERP-HFA-ITPC of correct vs wrong
% response trials, also calculates POS values to compare between these 2
% conditions and importantly, all of these is done for word onset and
% syllable onset locked events. Plot consists of 2 rows (different onset
% locks) and 5 columns (ERP-HFA-ITPC(correct)-ITPC(wrong response)-POS).
% Serdar AKKOL, HBML, May 2020

%% Preparations
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

%% arrange a little bit
% from fourierspectrum to powerspectrum for HFA
cfg                = [];
cfg.output         = 'abs';
cfg.keeptrials     = 'yes';
hfa_wlt            = ft_freqdescriptives(cfg,control_wlt);

% Baseline correct time-freq data for HFA
baseline = [-3.4 -3.1];
cfg              = [];
cfg.baseline     =  baseline;% seconds (prespeech part is 3.5455 seconds) [0.5sec + 3.0455sec]
cfg.baselinetype = 'relative';
cfg.parameter    = 'powspctrm';
[hfa_wlt]    = ft_freqbaseline(cfg, hfa_wlt);
cfg                = [];
cfg.frequency      = [70 150];
cfg.avgoverfreq    = 'yes';
hfa_wlt            = ft_selectdata(cfg, hfa_wlt);
hfa_wlt.powspctrm = smoothdata(hfa_wlt.powspctrm,4,'gaussian',10);

% to save space and to run stats only in low frequency
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
        elseif ~strcmp(control_events.word_info{t}.response{w},'0')
            word_trial_nos{2}(end+1,1) = t;
            word_onsets{2}(end+1,1) = control_wlt.time(nearest(control_wlt.time, control_events.word_info{t}.onset(w)));
        end
    end
end
accr_word_tf = SP_epoch_multi(control_wlt,word_trial_nos{1},word_onsets{1},-0.05,0.5);
wrong_word_tf = SP_epoch_multi(control_wlt,word_trial_nos{2},word_onsets{2},-0.05,0.5);
accr_word_erp = SP_epoch_multi(control_ERP,word_trial_nos{1},word_onsets{1},-0.05,0.5);
wrong_word_erp = SP_epoch_multi(control_ERP,word_trial_nos{2},word_onsets{2},-0.05,0.5);

%% find syllable onsets, group acc to accuracy and get TF data organized
syll_trial_nos{2}=[];syll_onsets{2}=[];
for t = 1:size(control_events,1)
    for w = 1:size(control_events.syllable_info{t},1)
        if strcmp(control_events.word_info{t}.response{control_events.syllable_info{t}.word_id(w)},'1')
            syll_trial_nos{1}(end+1,1) = t;
            syll_onsets{1}(end+1,1) = control_wlt.time(nearest(control_wlt.time, control_events.syllable_info{t}.onset(w)));
        elseif ~strcmp(control_events.word_info{t}.response{control_events.syllable_info{t}.word_id(w)},'0')
            syll_trial_nos{2}(end+1,1) = t;
            syll_onsets{2}(end+1,1) = control_wlt.time(nearest(control_wlt.time, control_events.syllable_info{t}.onset(w)));
        end
    end
end
[accr_syll_tf,freq_ITPC,time_ITPC] = SP_epoch_multi(control_wlt,syll_trial_nos{1},syll_onsets{1},-0.05,0.5);
wrong_syll_tf = SP_epoch_multi(control_wlt,syll_trial_nos{2},syll_onsets{2},-0.05,0.5);
[accr_syll_erp,~,time_ERP] = SP_epoch_multi(control_ERP,syll_trial_nos{1},syll_onsets{1},-0.05,0.5);
wrong_syll_erp = SP_epoch_multi(control_ERP,syll_trial_nos{2},syll_onsets{2},-0.05,0.5);


%% get HFA parts
[accr_word_hfa,~,time_hfa] = SP_epoch_multi(hfa_wlt,word_trial_nos{1},word_onsets{1},-0.05,0.5);
wrong_word_hfa = SP_epoch_multi(hfa_wlt,word_trial_nos{2},word_onsets{2},-0.05,0.5);
accr_syll_hfa = SP_epoch_multi(hfa_wlt,syll_trial_nos{1},syll_onsets{1},-0.05,0.5);
wrong_syll_hfa = SP_epoch_multi(hfa_wlt,syll_trial_nos{2},syll_onsets{2},-0.05,0.5);

%% Calculate POS statistics
p_zPOS{2}=[];
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
[~, ~, p_zPOS{1}] = PhaseOpposition(data1, data2, 1000, 3);


% equalize number of events(trials)
if size(accr_syll_tf,1) > size(wrong_syll_tf,1)
    nt = randperm(size(accr_syll_tf,1),size(wrong_syll_tf,1));
    accr_syll_tf = accr_syll_tf(nt,:,:,:);
else
    nt = randperm(size(wrong_syll_tf,1),size(accr_syll_tf,1));
    wrong_syll_tf = wrong_syll_tf(nt,:,:,:);
end
data1 = permute(accr_syll_tf,[2,3,4,1]);
data2 = permute(wrong_syll_tf,[2,3,4,1]);
% calculate different phase opposition values, best seems to be p_zPOS
[~, ~, p_zPOS{2}] = PhaseOpposition(data1, data2, 1000, 3);

clear data1 data2
%% Calculate ITPCz for each conditions
itpc = [];
% for accurate words
tmp      = accr_word_tf./abs(accr_word_tf);   % divide by amplitude
tmp      = sum(tmp,1);                            % sum angles across trials
tmp      = abs(tmp)/size(accr_word_tf,1);       % take the absolute value and normalize
% itpc(1,1,:,:,:) = squeeze(tmp); % remove the first singleton dimension
itpc(1,1,:,:,:) = size(accr_word_tf,1) * squeeze(tmp).^2; % remove the first singleton dimension and n*ITPC^2
% for wrong words
tmp      = wrong_word_tf./abs(wrong_word_tf);   % divide by amplitude
tmp      = sum(tmp,1);                            % sum angles across trials
tmp      = abs(tmp)/size(wrong_word_tf,1);       % take the absolute value and normalize
% itpc(1,2,:,:,:) = squeeze(tmp); % remove the first singleton dimension
itpc(1,2,:,:,:) = size(wrong_word_tf,1) * squeeze(tmp).^2; % remove the first singleton dimension and n*ITPC^2

% for accurate sylls
tmp      = accr_syll_tf./abs(accr_syll_tf);   % divide by amplitude
tmp      = sum(tmp,1);                            % sum angles across trials
tmp      = abs(tmp)/size(accr_syll_tf,1);       % take the absolute value and normalize
% itpc(2,1,:,:,:) = squeeze(tmp); % remove the first singleton dimension
itpc(2,1,:,:,:) = size(accr_syll_tf,1) * squeeze(tmp).^2; % remove the first singleton dimension and n*ITPC^2
% for wrong sylls
tmp      = wrong_syll_tf./abs(wrong_syll_tf);   % divide by amplitude
tmp      = sum(tmp,1);                            % sum angles across trials
tmp      = abs(tmp)/size(wrong_syll_tf,1);       % take the absolute value and normalize
% itpc(2,2,:,:,:) = squeeze(tmp); % remove the first singleton dimension
itpc(2,2,:,:,:) = size(wrong_syll_tf,1) * squeeze(tmp).^2; % remove the first singleton dimension and n*ITPC^2

clear tmp
%% The mini greatest plot
save_folder = fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'), 'WrdSyll_CorrWrng');
if ~exist(save_folder,'dir'),mkdir(save_folder),end
bwr = load('bwr_cmap.mat');

plh = 2; % plot height: word onset locked, syll onset locked
plw = 5; % plot width: ERP, HFA, ITPC-correct, ITPC-wrongresp, POS

for el = 1:length(control_wlt.label)
    curr_label = control_wlt.label{el};
    
    figure('Units','normalized','Position', [0 0  1 .8]);
    
    % plot average of corr&wrong - ERP - word
    subplot(plh,plw,subplotno(plw,1,1))
    for_avg=squeeze(accr_word_erp(:,el,:));
    hold on
    shadedErrorBar(time_ERP,mean(for_avg,1),stderr(for_avg),'lineprops','b');
    for_avg=squeeze(wrong_word_erp(:,el,:));
    shadedErrorBar(time_ERP,mean(for_avg,1),stderr(for_avg),'lineprops','r');
    title({'ERP (mean+stderr)';[num2str(size(accr_word_erp,1)) ' correct response - ' num2str(size(wrong_word_erp,1)) ' wrong response']})
    xlim([time_ERP(1) time_ERP(end)])
    set(gca, 'FontSize',12);
    plot([0 0], ylim,'k')
    plot(xlim,[0 0],'k')
    legend({'Correct','wrong resp'})
    
    % plot average of corr&wrong - ERP - syll
    subplot(plh,plw,subplotno(plw,2,1))
    for_avg=squeeze(accr_syll_erp(:,el,:));
    hold on
    shadedErrorBar(time_ERP,mean(for_avg,1),stderr(for_avg),'lineprops','b');
    for_avg=squeeze(wrong_syll_erp(:,el,:));
    shadedErrorBar(time_ERP,mean(for_avg,1),stderr(for_avg),'lineprops','r');
    title({'';[num2str(size(accr_syll_erp,1)) ' correct response - ' num2str(size(wrong_syll_erp,1)) ' wrong response']})
    xlim([time_ERP(1) time_ERP(end)])
    set(gca, 'FontSize',12);
    plot([0 0], ylim,'k')
    plot(xlim,[0 0],'k')
    legend({'Correct','wrong resp'})
    
    % plot average of corr&wrong - HFA - word
    subplot(plh,plw,subplotno(plw,1,2))
    for_avg=squeeze(accr_word_hfa(:,el,:));
    hold on
    shadedErrorBar(time_hfa,mean(for_avg,1),stderr(for_avg),'lineprops','b');
    for_avg=squeeze(wrong_word_hfa(:,el,:));
    shadedErrorBar(time_hfa,mean(for_avg,1),stderr(for_avg),'lineprops','r');
    title({'HFA (mean+stderr)';[num2str(size(accr_word_hfa,1)) ' correct response - ' num2str(size(wrong_word_hfa,1)) ' wrong response']})
    xlim([time_hfa(1) time_hfa(end)])
    ylims=[0 5];
    set(gca, 'FontSize',12);
    plot([0 0], ylims,'k')
    plot(xlim,[1 1],'k')
    legend({'Correct','wrong resp'})
    ylim(ylims)
    
    % plot average of corr&wrong - HFA - syll
    subplot(plh,plw,subplotno(plw,2,2))
    for_avg=squeeze(accr_syll_hfa(:,el,:));
    hold on
    shadedErrorBar(time_hfa,mean(for_avg,1),stderr(for_avg),'lineprops','b');
    for_avg=squeeze(wrong_syll_hfa(:,el,:));
    shadedErrorBar(time_hfa,mean(for_avg,1),stderr(for_avg),'lineprops','r');
    title({'';[num2str(size(accr_syll_hfa,1)) ' correct response - ' num2str(size(wrong_syll_hfa,1)) ' wrong response']})
    xlim([time_hfa(1) time_hfa(end)])
    set(gca, 'FontSize',12);
    plot([0 0], ylims,'k')
    plot(xlim,[1 1],'k')
    legend({'Correct','wrong resp'})
    ylim(ylims)
    
    % Loop to create the ITPC and POs plots
    for pp = 1:2 % first word onset, second syll onset
        for cond = 1:2 % loop columns: first correct, second wrong response
            
            subplot(plh, plw, subplotno(plw,pp,2+cond))
            surf(time_ITPC, freq_ITPC, squeeze(itpc(pp,cond,el,:,:)));set(gcf,'renderer','zbuffer');
            axis xy
            set(gca, 'FontSize',13,'FontWeight','bold');
            caxis([0 5])
            view(0,90); axis tight;
            
            hold on
            shading interp;
            plot3([0 0],ylim,[15 15],'k');
            
            if cond == 1 && pp == 1
                title({'ITPCz - Correct response events';['Random ' num2str(size(accr_word_tf,1)) ' events']});
            elseif cond == 2 && pp == 2
                title({'';['Random ' num2str(size(wrong_syll_tf,1)) ' events']});
            elseif cond == 2 && pp == 1
                title({'ITPCz - Wrong response events';['Random ' num2str(size(wrong_word_tf,1)) ' events']});
            elseif cond == 1 && pp == 2
                title({'';['Random ' num2str(size(accr_syll_tf,1)) ' events']});
            end
            
            if cond==1 && pp == 1
                ylabel({'Word onset locked';'Frequency (Hz)'});
            end
            if cond==1 && pp == 2
                ylabel({'Syll onset locked locked';'Frequency (Hz)'});
            end
            if pp==2
                xlabel('Time (s)');
            end
            
        end
        
        % plot raw p-values
        subplot(plh,plw,subplotno(plw,pp,5))
        pcolor(time_ITPC, freq_ITPC, -log10(squeeze(p_zPOS{pp}(el,:,:))));
        hold on;
        shading interp;
        contourf(time_ITPC, freq_ITPC, squeeze(p_zPOS{pp}(el,:,:)),[0.05 0.05],'Fill','off');
        set(gca, 'FontSize',13,'FontWeight','bold');
        if pp==1
            title({'Phase Opposition Sum';'p-values of z-scored POS'})
        else
            title('p-values of z-scored POS')
        end
        clim_p = [1 3];
        set(gca,'CLim',clim_p)  %for colorbar;
        plot3([0 0],ylim,[15 15],'k');
        if pp==1
            xlabel('Time (s)');
        end
        
    end
    
    % Create and delete new axes to plot colorbar of ITPC
    ax = axes;
    colormap(bwr.rgb_vals);
    cmaph = colorbar(ax);
    cmaph.Ticks = linspace(0,1,6);
    cmaph.TickLabels = num2cell(linspace(0,5,6));
    cmaph.FontSize = 13;cmaph.FontWeight='bold';
    cmaph.LineWidth = 1;
    colorTitleHandle = get(cmaph,'Title');
    set(colorTitleHandle ,'String','ITPCz values','FontSize',13,'FontWeight','bold','Position',[295 -35 0]);
    a=get(cmaph); %gets properties of colorbar
    a = a.Position; %gets the positon and size of the color bar
    set(cmaph,'Location','southoutside') % to change orientation
    set(cmaph,'Position',[a(1)/2+0.02 0.04 0.28 0.02]) % To change size
    ax.Visible = 'off';
    
    % Create and delete new axes to plot colorbar of p-values
    cmaph2 = colorbar(ax);
    xticksOI = [0.1 0.05 0.001];
    cmaph2.Ticks =  -(xticksOI-max(xticksOI))/range(xticksOI); %-(x-max_val)/ran : to fit it between zero and 1
    cmaph2.TickLabels = strsplit(num2str(xticksOI),' ');
    cmaph2.FontSize = 13;cmaph2.FontWeight='bold';
    cmaph2.LineWidth = 1;
    colorTitleHandle = get(cmaph2,'Title');
    set(colorTitleHandle ,'String','p-values','FontSize',13,'FontWeight','bold','Position',[145 -35 0]);
    a2=get(cmaph2); %gets properties of colorbar
    a2 = a2.Position; %gets the positon and size of the color bar
    set(cmaph2,'Location','southoutside') % to change orientation
    set(cmaph2,'Position',[8*a2(1)/9-0.01 0.04 0.13 0.02]) % To change size
    ax.Visible = 'off';
    
    % add the row names
    text(-0.07,0.6,'Word onset locked','Units','normalized','Rotation',90,'FontSize',18,'FontWeight','bold')
    text(-0.07,0.03,'Syllable onset locked','Units','normalized','Rotation',90,'FontSize',18,'FontWeight','bold')
    
    sgtitle(['Elec: ' curr_label ' - ERP-HFA-ITPC-Phase Opposition Sum plots'], 'FontSize',16,'FontWeight','bold')
    
    % Save the figure
    fprintf('\t-Saving electrode %s, out of %d\n',curr_label,size(channel_OI,1))
    print(fullfile(save_folder,[curr_label , '_WrdSyll_CorrWrng.jpg']),'-djpeg','-r300')
    close all
    
end
end
