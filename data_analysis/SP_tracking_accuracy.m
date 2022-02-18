function SP_tracking_accuracy(sbj_block,elecsOI)
% Idea: if neural tracking is higher during the attention sentence,
% accuracy of the target word detection would be higher.

%% some pre-adjustments
data_root = '/media/sakkol/HDD1/HBML/';
project_name = 'Speech_Perception';
if ~exist('elecsOI','var') || isempty(elecsOI)
    elecsOI = 'allElecs';
end
% steps = 1; % x10ms (Because it is sampled in 100Hz, 1 point is 10ms)
max_delay = 40; % x10ms = 400ms

% get the speech envelope of the attention sentence
[y,Fs] = audioread('Pre-stim-Attention-comma-M.wav');
[amp_env,~,~,~] = get_speech_2features(y,Fs,[]);
% resample that to 100 Hz
amp_env = resample(amp_env,100,Fs);
sbjs_elecs = cell(size(sbj_block,1),2);
best_delay = cell(size(sbj_block,1),2);
all_corr_res = cell(size(sbj_block,1),2);
all_pvals = cell(size(sbj_block,1),2);

%% run the correlations
for s = 1:size(sbj_block,1)
    
    sbj_ID = sbj_block{s,1};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    curr_block = sbj_block{s,2};
    fprintf('running at %dth of %d subjects\n',s,size(sbj_block,1))
    
    tmp=load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_info.mat']));info=tmp.info;clear tmp
    epoched_hfa = SP_getAttHFA(Sbj_Metadata, curr_block);
    
    % cross-correlate the envelope with the HFA of each trial
    % get both max correlation and correlation lag
    corr_res = zeros(length(info.channelinfo.Label),size(epoched_hfa.powspctrm,1),max_delay); % elec x trial x delay
    pval = zeros(length(info.channelinfo.Label),size(epoched_hfa.powspctrm,1),max_delay); % elec x trial x delay
    best_delay{s} = zeros(length(info.channelinfo.Label),2); % elec x 2 (best correlation and its delay)
    for el = 1:size(corr_res,1)          % loop electrode
        for d = 1:max_delay                 % loop delay
            % cut beginning of HFA and end of envelope
            hfa = squeeze(epoched_hfa.powspctrm(:,el,:,101+d+1:end));
            hfa = hfa(:,1:length(amp_env));hfa(isnan(hfa))=0;
            [corr_res(el,:,d),pval(el,:,d)] = corr(hfa',amp_env,'Type','Spearman');
        end
        
        % find the best overall delay
        [M,I] = max(squeeze(mean(corr_res(el,:,:),2)));
        best_delay{s}(el,1:2) = [M,I*10];
    end
    
    % remove bad channels
    [good_chans,good_chans_idx] = get_info_goodchans(info);
    % get the elecsOI
    if ~strcmp(elecsOI,'allElecs')
        ElecLoc=readtable(fullfile(data_root,'DERIVATIVES','freesurfer','ElecLoc_master.xlsx'));
        elecarea = ElecLoc(strcmp(ElecLoc.Subject,sbj_ID),:);
        elecs_OI = elecarea(ismember(elecarea.Area,elecsOI),:);
        good_chans_idx = good_chans_idx & ismember(info.channelinfo.Label,elecs_OI.Label);
        good_chans = good_chans(ismember(good_chans,elecs_OI.Label));
    end
        
    sbjs_elecs{s,1} = Sbj_Metadata.fsname;
    sbjs_elecs{s,2} = good_chans;
    
    all_corr_res{s,1} = corr_res(good_chans_idx,:,:);
    all_pvals{s,1} = pval(good_chans_idx,:,:);
    best_delay{s,1} = best_delay{s}(good_chans_idx,:);
    
end

%% Now the plotting: first with individual electrodes
entr_accr_corr = cell(size(sbj_block,1),2);all_entr_accr_corr=[];
for s = 1:size(sbj_block,1)
    fprintf('running at %dth of %d subjects\n',s,size(sbj_block,1))
    sbj_ID = sbj_block{s,1};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    curr_block = sbj_block{s,2};
    
    tmp=load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_info.mat']));info=tmp.info;clear tmp
    save_folder = fullfile(Sbj_Metadata.results,'Tracking_Accuracy_corr');
    if ~exist(save_folder,'dir'),mkdir(save_folder);end
    
    conds = unique(info.events.Cond_code);
    entr_accr_corr{s,1} = zeros(size(all_corr_res{s,1},1),length(conds),2);
    for el = 1:size(all_corr_res{s,1},1)
        bd = best_delay{s,1}(el,2)/10;
        figure('Units','normalized','Position', [0 0  1 .5]);
        % first plot overall
        subplot(1,length(conds)+1,1)
        yvals = all_corr_res{s,1}(el,:,bd)';
        xvals = info.events.Acc_word_count;
        [ccoorr_c,ccoorr_r] = corr(xvals,yvals,'Type','Spearman');
        
        plot(xvals,yvals,'o')
        hold on
        fit_res = fit(xvals,yvals, 'poly1');
        h=plot(fit_res);h.Color = [0 0 0];
        h.LineWidth = .75;legend('hide');
        xlim([0 6]);xticks(1:5)
        xlabel('Accuracy');ylabel('Entrainment strength')
        set(gca,'FontSize',14)
        title(['All conditions combined - rho=' num2str(ccoorr_c,2) ' (p=' num2str(ccoorr_r,2) ')']);
        
        % now plot each condition separately
        for c = 1:length(conds)
            subplot(1,length(conds)+1,c+1)
            yvals = all_corr_res{s,1}(el,info.events.Cond_code==conds(c),bd)';
            xvals = info.events.Acc_word_count(info.events.Cond_code==conds(c));
            [entr_accr_corr{s,1}(el,c,1),entr_accr_corr{s,1}(el,c,2)] = corr(xvals,yvals,'Type','Spearman');
            
            plot(xvals,yvals,'o')
            hold on
            fit_res = fit(xvals,yvals, 'poly1');
            h=plot(fit_res);h.Color = [0 0 0];
            h.LineWidth = .75;legend('hide');
            xlim([0 6]);xticks(1:5)
            xlabel('Accuracy');ylabel('Entrainment strength')
            set(gca,'FontSize',14)
            tmp=info.events.Condition(info.events.Cond_code==conds(c));
            title([tmp{1} '- rho=' num2str(entr_accr_corr{s,1}(el,c,1),2) ' (p=' num2str(entr_accr_corr{s,1}(el,c,2),2) ')']);clear tmp
        end
        all_entr_accr_corr = [all_entr_accr_corr;squeeze(entr_accr_corr{s,1}(el,:,1))];
        ax = axes;
        text(ax,.5,1.07,[sbjs_elecs{s,2}{el} ' - Correlation of accuracy with entrainment strength'],'Units','normalized','FontSize',18,'FontWeight','bold','HorizontalAlignment','center')
        ax.Visible = 'off';
        print(fullfile(save_folder,[sbjs_elecs{s,2}{el} '_plots.jpg']),'-djpeg')
        close all
    end
end


%% Now plotting the r-values of each electrode
% set the names according to number of subjects
if size(sbjs_elecs,1)==1
    coordName='PIAL';
    save_folder = fullfile(Sbj_Metadata.results,'Tracking_Accuracy_corr');
else
    coordName='FSAVERAGE';
    save_folder = fullfile(Sbj_Metadata.project_root,'Collected_Results','Tracking_Accuracy_corr');
    Sbj_Metadata = 'fsaverage';
end
if ~exist(save_folder,'dir'),mkdir(save_folder);end
if iscell(elecsOI),elecsOI = strjoin(elecsOI,'_');end
% gather electrodes
[AllSubElecNames,AllSubElecCoords] = gather_allelecsinfo(sbjs_elecs,coordName);

% plot the correlation values
color_elecs_wData(Sbj_Metadata,mean(all_entr_accr_corr,2),AllSubElecNames,AllSubElecCoords,'continuous',[],'r-values','inferno')
text(gca,.5,1.07,['Correlation of accuracy and max HFA-envelope correlation in all conditions'],'Units','normalized','FontSize',18,'FontWeight','bold','HorizontalAlignment','center')
print(fullfile(save_folder,[strjoin(sbjs_elecs(:,1),'_') '_AllConds_' elecsOI '.jpg']),'-djpeg','-r300')
close all

for c=1:length(conds)
    color_elecs_wData(Sbj_Metadata,all_entr_accr_corr(:,c),AllSubElecNames,AllSubElecCoords,'continuous',[],'r-values','inferno')
    tmp=info.events.Condition(info.events.Cond_code==conds(c));
    text(gca,.5,1.07,['Correlation of accuracy and max HFA-envelope correlation in ' tmp{1}],'Units','normalized','FontSize',18,'FontWeight','bold','HorizontalAlignment','center')
    print(fullfile(save_folder,[strjoin(sbjs_elecs(:,1),'_') '_' erase(tmp{1},' ') '_' elecsOI '.jpg']),'-djpeg','-r300')
end
close all

save(fullfile(save_folder,[strjoin(sbjs_elecs(:,1),'_') '_' elecsOI '_output.mat']),'entr_accr_corr','sbjs_elecs','best_delay','all_corr_res','all_pvals')

end