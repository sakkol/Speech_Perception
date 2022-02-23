function SP_TrackingTarget_Accuracy(sbj_block,elecsOI)
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

% resample that to 100 Hz
sbjs_elecs = cell(size(sbj_block,1),2);
all_corr_res = cell(size(sbj_block,1),1);
all_pvals = cell(size(sbj_block,1),1);

%% run the correlations
for s = 1:size(sbj_block,1)
    
    sbj_ID = sbj_block{s,1};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    curr_block = sbj_block{s,2};
    fprintf('running at %dth of %d subjects\n',s,size(sbj_block,1))
    
    tmp=load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_info.mat']));info=tmp.info;clear tmp
    [epoched_hfa,~,events_sel] = SP_getTargetHFA(Sbj_Metadata, curr_block);
            
    % get the best delay from previous analysis for the electrodes that the
    % analysis was run on
    save_folder = fullfile(Sbj_Metadata.results,'NoiseInvariance');
    tmp = load(fullfile(Sbj_Metadata.results,'NoiseInvariance',[replace(sbj_ID,'_2','_02') '_allElecs_output.mat']));
    
    % cross-correlate the envelope with the HFA of each trial
    % get both max correlation and correlation lag
    corr_res = zeros(size(tmp.all_collapsed_comb,1),size(epoched_hfa.powspctrm,1),max_delay); % elec x trial x delay
    pval = zeros(size(tmp.all_collapsed_comb,1),size(epoched_hfa.powspctrm,1),max_delay); % elec x trial x delay
    if strcmp(events_sel.lang{1},'Spanish')
        audio_dir = '/home/sakkol/Documents/TASKS/Matrix_Speech_Task/Spanish_main_stim_loc';
    else
        audio_dir = '/home/sakkol/Documents/TASKS/Matrix_Speech_Task/English_main_stim_loc';
    end
    
    sbjs_elecs{s,1} = Sbj_Metadata.fsname;
    sbjs_elecs{s,2} = tmp.sbjs_elecs{1,2};
    select_chans_idx = ones(size(tmp.all_collapsed_comb,1),1);
    % get the elecsOI
    if ~strcmp(elecsOI,'allElecs')
        ElecLoc=readtable(fullfile(data_root,'DERIVATIVES','freesurfer','ElecLoc_master.xlsx'));
        elecarea = ElecLoc(strcmp(ElecLoc.Subject,sbj_ID),:);
        elecs_OI = elecarea(ismember(elecarea.Area,elecsOI),:);
        select_chans_idx = ismember(sbjs_elecs{s,2},elecs_OI.Label);
        sbjs_elecs{s,2} = sbjs_elecs{s,2}(select_chans_idx);
    end
    
    for ev = 1:size(events_sel,1)
        % get the speech envelope of the attention sentence
        filename = find_sentence(events_sel.sentence{ev},audio_dir,'0.9');
        [y,Fs] = audioread(filename);
        [amp_env,~,~,~] = get_speech_2features(y,Fs,[]);
        amp_env = resample(amp_env,100,Fs);
        for el = 1:size(corr_res,1)          % loop electrode
            el_idx = strcmp(epoched_hfa.label,tmp.sbjs_elecs{1,2}{el});
            for d = 1:max_delay                 % loop delay
                % cut beginning of HFA and end of envelope
                hfa = squeeze(epoched_hfa.powspctrm(ev,el_idx,:,101+d+1:end));
                hfa = hfa(1:length(amp_env));hfa(isnan(hfa))=0;
                [corr_res(el,ev,d),pval(el,ev,d)] = corr(hfa,amp_env,'Type','Spearman');
            end
        end
    end
    
    all_corr_res{s,1} = corr_res(select_chans_idx,:,:);
    all_pvals{s,1} = pval(select_chans_idx,:,:);
    
end

%% Now the plotting: first with individual electrodes
    
all_entr_accr_corr=[];all_entr_accr_corr_pval=[];
all_bestdelays=[];
for s = 1:size(sbj_block,1)
    fprintf('running at %dth of %d subjects\n',s,size(sbj_block,1))
    sbj_ID = sbj_block{s,1};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    curr_block = sbj_block{s,2};
    
    tmp=load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_info.mat']));info=tmp.info;clear tmp
    [~,~,events_sel] = SP_getTargetHFA(Sbj_Metadata, curr_block);
    save_folder = fullfile(Sbj_Metadata.results,'TrackingTarget_Accuracy');
    if ~exist(save_folder,'dir'),mkdir(save_folder);end
    
    if ~isempty(all_corr_res{s,1})
    
    % correlate the accuracy with entrainment correlation for each electrode
    for el = 1:size(all_corr_res{s,1},1)
        % get the best delay from previous analysis:
        tmp = load(fullfile(Sbj_Metadata.results,'NoiseInvariance',[replace(sbj_ID,'_2','_02') '_allElecs_output.mat']));
        bd = tmp.all_collapsed_comb(el,2)/10;
        all_bestdelays = [all_bestdelays;bd*10];
        
        figure('Units','normalized','Position', [0 0  .4 .5]);
        yvals = all_corr_res{s,1}(el,:,bd)';
        xvals = events_sel.Acc_word_count;
        [ccoorr_c,ccoorr_r] = corr(xvals,yvals,'Type','Spearman');
        
        plot(xvals,yvals,'o')
        hold on
        fit_res = fit(xvals,yvals, 'poly1');
        h=plot(fit_res);h.Color = [0 0 0];
        h.LineWidth = .75;legend('hide');
        xlim([0 6]);xticks(1:5)
        xlabel('Accuracy');ylabel('Entrainment strength')
        set(gca,'FontSize',14)
        title(['rho=' num2str(ccoorr_c,2) ' (p=' num2str(ccoorr_r,2) ')']);
        
        all_entr_accr_corr = [all_entr_accr_corr;ccoorr_c];
        all_entr_accr_corr_pval = [all_entr_accr_corr_pval;ccoorr_r];
        ax = axes;
        text(ax,.5,1.07,[sbjs_elecs{s,2}{el} ' - Correlation of accuracy with entrainment strength'],'Units','normalized','FontSize',18,'FontWeight','bold','HorizontalAlignment','center')
        ax.Visible = 'off';
        print(fullfile(save_folder,[sbjs_elecs{s,2}{el} '_plots.jpg']),'-djpeg')
        close all
    end
    end
end


%% Now plotting the r-values of each electrode
% set the names according to number of subjects
if size(sbjs_elecs,1)==1
    coordName='PIAL';
    save_folder = fullfile(Sbj_Metadata.results,'TrackingTarget_Accuracy');
else
    coordName='FSAVERAGE';
    save_folder = fullfile(Sbj_Metadata.project_root,'Collected_Results','TrackingTarget_Accuracy');
    Sbj_Metadata = 'fsaverage';
end
if ~exist(save_folder,'dir'),mkdir(save_folder);end
if iscell(elecsOI),elecsOI = strjoin(elecsOI,'_');end
% gather electrodes
[AllSubElecNames,AllSubElecCoords] = gather_allelecsinfo(sbjs_elecs,coordName);

if ~(length(all_entr_accr_corr) < 3)
% plot the correlation values
color_elecs_wData(Sbj_Metadata,all_entr_accr_corr,AllSubElecNames,AllSubElecCoords,'continuous',[],'r-values','inferno')
text(gca,.5,1.07,['Correlation of accuracy and max HFA-envelope correlation in all conditions'],'Units','normalized','FontSize',18,'FontWeight','bold','HorizontalAlignment','center')
print(fullfile(save_folder,[strjoin(sbjs_elecs(:,1),'_') '_' elecsOI '_brain.jpg']),'-djpeg','-r300')
close all

% plot correlation of ROL with how predictive electrode is
figure('Units','normalized','Position', [0 0  .2 .5]);
yvals = all_entr_accr_corr;
xvals = all_bestdelays;
[ccoorr_c,ccoorr_r] = corr(xvals,yvals,'Type','Spearman');

plot(xvals,yvals,'o')
hold on
fit_res = fit(xvals,yvals, 'poly1');
h=plot(fit_res);h.Color = [0 0 0];
h.LineWidth = .75;legend('hide');
xlabel('ROL');ylabel('Entrainment strength')
set(gca,'FontSize',14)
title(['rho=' num2str(ccoorr_c,2) ' (p=' num2str(ccoorr_r,2) ')']);
print(fullfile(save_folder,[strjoin(sbjs_elecs(:,1),'_') '_' elecsOI '_corrplot.jpg']),'-djpeg')

end

save(fullfile(save_folder,[strjoin(sbjs_elecs(:,1),'_') '_' elecsOI '_output.mat']),'all_entr_accr_corr','all_entr_accr_corr_pval','sbjs_elecs','all_bestdelays','all_corr_res','all_pvals')



end