function SP_noiseInvariance(sbj_block,elecsOI)
% Idea: looking at the electrodes that respond solely to speech vs sounds

%% some pre-adjustments
data_root = '/media/sakkol/HDD1/HBML/';
project_name = 'Speech_Perception';
if ~exist('elecsOI','var') || isempty(elecsOI)
    elecsOI = 'allElecs';
end
if ~strcmp(elecsOI,'allElecs')
    ElecLoc=readtable(fullfile(data_root,'DERIVATIVES','freesurfer','ElecLoc_master.xlsx'));
    area_OI = {'HG','STG'};
end
% steps = 1; % x10ms (Because it is sampled in 100Hz, 1 point is 10ms)
max_delay = 40; % x10ms = 400ms

% get the speech envelope of the attention sentence
[y,Fs] = audioread('Pre-stim-Attention-comma-M.wav');
[amp_env,~,~,~] = get_speech_2features(y,Fs,[]);
% resample that to 100 Hz
amp_env = resample(amp_env,100,Fs);

%% running the correlations
all_corr_res = cell(size(sbj_block,1),1);
all_pvals = cell(size(sbj_block,1),1);
all_corrected_pval = cell(size(sbj_block,1),1);
all_collapsed = cell(size(sbj_block,1),1);
sbjs_elecs = cell(size(sbj_block,1),2);
all_collapsed_comb=[];
for s = 1:size(sbj_block,1)
    
    sbj_ID = sbj_block{s,1};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    curr_block = sbj_block{s,2};
    
    tmp=load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_info.mat']));info=tmp.info;clear tmp
    epoched_hfa = SP_getAttHFA(Sbj_Metadata, curr_block);
    
    % cross-correlate the envelope with the HFA
    % get both max correlation and correlation lag
    corr_res = zeros(length(info.channelinfo.Label),max_delay); % elec x delay
    pval = zeros(length(info.channelinfo.Label),max_delay); % elec x delay
    for el = 1:size(corr_res,1)          % loop electrode
        for d = 1:max_delay                 % loop delay
            % cut beginning of HFA and end of envelope
            hfa = squeeze(mean(epoched_hfa.powspctrm(:,el,:,101+d+1:end),1));
            hfa = hfa(1:length(amp_env));hfa(isnan(hfa))=0;
            [corr_res(el,d),pval(el,d)] = corr(hfa,amp_env,'Type','Spearman');
            if pval(el,d)>0.05/numel(pval)
                corr_res(el,d)=0;
            end
        end
    end
    % correct for family-wise error
    corrected_pval=bonf_holm(pval);
    for el = 1:size(corr_res,1)          % loop electrode
        for d = 1:max_delay                 % loop delay
            if corrected_pval>0.05
                corr_res(el,d)=0;
            end
        end
    end
    
    % remove bad channels
    [good_chans,good_chans_idx] = get_info_goodchans(info);
    % get the elecsOI
    if ~strcmp(elecsOI,'allElecs')
        elecarea = ElecLoc(strcmp(ElecLoc.Subject,sbj_ID),:);
        elecs_OI = elecarea(ismember(elecarea.Area,area_OI),:);
        good_chans_idx = good_chans_idx & ismember(info.channelinfo.Label,elecs_OI.Label);
        good_chans = good_chans(ismember(good_chans,elecs_OI.Label));
    end
    
    % find the max in each correlation delay
    collapsed=zeros(length(good_chans_idx),2);
    [M,I] = max(corr_res,[],2);
    collapsed(:,1) = M;
    collapsed(:,2) = I*10;
    % select the necessary channels
    good_chans_idx = good_chans_idx & collapsed(:,1)>0.1;
    sbjs_elecs{s,1} = Sbj_Metadata.fsname;
    sbjs_elecs{s,2} = info.channelinfo.Label(good_chans_idx);
    all_corr_res{s,1} = corr_res(good_chans_idx,:);
    all_pvals{s,1} = pval(good_chans_idx,:);
    all_corrected_pval{s,1} = corrected_pval(good_chans_idx,:);
    all_collapsed{s,1} = collapsed(good_chans_idx,:);
    all_collapsed_comb = [all_collapsed_comb;all_collapsed{s,1}];
    
end


%% Now the plotting
% set the names according to number of subjects
if size(sbjs_elecs,1)==1
    coordName='PIAL';
    save_folder = fullfile(Sbj_Metadata.results,'NoiseInvariance');
else
    coordName='FSAVERAGE';
    save_folder = fullfile(Sbj_Metadata.project_root,'Collected_Results','NoiseInvariance');
    Sbj_Metadata = 'fsaverage';
end
if ~exist(save_folder,'dir'),mkdir(save_folder);end
% gather electrodes
[AllSubElecNames,AllSubElecCoords] = gather_allelecsinfo(sbjs_elecs,coordName);

if ~(sum(all_collapsed_comb(:,1)==0)==length(all_collapsed_comb(:,1))) % this is to skip those cases when there are no significant following
    
%     % just to remove outliers for once:
%     rm_idx = all_collapsed_comb(:,2)<250;
%     all_collapsed_comb = all_collapsed_comb(rm_idx,:);
%     AllSubElecNames = AllSubElecNames(rm_idx,:);
%     AllSubElecCoords = AllSubElecCoords(rm_idx,:);


% first plot the correlation values
color_elecs_wData(Sbj_Metadata,all_collapsed_comb(:,1),AllSubElecNames,AllSubElecCoords,'continuous',[],'r-values','inferno')
text(gca,.5,1.07,'Max correlation values','Units','normalized','FontSize',18,'FontWeight','bold','HorizontalAlignment','center')
print(fullfile(save_folder,[strjoin(sbjs_elecs(:,1),'_') '_' elecsOI '_rvals_outlierRM.jpg']),'-djpeg','-r300')

% second plot the correlation lag
color_elecs_wData(Sbj_Metadata,all_collapsed_comb(:,2),AllSubElecNames,AllSubElecCoords,'continuous',[],'r-values','inferno')
text(gca,.5,1.07,'Max correlation lag','Units','normalized','FontSize',18,'FontWeight','bold','HorizontalAlignment','center')
print(fullfile(save_folder,[strjoin(sbjs_elecs(:,1),'_') '_' elecsOI '_ROL_outlierRM.jpg']),'-djpeg','-r300')
close all

% lastly: correlate the delay and lag
figure('Position',[0 0 1400 1000])
[ccoorr_c,ccoorr_r] = corr(all_collapsed_comb(:,1),all_collapsed_comb(:,2),'Type','Spearman');
plot(all_collapsed_comb(:,1),all_collapsed_comb(:,2),'o')
hold on
fit_res = fit(all_collapsed_comb(:,1),all_collapsed_comb(:,2), 'poly1');
h=plot(fit_res);h.Color = [0 0 0];
h.LineWidth = .75;legend('hide');
xlabel('Correlation');ylabel('Delay')
set(gca,'FontSize',14)
title(['Correlating delay with max HFA-envelope correlation across electrodes - rho=' num2str(ccoorr_c) ' (p=' num2str(ccoorr_r) ')']);
print(fullfile(save_folder,[strjoin(sbjs_elecs(:,1),'_') '_' elecsOI '_ROL_rval_corr_outlierRM.jpg']),'-djpeg','-r300')
close all
end

% save the output
save(fullfile(save_folder,[strjoin(sbjs_elecs(:,1),'_') '_' elecsOI '_output.mat']),'all_collapsed_comb','all_corr_res','all_pvals','sbjs_elecs')


end


