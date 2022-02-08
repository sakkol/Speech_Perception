function SP_noiseInvariance(sbj_block)

% some pre-adjustments
data_root = '/media/sakkol/HDD1/HBML/';
project_name = 'Speech_Perception';
freq = 70:5:200;freq(ismember(freq,[120 180]))=[];
fourier_pow = 'pow';
steps = 1; % x10ms (Because it is sampled in 100Hz, 1 point is 10ms)
max_delay = 40; % x10ms = 400ms

% get the speech envelope of the attention sentence
[y,Fs] = audioread('Pre-stim-Attention-comma-M.wav');
[amp_env,~,~,~] = get_speech_2features(y,Fs,[]);
% resample that to 100 Hz
amp_env = resample(amp_env,100,Fs);

all_corr_res = cell(size(sbj_block,1),1);
all_collapsed = cell(size(sbj_block,1),1);
sbjs_elecs = cell(size(sbj_block,1),1);
all_collapsed_comb=[];
for s = 1:size(sbj_block,1)
    
    sbj_ID = sbj_block{s,1};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    curr_block = sbj_block{s,2};
    
    load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_info.mat']))
    load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_avg.mat']))
    
    % run wavelet
    pre  = 1;
    post = 3.5;
    [~, epoched_wlt] = master_waveletTF(ecog_avg.ftrip,info.events.att_sent_onset,pre,post,freq,fourier_pow);
    
    % run HFA
    cfg=[];
    cfg.frequency=[69 201];
    cfg.avgoverfreq='no';
    epoched_hfa = ft_selectdata(cfg, epoched_wlt);
    cfg              = [];
    cfg.baseline     = [-.95 -.55];
    cfg.baselinetype = 'relative';
    cfg.parameter    = 'powspctrm';
    [epoched_hfa]    = ft_freqbaseline(cfg, epoched_hfa);
    cfg=[];
    cfg.frequency=[70 200];
    cfg.avgoverfreq='yes';
    epoched_hfa = ft_selectdata(cfg, epoched_hfa);
%     epoched_hfa = master_smoothData(epoched_hfa);
    
    % cross-correlate the envelope with the HFA
    % get both max correlation and correlation lag
    corr_res = zeros(length(info.channelinfo.Label),max_delay); % elec x delay
    pval = zeros(length(info.channelinfo.Label),max_delay); % elec x delay
    for el = 1:size(corr_res,1)          % loop electrode
        fprintf('running at %dth of %d electrodes\n',el,size(corr_res,1))
        for d = 1:max_delay                 % loop delay
            % cut beginning of HFA and end of envelope
            hfa = squeeze(mean(epoched_hfa.powspctrm(:,el,:,101+d+1:end),1));
            hfa = hfa(1:length(amp_env));hfa(isnan(hfa))=0;
            [corr_res(el,d),pval(el,d)] = corr(hfa,amp_env,'Type','Spearman');
            if pval(el,d)>0.05
                corr_res(el,d)=0;
            end
        end
    end
    
    % remove bad channels
    [good_chans,good_chans_idx] = get_info_goodchans(info);
    sbjs_elecs{s,1} = Sbj_Metadata.fsname;
    sbjs_elecs{s,2} = good_chans;
    
    all_corr_res{s,1} = corr_res(good_chans_idx,:);
    
    % find the collapsed results
    collapsed=zeros(sum(good_chans_idx),2);
    [M,I] = max(all_corr_res{s,1},[],2);
    collapsed(:,1) = M;
    collapsed(:,2) = I*10;
    
    all_collapsed{s,1} = collapsed;
    all_collapsed_comb = [all_collapsed_comb;collapsed];
    
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

% first plot the correlation values
figure('Units','normalized','Position', [0 0  1 1]);
color_elecs_wData(Sbj_Metadata,all_collapsed_comb(:,1),AllSubElecNames,AllSubElecCoords,'continuous',[],'r-values','inferno')
text(gca,.5,1.07,'Max correlation','Units','normalized','FontSize',18,'FontWeight','bold','HorizontalAlignment','center')
print(fullfile(save_folder,[strjoin(sbjs_elecs(:,1),'_') '_rvals.jpg']),'-djpeg','-r300')

% second plot the correlation lag
figure('Units','normalized','Position', [0 0  1 1]);
color_elecs_wData(Sbj_Metadata,all_collapsed_comb(:,2),AllSubElecNames,AllSubElecCoords,'continuous',[],'r-values','inferno')
text(gca,.5,1.07,'Max correlation lag','Units','normalized','FontSize',18,'FontWeight','bold','HorizontalAlignment','center')
print(fullfile(save_folder,[strjoin(sbjs_elecs(:,1),'_') '_ROL.jpg']),'-djpeg','-r300')


end


