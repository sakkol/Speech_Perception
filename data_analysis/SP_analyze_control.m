%% TO-DO:
% Add Rayleigh stats, because it is important only for control
control_wlt = addRayleigh_to_ft(control_wlt); % ?????????

%% Prepare
data_root = '/media/sakkol/HDD1/HBML/';
project_name = 'Speech_Perception';
sbj_ID = 'NS148';
Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'

% Select blocks to import
blocklistsheet = [char(Sbj_Metadata.project_root) filesep char(Sbj_Metadata.project_name) '_BlockInfo.xlsx']; % "F:\HBML\PROJECTS_DATA\CL_Train\CL_Train_BlockLists.xlsx";
blocklistall = readtable(blocklistsheet);
control_blocks = blocklistall.BlockList(strcmpi(blocklistall.sbj_ID,Sbj_Metadata.sbj_ID) & contains(blocklistall.conditions_code,'1'));

clear blocklistall blocklistsheet

%% bring in these blocks and combine only the control events

for b = 1:length(control_blocks)
    curr_block = control_blocks{b};
    
    % Load iEEG
    load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_wlt.mat']))
    events = epoched_wlt.events;
    
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

%% Separate fourierspectrums of correct responses and others
corr_rspn_fft = {0};
no_rspn_fft = {0};
wrng_rspn_fft = {0};

% Which freq bands
freq_bands = {[1 4], [4 8],[8 12]};


for f = 1:length(freq_bands)
    cl_freqs = [nearest(control_wlt.freq, freq_bands{f}(1)), ...
        nearest(control_wlt.freq,freq_bands{f}(2))];
    
    corr_ind = 1;
    no_ind = 1;
    wrng_ind = 1;
    
    
    for t = 1:size(control_events,1)
        
        for w = 1:5
            % Separate words
            % Check if correct
            % Collect fft in a cell structure?
            if strcmp(control_events.word_info{t}.response{w},'1')
                % find closest timepoint
                cl_times = [nearest(control_wlt.time, control_events.word_info{t}.onset(w))...
                    nearest(control_wlt.time,control_events.word_info{t}.offset(w))];
                corr_rspn_fft{corr_ind,f} = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)));
                corr_ind = corr_ind+1;
            elseif strcmp(control_events.word_info{t}.response{w},'0')
                % find closest timepoint
                cl_times = [nearest(control_wlt.time, control_events.word_info{t}.onset(w))...
                    nearest(control_wlt.time,control_events.word_info{t}.offset(w))];
                no_rspn_fft{no_ind,f} = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)));
                no_ind = no_ind+1;
                
            else
                % find closest timepoint
                cl_times = [nearest(control_wlt.time, control_events.word_info{t}.onset(w))...
                    nearest(control_wlt.time,control_events.word_info{t}.offset(w))];
                wrng_rspn_fft{wrng_ind,f} = squeeze(control_wlt.fourierspctrm(t,:,cl_freqs(1):cl_freqs(2),cl_times(1):cl_times(2)));
                wrng_ind = wrng_ind+1;
                
            end
            
        end
    end
end

fft_of_words = [];

fft_of_words.freq_bands = freq_bands;
fft_of_words.freq_band_dtls{1} = control_wlt.freq(nearest(control_wlt.freq, freq_bands{1}(1)):...
    nearest(control_wlt.freq,freq_bands{1}(2)));
fft_of_words.freq_band_dtls{2} = control_wlt.freq(nearest(control_wlt.freq, freq_bands{2}(1)):...
    nearest(control_wlt.freq,freq_bands{2}(2)));
fft_of_words.freq_band_dtls{3} = control_wlt.freq(nearest(control_wlt.freq, freq_bands{3}(1)):...
    nearest(control_wlt.freq,freq_bands{3}(2)));

fft_of_words.corr_rspn_fft = corr_rspn_fft;
fft_of_words.no_rspn_fft = no_rspn_fft;
fft_of_words.wrng_rspn_fft = wrng_rspn_fft;

save(fullfile(Sbj_Metadata.results, [strjoin(control_blocks,'_') '_ctrl_word_fft.mat']),'fft_of_words','-v7.3');

%% Compare angles
comp_angle_res{size(corr_rspn_fft{1},1),size(corr_rspn_fft,2)} = [];
for el = 1:size(corr_rspn_fft{1},1) % loop electrode
    for f = 1:size(corr_rspn_fft,2) % loop frequency band
        % average frequency and time dimension per trial
        for i=1:size(corr_rspn_fft,1)
            corr_angles(i) = squeeze(mean(angle(corr_rspn_fft{i,f}(1,:,:)),[2 3]));
        end
        
        for i=1:size(no_rspn_fft,1)
            noncorr_angles(i) = squeeze(mean(angle(no_rspn_fft{i,f}(1,:,:)),[2 3]));
        end
        for i=1:size(wrng_rspn_fft,1)
            
            noncorr_angles(end+1) = squeeze(mean(angle(wrng_rspn_fft{i,f}(1,:,:)),[2 3]));
        end
        
        % Test for difference in angles
        [pval, table] = circ_wwtest(corr_angles', noncorr_angles');
        
        comp_angle_res{el,f}.pval = pval;
        comp_angle_res{el,f}.table = table;
        
        % plot rose plot
        figure('Units','normalized','Position', [0 0  .5 .5]);
        subplot(1,2,1)
        circ_plot(corr_angles,'hist',[],90,true,true,'linewidth',2,'color','r');
        subplot(1,2,2)
        circ_plot(noncorr_angles,'hist',[],90,true,true,'linewidth',2,'color','r');
        
        if pval > 0.05
            sgtitle({['Blocks: ' strjoin(control_blocks,' & ')];...
                ['Correct vs Non-correct (wrong+unheard) word report circular stats p-value: ' num2str(pval)]})
        else
            sgtitle({['Blocks: ' strjoin(control_blocks,' & ')];...
                ['Correct vs Non-correct (wrong+unheard) word report circular stats p-value: ' num2str(pval)]},'Color','r')
        end
        
        print(fullfile(Sbj_Metadata.iEEG_data,curr_block,'PICS','events.jpg'),'-djpeg','-r300')
    end
end


% Stats:
% First normality:
[p,~] = rayleigh(angle(bl_data(chn_i,:))');
eval([plotconds_toname{conds}.BlockList{plot_i} ' = [];']);
eval([plotconds_toname{conds}.BlockList{plot_i} '.bl_rayleigh.p = p;']);

[p,~] = rayleigh(angle(post_data(chn_i,:))');
eval([plotconds_toname{conds}.BlockList{plot_i} ' = [];']);
eval([plotconds_toname{conds}.BlockList{plot_i} '.post_rayleigh.p = p;']);

[pval, table] = circ_wwtest(angle(bl_data(chn_i,:))', angle(post_data(chn_i,:))');
eval([plotconds_toname{conds}.BlockList{plot_i} '.ttest.p = pval;']);
eval([plotconds_toname{conds}.BlockList{plot_i} '.ttest.table = table;']);


% 2 roseplots for PLVs
subtightplot(3,plot_width,[2*plot_i-1, 2*plot_i],[0.1 0.01],0.1,0.02);

circ_plot(angle(bl_data(chn_i,:)),'hist',[],90,true,true,'linewidth',2,'color','r');

if pval > 0.05
    title({['Block: ' plotconds{conds}.BlockList{plot_i} '; Stim:' num2str(plotconds{conds}.trainFreq(plot_i)) 'Hz'];...
        ['Pre vs post circular stats p-value: ' num2str(pval)];'';'Prestim PLV roseplot'})
else
    title({['Block: ' plotconds{conds}.BlockList{plot_i} '; Stim:' num2str(plotconds{conds}.trainFreq(plot_i)) 'Hz'];...
        ['Pre vs post circular stats p-value: ' num2str(pval)];'';'Prestim PLV roseplot'},'Color','r')
end

subtightplot(3,plot_width,[2*plot_i+plot_width-1, 2*plot_i+plot_width],[0.1 0.01],0.1,0.02);

circ_plot(angle(post_data(chn_i,:)),'hist',[],90,true,true,'linewidth',2,'color','r');

title({'Prestim PLV roseplot'})








