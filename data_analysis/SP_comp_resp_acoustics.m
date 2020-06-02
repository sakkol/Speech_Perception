function SP_comp_resp_acoustics(Sbj_Metadata,control_blocks)
% Idea: During control condition, words that were perceived might differ
% not because of the underlying phase (alone), but the acoustics of the
% word may play a role. For example, peak of the amplitude may be pretty
% high that it exceeds the noise. Then, patient may fill the underperceived
% syllables from memory/somehow.
% So this function: loads control events for given subject, calculates
% amplitude envelope and its first derivative, groups peaks of these 2
% features based on correct vs wrong vs no responses, then plots and
% compares these.

%% Select blocks to import
vars=who;
if ~ismember(vars,'control_blocks')
    control_blocks = select_cont_blocks(Sbj_Metadata);
elseif isempty(control_blocks)
    control_blocks = select_cont_blocks(Sbj_Metadata);
end

save_folder = fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),'SP_comp_resp_acoustics');
if ~exist(save_folder,'dir'),mkdir(save_folder),end
clear vars

%% bring in these control events
% load control events
fprintf('Loading control events from:\n\t->%s\n',fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),[strjoin(control_blocks,'_') '_control_wltERP.mat']))
load(fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),[strjoin(control_blocks,'_') '_control_wltERP.mat']),'control_events')

%% Get the speech
Fs = 24000;
main_stim_loc = '/home/sakkol/Documents/BACKUPS/Speech_Perception_BACKUP(2019.09.11)/SpeechPerception_Task(backup-11.09.2019)/SpeechPerception_Task/main_stim_loc';

% speech rate
[filepath,~,~] = fileparts(control_events.trial_details{1}.speech.file);
name = strsplit(filepath,filesep);
speech_rate = name{end}(end-2:end);

% configuration
cfg = [];
cfg.SNR = -4;
cfg.speech.noise = 'silence';
cfg.LvsR = 'L';
cfg.delay = 0;

corrRate=[];corrEnv=[];
noRate=[];noEnv=[];
wrngRate=[];wrngEnv=[];

for c = 1:size(control_events,1)
    % get speech
    cfg.speech.file = find_sentence(control_events.sentence{c},main_stim_loc,speech_rate);
    [speech,~]=stim_creatorv2(cfg);
    
    % get envelope and peaks events
    [peakRate, peakEnv, amp_env, deriv_amp_env] = get_speech_peaks(speech(:,1),Fs,0);
    
    % loop words
    for w=1:5
        
        begsec = control_events.word_info{c}.onset(w);
        endsec = control_events.word_info{c}.offset(w);
        
        maxRate = max(deriv_amp_env(peakRate(peakRate>begsec*Fs&peakRate<endsec*Fs)));
        maxEnv = max(amp_env(peakEnv(peakEnv>begsec*Fs&peakEnv<endsec*Fs)));
        
        if isempty(maxRate) || isempty(maxEnv)
            maxRate = max(deriv_amp_env(peakRate(peakRate>(begsec-0.1)*Fs&peakRate<(endsec+0.1)*Fs)));
            maxEnv = max(amp_env(peakEnv(peakEnv>(begsec-0.1)*Fs&peakEnv<(endsec+0.1)*Fs)));
        end
        
        if maxRate==0 || maxEnv == 0
            fprintf('Trial:%d;Word:%d\n',c,w)
        end
        
        % check what is the response when that word was heard
        wword_resp = control_events.word_info{c}.response{w};

        if strcmp(wword_resp,'1')
            corrRate(end+1) = maxRate;
            corrEnv(end+1) = maxEnv;
        elseif strcmp(wword_resp,'0')
            noRate(end+1) = maxRate;
            noEnv(end+1) = maxEnv;
        else % wrong responses
            wrngRate(end+1) = maxRate;
            wrngEnv(end+1) = maxEnv;
        end        
    end
    
end


%% Plot
figure('Units','normalized','Position', [0 0  .3 .5]);

% first for Envelope
subplot(1,2,1)
arr_data = [corrEnv,noEnv,wrngEnv];
arr_loc = [ones(1,length(corrEnv)),ones(1,length(noEnv))*2,ones(1,length(wrngEnv))*3];
boxplot(arr_data,arr_loc)
hold on
for i=1:3
    scatter(arr_loc(arr_loc==i)',arr_data(arr_loc==i)','filled','jitter','on','MarkerFaceColor','k')
end

[~,~,stats] = anova1(arr_data,arr_loc,'off');
[c,~,~,~] = multcompare(stats,'Display','off');

% plot significance
y=ylim;
ylim([y(1) y(2)+y(2)*0.15]);y=ylim;
norm_y=(y(2)-y(1))/20;
for i = 1:size(c,1)
    if c(i,end)<0.05
        ccol = 'r';
    else
        ccol = 'k';
    end
    plot([c(i,1) c(i,2)],[y(2) - i*norm_y, y(2) - i*norm_y],ccol)
    text((c(i,1)+c(i,2))/2,y(2) - i*norm_y + norm_y/2,num2str(c(i,end)),'Color',ccol,'HorizontalAlignment','center')
end
set(gca,'FontSize',10)

title('Max envelope events','FontSize',16,'FontWeight','bold')
ylabel('Norm amplitude','FontSize',16,'FontWeight','bold')
xticklabels({'Correct-resp','No-resp','Wrong-resp'})


% second for Rate
subplot(1,2,2)
arr_data = [corrRate,noRate,wrngRate];
arr_loc = [ones(1,length(corrRate)),ones(1,length(noRate))*2,ones(1,length(wrngRate))*3];
boxplot(arr_data,arr_loc)
hold on
for i=1:3
    scatter(arr_loc(arr_loc==i)',arr_data(arr_loc==i)','filled','jitter','on','MarkerFaceColor','k')
end

[~,~,stats] = anova1(arr_data,arr_loc,'off');
[c,~,~,~] = multcompare(stats,'Display','off');

% plot significance
y=ylim;
ylim([y(1) y(2)+y(2)*0.15]);y=ylim;
norm_y=(y(2)-y(1))/20;
for i = 1:size(c,1)
    if c(i,end)<0.05
        ccol = 'r';
    else
        ccol = 'k';
    end
    plot([c(i,1) c(i,2)],[y(2) - i*norm_y, y(2) - i*norm_y],ccol)
    text((c(i,1)+c(i,2))/2,y(2) - i*norm_y + norm_y/2,num2str(c(i,end)),'Color',ccol,'HorizontalAlignment','center')
end
set(gca,'FontSize',10)

title('Max rate events','FontSize',16,'FontWeight','bold')
ylabel('Norm amplitude','FontSize',16,'FontWeight','bold')
xticklabels({'Correct-resp','No-resp','Wrong-resp'})

% save the figure
fprintf('Finished acoustics comparison for %s\n',Sbj_Metadata.sbj_ID)
print(fullfile(save_folder,'Compare_acoustics.jpg'),'-djpeg','-r300')


end