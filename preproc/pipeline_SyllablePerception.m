%% Admin
% Prepare
data_root = '/media/sakkol/HDD1/HBML/';
project_name = 'SyllablePerception';
sbj_ID = 'NS163';
Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'

% Get params directly from BlockList excel sheet
curr_block = Sbj_Metadata.BlockLists{1}
params = AllBlockInfo2params(Sbj_Metadata,curr_block)

%% Prepare ecog
if ~strcmpi(params.CurrBlockInfo.EEGDAT,'edf')
    [ecog] = TDT2ecog(params);
else
    % for edfs
    params.edf = fullfile(params.directory,[params.filename, '.edf']);
    params.analog.ttl            = 'DC1';
    % params.analog.audio          = 'DC7';
    % params.analog.micro          = 'C210';
    ecog=edf2ecog(params);
    tmp = ecog.ftrip.time;ecog.ftrip.time={[]};ecog.ftrip.time{1}=tmp;clear tmp
    
    % not needed, just can stay here as a reminder: data = edf2fieldtrip(params.edf)
end

%% Find bad channels using PSD
% find_bad_chans(ecog);
ecog = bad_chan_GUI([],ecog);
if size(ecog.bad_chans,2)>1,ecog.bad_chans = ecog.bad_chans';end % if bad_chans are in row order.

% Look at EEG, write bad channels to Xls-sheet
% cfg=[];
% cfg.channel = ecog.bad_chans;
% ecog_databrowser(ecog)

%% Check the good channels (write out bad channels, manually add them to ecog.bad_chans)
cfg = [];
cfg.viewmode = 'vertical';
cfg.blocksize = 20;
cfg.preproc.bsfilter       = 'yes';
cfg.preproc.bsfiltord      = 3;
cfg.preproc.bsfreq         = [59 61; 119 121; 179 181];
% cfg.preproc.bsfreq         = [59 61; 119 121; 179 181; 200 1000]; %sabina (visualisation purposes only, for noisy data)
cfg.preproc.demean         = 'yes';
good_chns = get_good_chans(ecog,2);
temp = ecog.ftrip;
temp.label = ecog.ftrip.label(good_chns);
temp.trial{1} = ecog.ftrip.trial{1}(good_chns,:);
temp.nChans = length(good_chns);
cfg = ft_databrowser(cfg, temp);
clear temp good_chns

%% If you added bad (or SOZ, spikey, out) chans to xls. Read xls in again
% this will overwrite several fields in the ecog structure
% [labelfile,ecog] = read_labelfile(params.labelfile,ecog);
% if you don't want to be asked if you'd want to overwrite.
[labelfile,ecog] = read_labelfile(params.labelfile,ecog,1); 

%% Write the bad channels seen in this block
% ecog.bad_chans = [ecog.bad_chans; {'RPs15';'RFp2';'RFp3';'RFp4';'RHs11';'RHs13';'RPc13'}];%;'Second';'In Column Order'}];
select_bad_chans

save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog.mat']),'ecog','-v7.3');

%% Get events
% first select edf or tdt analog channel
if isfield(ecog.analog,'trial')
    analog_struct = ecog.analog;
else % for edfs, where 
    analog_struct.trial{1} = ecog.analog.ttl;
    analog_struct.fs = ecog.analog.fs;
end

% If needed: check for threshold
figure; plot(analog_struct.trial{1}(1,:));
title('Check for threshold');

% Analog2digital of the noise channel
if isfield(ecog.analog,'trial') % amplitude threshold
    thr_ampl = 0.01; % amplitude threshold
    noise_ch = analog_struct.trial{1}(1,:);
else % for edfs 
    thr_ampl = 500000;
    noise_ch = demean(analog_struct.trial{1}(1,:));
end
refract_tpts = floor(.5*analog_struct.fs); % .5 second only
digital_trig_chan=analog2digital_trig(noise_ch,thr_ampl,refract_tpts,0);

trial_onsets_tpts=find(digital_trig_chan==1);

% plot events
figure('Units','normalized','Position', [0 0  1 .5]);
plot(noise_ch); hold on
for i=1:length(trial_onsets_tpts)
    plot([trial_onsets_tpts(i) trial_onsets_tpts(i)],ylim)
end
title([num2str(length(trial_onsets_tpts)) ' onsets found']);
xlim([1 length(noise_ch)])
print(fullfile(Sbj_Metadata.iEEG_data,curr_block,'PICS','events.jpg'),'-djpeg','-r300')

response_table = readtable(fullfile(Sbj_Metadata.behavioral_root, curr_block, [curr_block '_response_sheet.xlsx']));
tmp = load(fullfile(Sbj_Metadata.behavioral_root, curr_block, [curr_block '.mat']));

onset = zeros(length(trial_onsets_tpts),1);
duration = zeros(length(trial_onsets_tpts),1);
SyllablePresented = cell(length(trial_onsets_tpts),1);
Response = cell(length(trial_onsets_tpts),1);

for i=1:length(trial_onsets_tpts)
    onset(i) = trial_onsets_tpts(i)/ecog.analog.fs;
    SyllablePresented{i,1} = response_table.SyllablePresented{i};
    Response{i,1} = response_table.Response{i};
end

events = [array2table(onset) array2table(duration) cell2table(SyllablePresented) cell2table(Response)];
events = [events,tmp.events_table(1:length(trial_onsets_tpts),:)];
for i=1:size(events.duration)
    events.duration(i) = length(events.trials{i})/44100;
end

%% trial based rejection
fs = ecog.ftrip.fsample;
% Notch filter, demean
cfg=[];
cfg.demean         = 'yes';
cfg.bsfilter       = 'yes';
cfg.bsfiltord      = 3;
cfg.bsfreq         = [59 61; 119 121; 179 181];
cont_notched = ft_preprocessing(cfg,ecog.ftrip);

% speech onset locked
pre = 0; % seconds 
post = 5; % seconds 
trl_trg           = [];
trl_trg(:,1)      = floor( events.onset*fs - fs*pre );
trl_trg(:,2)      = floor( events.onset*fs + fs*post );
% trl_trg(:,2)      = floor( events.response_onset*fs );
trl_trg(:,3)      = floor( -pre*fs );

% Epoch
cfg      = [];
cfg.trl  = trl_trg;
trials_all  = ft_redefinetrial(cfg,cont_notched);

cfg=[];
cfg.method = 'summary';
cfg.keepchannel = 'yes';
good_chns = get_good_chans(ecog,2);
cfg.channel = good_chns; % assuming bad channels are the same across blocks!!
cfg.metric = 'zvalue';
cfg.keeptrial = 'no';
trials_clean = ft_rejectvisual(cfg,trials_all);

%% if there were bad trials, take note here and remove them from further saving
% bad_trials_idx = [3, 5];
% all_idx = (1:length(trl_trg))';
% good_trials_idx = all_idx(~ismember(all_idx,bad_trials_idx));
% events = events(good_trials_idx,:);

%% Create info.mat file
cfg=[];
cfg.events = events;
cfg.ecog = ecog;
cfg.elecInfo = ecog.params.labelfile;
info = create_info(cfg);
save(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']),'info');

%% 1.Save raw (unaveraged) ecog 2.average and bipolarize 3.epoch into trials
save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog.mat']),'ecog','-v7.3');

% Re-reference - epoch - save
% Average ref
ecog.ftrip = cont_notched; % nothched or not-noteched
plot_stuff=0;
ignore_szr_chans=1;
ecog_avg=ecog_avg_ref(ecog,plot_stuff,ignore_szr_chans);
save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_avg.mat']),'ecog_avg','-v7.3');

% Also BIPOLAR reference
ecog_bp=ecog_bipolarize(ecog);
save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_bp.mat']),'ecog_bp','-v7.3');

%% Wavelet analysis
clearvars -except Sbj_Metadata curr_block
load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_avg.mat']))
load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']))
events = info.events; clear info

% speech onset locked
pre  = 1; % seconds 
post = 8; % seconds 
[epoched_data, epoched_wlt] = master_waveletTF(ecog_avg.ftrip,events.onset,pre,post);

% save wlt and data
save(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_wlt.mat']),'epoched_wlt','epoched_data','-v7.3');

%% Behavioral analysis: accuracy of each block (raw trial numbers)
load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']))
syllable_couples = {'pa','ba';'ta','da';'fa','va';'sa','za'};
unvoiced_syll = {'pa';'ta';'fa';'sa'};
voiced_syll = {'ba';'da';'va';'za'};

info.events.Accr = num2cell(strcmp(info.events.SyllablePresented,info.events.Response));
info.events.Accr(strcmp(info.events.Response,'')) = {'noresp'};

consonant = cell(size(info.events,1),1);
for t=1:size(info.events,1)
    if ismember(info.events.SyllablePresented{t},voiced_syll)
        consonant{t} = 'voiced';
    else
        consonant{t} = 'unvoiced';
    end
    if strcmp(info.events.SyllablePresented{t},info.events.Response{t})
        info.events.Accr{t} = num2str(strcmp(info.events.SyllablePresented{t},info.events.Response{t}));
    elseif strcmp(info.events.Response{t},'')
        info.events.Accr{t} = 'noresp';
    else
        info.events.Accr{t} = '0';
    end
end
info.events = [info.events,cell2table(consonant)];

if strcmp(Sbj_Metadata.sbj_ID,'NS163') && strcmp(curr_block,'SyllablePerception_EStim')
    info.events = info.events(~ismember(1:64,[1,11,14,16,21,22,25,31,32,33,37,41,43,51,54,58,61,62,63,64]),:);
end

voiced_corr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'1'));
voiced_incorr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'0'));
voiced_noresp = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'noresp'));

unvoiced_corr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'1'));
unvoiced_incorr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'0'));
unvoiced_noresp = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'noresp'));

figure('Position',[100 100 1000 1000])
to_bar = [voiced_corr voiced_incorr voiced_noresp;...
          unvoiced_corr unvoiced_incorr unvoiced_noresp];
hBar = bar(to_bar,'stacked');
set(gca,'xticklabels',{'Voiced','Unvoiced'},'FontSize',15)
grid on
legend('Correct response','Incorrect Response','No-response')                                                  % Return ‘bar’ Handle
for k1 = 1:size(to_bar,2)
    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, hBar(k1).XOffset');    % Note: ‘XOffset’ Is An Undocumented Feature, This Selects The ‘bar’ Centres
    ydt(k1,:) = hBar(k1).YData;                                     % Individual Bar Heights
    if k1==1
        text(ctr(k1,:), ydt(k1,:), sprintfc('%.0f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',18, 'Color','k')
    else
        text(ctr(k1,:), ydt(k1,:)+sum(ydt(1:k1-1,:),1), sprintfc('%.0f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',18, 'Color','k')
    end
end

% Observed data
n1 = voiced_corr; N1 = voiced_incorr+voiced_noresp+voiced_corr;
n2 = unvoiced_corr; N2 = unvoiced_incorr+unvoiced_noresp+unvoiced_corr;
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval] = crosstab(x1,x2)
title(['Comparing accuracy between voiced vs unvoiced, p-val: ' num2str(pval)])

print(fullfile(Sbj_Metadata.results,[curr_block '_accuracy_onlyclean.jpg']),'-djpeg','-r300')

%% Behavioral analysis: accuracy of each block (percentage)
load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']))
syllable_couples = {'pa','ba';'ta','da';'fa','va';'sa','za'};
unvoiced_syll = {'pa';'ta';'fa';'sa'};
voiced_syll = {'ba';'da';'va';'za'};

info.events.Accr = num2cell(strcmp(info.events.SyllablePresented,info.events.Response));
info.events.Accr(strcmp(info.events.Response,'')) = {'noresp'};

consonant = cell(size(info.events,1),1);
for t=1:size(info.events,1)
    if ismember(info.events.SyllablePresented{t},voiced_syll)
        consonant{t} = 'voiced';
    else
        consonant{t} = 'unvoiced';
    end
    if strcmp(info.events.SyllablePresented{t},info.events.Response{t})
        info.events.Accr{t} = num2str(strcmp(info.events.SyllablePresented{t},info.events.Response{t}));
    elseif strcmp(info.events.Response{t},'')
        info.events.Accr{t} = 'noresp';
    else
        info.events.Accr{t} = '0';
    end
end
info.events = [info.events,cell2table(consonant)];

if strcmp(Sbj_Metadata.sbj_ID,'NS163') && strcmp(curr_block,'SyllablePerception_EStim')
    info.events = info.events(~ismember(1:64,[1,11,14,16,21,22,25,31,32,33,37,41,43,51,54,58,61,62,63,64]),:);
end

voiced_corr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'1'));
voiced_incorr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'0'));
voiced_noresp = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'noresp'));
voiced_total = voiced_corr +voiced_incorr+voiced_noresp;
unvoiced_corr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'1'));
unvoiced_incorr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'0'));
unvoiced_noresp = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'noresp'));
unvoiced_total = unvoiced_corr + unvoiced_incorr+unvoiced_noresp;

figure('Position',[100 100 1000 1000])
to_bar = 100*[voiced_corr/voiced_total voiced_incorr/voiced_total voiced_noresp/voiced_total;...
          unvoiced_corr/unvoiced_total unvoiced_incorr/unvoiced_total unvoiced_noresp/unvoiced_total];
hBar = bar(to_bar,'stacked');
set(gca,'xticklabels',{'Voiced','Unvoiced'},'FontSize',15)
grid on
legend('Correct response','Incorrect Response','No-response','Location','southoutside','NumColumns',3)
ctr=[];ydt=[];
for k1 = 1:size(to_bar,2)
    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, hBar(k1).XOffset');    % Note: ‘XOffset’ Is An Undocumented Feature, This Selects The ‘bar’ Centres
    ydt(k1,:) = hBar(k1).YData;                                     % Individual Bar Heights
    if k1==1
        text(ctr(k1,:), ydt(k1,:), sprintfc('%.0f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',18, 'Color','k')
    else
        text(ctr(k1,:), ydt(k1,:)+sum(ydt(1:k1-1,:),1), sprintfc('%.0f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',18, 'Color','k')
    end
end

% Observed data
n1 = voiced_corr; N1 = voiced_incorr+voiced_noresp+voiced_corr;
n2 = unvoiced_corr; N2 = unvoiced_incorr+unvoiced_noresp+unvoiced_corr;
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval] = crosstab(x1,x2)
title(['Comparing accuracy between voiced vs unvoiced, p-val: ' num2str(pval)])
ylim([0 100])

print(fullfile(Sbj_Metadata.results,[curr_block '_accuracy_onlyclean_percentage.jpg']),'-djpeg','-r300')

%% Comparing between blocks raw trial numbers
syllable_couples = {'pa','ba';'ta','da';'fa','va';'sa','za'};
unvoiced_syll = {'pa';'ta';'fa';'sa'};
voiced_syll = {'ba';'da';'va';'za'};
block=cell(1,2);
for i=1:2
    load(fullfile(Sbj_Metadata.iEEG_data,Sbj_Metadata.BlockLists{i},[Sbj_Metadata.BlockLists{i} '_info.mat']));
    info.events.Accr = num2cell(strcmp(info.events.SyllablePresented,info.events.Response));
    info.events.Accr(strcmp(info.events.Response,'')) = {'noresp'};
    
    consonant = cell(size(info.events,1),1);
    for t=1:size(info.events,1)
        if ismember(info.events.SyllablePresented{t},voiced_syll)
            consonant{t} = 'voiced';
        else
            consonant{t} = 'unvoiced';
        end
        if strcmp(info.events.SyllablePresented{t},info.events.Response{t})
            info.events.Accr{t} = num2str(strcmp(info.events.SyllablePresented{t},info.events.Response{t}));
        elseif strcmp(info.events.Response{t},'')
            info.events.Accr{t} = 'noresp';
        else
            info.events.Accr{t} = '0';
        end
    end
    info.events = [info.events,cell2table(consonant)];
    if strcmp(Sbj_Metadata.sbj_ID,'NS163') && strcmp(Sbj_Metadata.BlockLists{i},'SyllablePerception_EStim')
        info.events = info.events(~ismember(1:64,[1,11,14,16,21,22,25,31,32,33,37,41,43,51,54,58,61,62,63,64]),:);
    end
    
    voiced_corr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'1'));
    voiced_incorr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'0'));
    voiced_noresp = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'noresp'));
    
    unvoiced_corr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'1'));
    unvoiced_incorr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'0'));
    unvoiced_noresp = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'noresp'));
    
    block{i}.voiced_corr = voiced_corr;
    block{i}.voiced_incorr = voiced_incorr;
    block{i}.voiced_noresp = voiced_noresp;
    block{i}.unvoiced_corr = unvoiced_corr;
    block{i}.unvoiced_incorr = unvoiced_incorr;
    block{i}.unvoiced_noresp = unvoiced_noresp;
end


figure('Position',[100 100 1200 1000])
to_bar = [block{1}.voiced_corr block{1}.voiced_incorr block{1}.voiced_noresp;...
          block{2}.voiced_corr block{2}.voiced_incorr block{2}.voiced_noresp;...
          block{1}.unvoiced_corr block{1}.unvoiced_incorr block{1}.unvoiced_noresp;...
          block{2}.unvoiced_corr block{2}.unvoiced_incorr block{2}.unvoiced_noresp];
hBar = bar(to_bar,'stacked');
set(gca,'xticklabels',{'Voiced-NoEStim','Voiced-EStim','Unvoiced-NoEStim','Unvoiced-EStim'},'FontSize',15)
grid on
legend('Correct response','Incorrect Response','No-response')                                                  % Return ‘bar’ Handle
ctr=[];ydt=[];
for k1 = 1:size(to_bar,2)
    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, hBar(k1).XOffset');    % Note: ‘XOffset’ Is An Undocumented Feature, This Selects The ‘bar’ Centres
    ydt(k1,:) = hBar(k1).YData;                                     % Individual Bar Heights
    if k1==1
        text(ctr(k1,:), ydt(k1,:), sprintfc('%.0f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',18, 'Color','k')
    else
        text(ctr(k1,:), ydt(k1,:)+sum(ydt(1:k1-1,:),1), sprintfc('%.0f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',18, 'Color','k')
    end
end

% Observed data
n1 = block{1}.voiced_corr; N1 = block{1}.voiced_corr+block{1}.voiced_incorr+block{1}.voiced_noresp;
n2 = block{2}.voiced_corr; N2 = block{2}.voiced_corr+block{1}.voiced_incorr+block{1}.voiced_noresp;
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl1,chi2stat1,pval1] = crosstab(x1,x2)

n1 = block{1}.unvoiced_corr; N1 = block{1}.unvoiced_corr+block{1}.unvoiced_incorr+block{1}.unvoiced_noresp;
n2 = block{2}.unvoiced_corr; N2 = block{2}.unvoiced_corr+block{1}.unvoiced_incorr+block{1}.unvoiced_noresp;
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl2,chi2stat2,pval2] = crosstab(x1,x2)
title({['Comparing accuracy of voiced between no-EStim vs EStim, p-val: ' num2str(pval1)];...
       ['Comparing accuracy of unvoiced between no-EStim vs EStim, p-val: ' num2str(pval2)]})

print(fullfile(Sbj_Metadata.results,[curr_block '_block_comp_onlyclean.jpg']),'-djpeg','-r300')
% print(fullfile(Sbj_Metadata.results,[curr_block '_block_comp.jpg']),'-djpeg','-r300')

%% Comparing between blocks with percentages, NOT raw trial numbers
syllable_couples = {'pa','ba';'ta','da';'fa','va';'sa','za'};
unvoiced_syll = {'pa';'ta';'fa';'sa'};
voiced_syll = {'ba';'da';'va';'za'};
block=cell(1,2);
for i=1:2
    load(fullfile(Sbj_Metadata.iEEG_data,Sbj_Metadata.BlockLists{i},[Sbj_Metadata.BlockLists{i} '_info.mat']));
    info.events.Accr = num2cell(strcmp(info.events.SyllablePresented,info.events.Response));
    info.events.Accr(strcmp(info.events.Response,'')) = {'noresp'};
    
    consonant = cell(size(info.events,1),1);
    for t=1:size(info.events,1)
        if ismember(info.events.SyllablePresented{t},voiced_syll)
            consonant{t} = 'voiced';
        else
            consonant{t} = 'unvoiced';
        end
        if strcmp(info.events.SyllablePresented{t},info.events.Response{t})
            info.events.Accr{t} = num2str(strcmp(info.events.SyllablePresented{t},info.events.Response{t}));
        elseif strcmp(info.events.Response{t},'')
            info.events.Accr{t} = 'noresp';
        else
            info.events.Accr{t} = '0';
        end
    end
    info.events = [info.events,cell2table(consonant)];
    % remove bad trials in e-stim block (the ones that extend into speaking or the ones that are not stimulating during the stim presentation)
    if strcmp(Sbj_Metadata.sbj_ID,'NS163') && strcmp(Sbj_Metadata.BlockLists{i},'SyllablePerception_EStim')
        info.events = info.events(~ismember(1:64,[1,11,14,16,21,22,25,31,32,33,37,41,43,51,54,58,61,62,63,64]),:);
    end
    
    voiced_corr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'1'));
    voiced_incorr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'0'));
    voiced_noresp = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'noresp'));
    
    unvoiced_corr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'1'));
    unvoiced_incorr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'0'));
    unvoiced_noresp = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'noresp'));
    
    voiced_total = voiced_corr+voiced_incorr+voiced_noresp;
    unvoiced_total = unvoiced_corr+unvoiced_incorr+unvoiced_noresp;
    block{i}.voiced_corr = voiced_corr;
    block{i}.voiced_incorr = voiced_incorr;
    block{i}.voiced_noresp = voiced_noresp;
    block{i}.voiced_total = voiced_total;
    block{i}.unvoiced_corr = unvoiced_corr;
    block{i}.unvoiced_incorr = unvoiced_incorr;
    block{i}.unvoiced_noresp = unvoiced_noresp;
    block{i}.unvoiced_total = unvoiced_total;
end


figure('Position',[100 100 1200 1000])
to_bar = [block{1}.voiced_corr/block{1}.voiced_total block{1}.voiced_incorr/block{1}.voiced_total block{1}.voiced_noresp/block{1}.voiced_total;...
          block{2}.voiced_corr/block{2}.voiced_total block{2}.voiced_incorr/block{2}.voiced_total block{2}.voiced_noresp/block{2}.voiced_total;...
          block{1}.unvoiced_corr/block{1}.unvoiced_total block{1}.unvoiced_incorr/block{1}.unvoiced_total block{1}.unvoiced_noresp/block{1}.unvoiced_total;...
          block{2}.unvoiced_corr/block{2}.unvoiced_total block{2}.unvoiced_incorr/block{2}.unvoiced_total block{2}.unvoiced_noresp/block{2}.unvoiced_total]*100;
hBar = bar(to_bar,'stacked');
set(gca,'xticklabels',{'Voiced-NoEStim','Voiced-EStim','Unvoiced-NoEStim','Unvoiced-EStim'},'FontSize',15)
grid on
legend('Correct response','Incorrect Response','No-response','Location','southoutside','NumColumns',3)
ctr=[];ydt=[];
for k1 = 1:size(to_bar,2)
    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, hBar(k1).XOffset');    % Note: ‘XOffset’ Is An Undocumented Feature, This Selects The ‘bar’ Centres
    ydt(k1,:) = hBar(k1).YData;                                     % Individual Bar Heights
    if k1==1
        text(ctr(k1,:), ydt(k1,:), sprintfc('%.0f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',18, 'Color','k')
    else
        text(ctr(k1,:), ydt(k1,:)+sum(ydt(1:k1-1,:),1), sprintfc('%.0f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',18, 'Color','k')
    end
end

% Observed data
n1 = block{1}.voiced_corr; N1 = block{1}.voiced_corr+block{1}.voiced_incorr+block{1}.voiced_noresp;
n2 = block{2}.voiced_corr; N2 = block{2}.voiced_corr+block{1}.voiced_incorr+block{1}.voiced_noresp;
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl1,chi2stat1,pval1] = crosstab(x1,x2)

n1 = block{1}.unvoiced_corr; N1 = block{1}.unvoiced_corr+block{1}.unvoiced_incorr+block{1}.unvoiced_noresp;
n2 = block{2}.unvoiced_corr; N2 = block{2}.unvoiced_corr+block{1}.unvoiced_incorr+block{1}.unvoiced_noresp;
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl2,chi2stat2,pval2] = crosstab(x1,x2)
title({['Comparing accuracy of voiced between no-EStim vs EStim, p-val: ' num2str(pval1)];...
       ['Comparing accuracy of unvoiced between no-EStim vs EStim, p-val: ' num2str(pval2)]})

print(fullfile(Sbj_Metadata.results,[curr_block '_block_comp_onlyclean_percent.jpg']),'-djpeg','-r300')
% print(fullfile(Sbj_Metadata.results,[curr_block '_block_comp.jpg']),'-djpeg','-r300')


%% Plotting wavelet plot (with bad trials)
load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_wlt.mat']),'epoched_wlt','epoched_data');
%  Baseline correct time-freq data
load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']))
syllable_couples = {'pa','ba';'ta','da';'fa','va';'sa','za'};
unvoiced_syll = {'pa';'ta';'fa';'sa'};
voiced_syll = {'ba';'da';'va';'za'};

info.events.Accr = num2cell(strcmp(info.events.SyllablePresented,info.events.Response));
info.events.Accr(strcmp(info.events.Response,'')) = {'noresp'};
consonant = cell(size(info.events,1),1);
for t=1:size(info.events,1)
    if ismember(info.events.SyllablePresented{t},voiced_syll)
        consonant{t} = 'voiced';
    else
        consonant{t} = 'unvoiced';
    end
    if strcmp(info.events.SyllablePresented{t},info.events.Response{t})
        info.events.Accr{t} = num2str(strcmp(info.events.SyllablePresented{t},info.events.Response{t}));
    elseif strcmp(info.events.Response{t},'')
        info.events.Accr{t} = 'noresp';
    else
        info.events.Accr{t} = '0';
    end
end
info.events = [info.events,cell2table(consonant)];

% shorten a bit
cfg=[];
cfg.latency = [-1 4];
epoched_wlt = ft_selectdata(cfg,epoched_wlt);
epoched_data = ft_selectdata(cfg,epoched_data);
% from fourierspectrum to powerspectrum
cfg = [];
cfg.output='abs';
cfg.keeptrials = 'yes';
epoched_wlt=ft_freqdescriptives(cfg,epoched_wlt);
epoched_wlt.dimord = 'rpt_chan_freq_time';
cfg              = [];
cfg.baseline     = [-.5 -0.1];
cfg.baselinetype = 'relative';
cfg.parameter    = 'powspctrm';
[epoched_wlt]         = ft_freqbaseline(cfg, epoched_wlt);

% HFA calc
cfg=[];
cfg.frequency = [70 150];
cfg.avgoverfreq = 'yes';
cfg.nanmean = 'yes';
epoched_HFA_bc = ft_selectdata(cfg, epoched_wlt);
[epoched_HFA_bc]   = master_smoothData(epoched_HFA_bc,[]);

voiced_events = find(ismember(info.events.SyllablePresented,voiced_syll) & ismember(info.events.Accr,'1'))';
unvoiced_events = find(ismember(info.events.SyllablePresented,unvoiced_syll) & ismember(info.events.Accr,'1'))';
% remove the first 38 trials
voiced_events = voiced_events(~ismember(voiced_events,[1:43,69,73]));
unvoiced_events = unvoiced_events(~ismember(unvoiced_events,[1:43,69,73]));

totitle1 = 'Voiced syllables';
totitle2 = 'Unvoiced syllables';
syll_time = 0.5;
response_time = mean(info.events.duration);

cfg=[];
cfg.trials = voiced_events;
voiced_HFA_bc = ft_selectdata(cfg, epoched_HFA_bc);
voiced_ERP = ft_selectdata(cfg, epoched_data);

cfg=[];
cfg.trials = unvoiced_events;
unvoiced_HFA_bc = ft_selectdata(cfg, epoched_HFA_bc);
unvoiced_ERP = ft_selectdata(cfg, epoched_data);


%% Plotting wavelet plot (withOUT bad trials)
bwr = load('bwr_cmap.mat');
save_folder = fullfile(Sbj_Metadata.results,'Quick_results_38removed',curr_block);
if ~exist(save_folder,'dir'),mkdir(save_folder),end


for el = 1:length(epoched_data.label)
    
    % 2x2 plot with spectrogram on bottom and HFA in the R-top, ERP in
    % L-top
    figure('Units','normalized','Position', [0 0  .7 1]);
    
    % plot ERP
    subplot(221)
    for_avg=smoothdata(squeeze(voiced_ERP.trial(:,el,:)),2,'gaussian',20);
    hold on
    shadedErrorBar(voiced_ERP.time,mean(for_avg,1),stderr(for_avg),'lineprops','b');
    for_avg=smoothdata(squeeze(unvoiced_ERP.trial(:,el,:)),2,'gaussian',20);
    shadedErrorBar(voiced_ERP.time,mean(for_avg,1),stderr(for_avg),'lineprops','r');
    title({'ERP (mean+SEM)';[num2str(size(unvoiced_ERP.trial,1)) ' ' totitle1 ' - ' num2str(size(voiced_ERP.trial,1)) ' ' totitle2]})
    legend({totitle1,totitle2},'Location','southeast')
    xlim([epoched_data.time(1) epoched_data.time(end)])
    p=plot(xlim,[0 0],'k');set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(gca, 'FontSize',12);
    ylims = ylim;
    
    p=plot([0 0], ylim,'k');set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    text(0,ylims(2)-(ylims(2)-ylims(1))/13,'Noise starts')
    p=plot([syll_time syll_time], ylim,'k');set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    text(syll_time,ylims(2)-(ylims(2)-ylims(1))/10,'Syllable onset')
    p=plot([response_time response_time], ylim,'k');set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    text(response_time,ylims(2)-(ylims(2)-ylims(1))/10,'Avg.Resp.Onset')
    set(gca, 'FontSize',13,'FontWeight','bold');
    xlabel('Time (s)');
    ylim(ylims)
    ylabel('Voltage (uV)')
    
    % plot HFA
    subplot(222)
    for_avg=smoothdata(squeeze(voiced_HFA_bc.powspctrm(:,el,:,:)),2,'gaussian',15);
    hold on
    shadedErrorBar(voiced_HFA_bc.time,mean(for_avg,1),stderr(for_avg),'lineprops','b');
    for_avg=smoothdata(squeeze(unvoiced_HFA_bc.powspctrm(:,el,:,:)),2,'gaussian',15);
    shadedErrorBar(voiced_HFA_bc.time,mean(for_avg,1),stderr(for_avg),'lineprops','r');
    title({'HFA (mean+SEM)';[num2str(size(unvoiced_HFA_bc.powspctrm,1)) ' ' totitle1 ' - ' num2str(size(voiced_HFA_bc.powspctrm,1)) ' ' totitle2]})
    plot(xlim,[1 1],'k')
    legend({totitle1,totitle2},'Location','southeast')
    xlim([epoched_data.time(1) epoched_data.time(end)])
    set(gca, 'FontSize',12);
    ylim([0.5 3])
    ylims = ylim;
    
    p=plot([0 0], ylim,'k');set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    text(0,ylims(2)-(ylims(2)-ylims(1))/13,'Noise starts')
    p=plot([syll_time syll_time], ylim,'k');set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    text(syll_time,ylims(2)-(ylims(2)-ylims(1))/10,'Syllable onset')
    p=plot([response_time response_time], ylim,'k');set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    text(response_time,ylims(2)-(ylims(2)-ylims(1))/10,'Avg.Resp.Onset')
    set(gca, 'FontSize',13,'FontWeight','bold');
    xlabel('Time (s)');
    ylim(ylims)
    ylabel('HFA power (relative)')
    
    % plot spectrograms
    for s=1:2
        subplot(2,2,s+2)
        if s==1
            avg_spec = squeeze(nanmean(epoched_wlt.powspctrm(voiced_events,el,:,:),1));
        else
            avg_spec = squeeze(nanmean(epoched_wlt.powspctrm(unvoiced_events,el,:,:),1));
        end
        h = pcolor(epoched_wlt.time,epoched_wlt.freq,avg_spec);
        h.EdgeColor = 'none';
        axis xy
        set(gca,'YScale','log')
        yticks([4,8,12,20,50,70,100,150,200])
        set(gca, 'FontSize',13,'FontWeight','bold');
        caxis([0 3]);
        hold on
        xlim([epoched_wlt.time(1) epoched_wlt.time(end)])
        hold on
        ylim([3 200])
        ylims = ylim;
        plot3([0, 0],ylims, [150 150],'k')
        text(0,ylims(2)-(ylims(2)-ylims(1))/13,'Noise starts')
        
        plot3([syll_time, syll_time],ylims, [150 150],'k')
        text(syll_time,ylims(2)-(ylims(2)-ylims(1))/10,'Syllable onset')
        
        plot3([response_time, response_time],ylims, [150 150],'k')
        text(response_time,ylims(2)-(ylims(2)-ylims(1))/10,'Avg.Resp.Onset')
        xlabel('Time (s)');
        ylabel('Frequency (Hz-log scale)');
        
        if s==1
            title(['Spectrogram for ' totitle1])
        else
            title(['Spectrogram for ' totitle2])
        end
        c=colorbar('southoutside');
        c.Label.String = 'Power (relative)';
    end
    
     colormap(bwr.rgb_vals);
    
     % main title
    sgtitle(['Elec: ' epoched_data.label{el} ' - Location: ' info.channelinfo.Desikan_Killiany{el}], 'FontSize',16,'FontWeight','bold')
    
    % Save the figure
    fprintf('\t-Saving electrode %s, out of %d\n',epoched_data.label{el},length(epoched_data.label))
    print(fullfile(save_folder,[epoched_data.label{el} , '_ERP-HFA-spectr.jpg']),'-djpeg')
    close all
end



%% SECOND HALF OF FIRST BLOCK
%% SECOND HALF OF FIRST BLOCK::

%% Behavioral analysis: accuracy of each block (raw trial numbers)
load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']))
syllable_couples = {'pa','ba';'ta','da';'fa','va';'sa','za'};
unvoiced_syll = {'pa';'ta';'fa';'sa'};
voiced_syll = {'ba';'da';'va';'za'};

info.events.Accr = num2cell(strcmp(info.events.SyllablePresented,info.events.Response));
info.events.Accr(strcmp(info.events.Response,'')) = {'noresp'};

consonant = cell(size(info.events,1),1);
for t=1:size(info.events,1)
    if ismember(info.events.SyllablePresented{t},voiced_syll)
        consonant{t} = 'voiced';
    else
        consonant{t} = 'unvoiced';
    end
    if strcmp(info.events.SyllablePresented{t},info.events.Response{t})
        info.events.Accr{t} = num2str(strcmp(info.events.SyllablePresented{t},info.events.Response{t}));
    elseif strcmp(info.events.Response{t},'')
        info.events.Accr{t} = 'noresp';
    else
        info.events.Accr{t} = '0';
    end
end
info.events = [info.events,cell2table(consonant)];

if strcmp(Sbj_Metadata.sbj_ID,'NS163') && strcmp(curr_block,'SyllablePerception_EStim')
    info.events = info.events(~ismember(1:64,[1,11,14,16,21,22,25,31,32,33,37,41,43,51,54,58,61,62,63,64]),:);
end
if strcmp(Sbj_Metadata.sbj_ID,'NS163') && strcmp(curr_block,'SyllablePerception')
    info.events = info.events(~ismember(1:92,1:46),:);
end

voiced_corr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'1'));
voiced_incorr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'0'));
voiced_noresp = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'noresp'));

unvoiced_corr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'1'));
unvoiced_incorr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'0'));
unvoiced_noresp = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'noresp'));

figure('Position',[100 100 1000 1000])
to_bar = [voiced_corr voiced_incorr voiced_noresp;...
          unvoiced_corr unvoiced_incorr unvoiced_noresp];
hBar = bar(to_bar,'stacked');
set(gca,'xticklabels',{'Voiced','Unvoiced'},'FontSize',15)
grid on
legend('Correct response','Incorrect Response','No-response')                                                  % Return ‘bar’ Handle
for k1 = 1:size(to_bar,2)
    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, hBar(k1).XOffset');    % Note: ‘XOffset’ Is An Undocumented Feature, This Selects The ‘bar’ Centres
    ydt(k1,:) = hBar(k1).YData;                                     % Individual Bar Heights
    if k1==1
        text(ctr(k1,:), ydt(k1,:), sprintfc('%.0f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',18, 'Color','k')
    else
        text(ctr(k1,:), ydt(k1,:)+sum(ydt(1:k1-1,:),1), sprintfc('%.0f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',18, 'Color','k')
    end
end

% Observed data
n1 = voiced_corr; N1 = voiced_incorr+voiced_noresp+voiced_corr;
n2 = unvoiced_corr; N2 = unvoiced_incorr+unvoiced_noresp+unvoiced_corr;
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval] = crosstab(x1,x2)
title(['Comparing accuracy between voiced vs unvoiced, p-val: ' num2str(pval)])

print(fullfile(Sbj_Metadata.results,['2ndhalf_' curr_block '_accuracy_onlyclean.jpg']),'-djpeg','-r300')

%% Behavioral analysis: accuracy of each block (percentage)
load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']))
syllable_couples = {'pa','ba';'ta','da';'fa','va';'sa','za'};
unvoiced_syll = {'pa';'ta';'fa';'sa'};
voiced_syll = {'ba';'da';'va';'za'};

info.events.Accr = num2cell(strcmp(info.events.SyllablePresented,info.events.Response));
info.events.Accr(strcmp(info.events.Response,'')) = {'noresp'};

consonant = cell(size(info.events,1),1);
for t=1:size(info.events,1)
    if ismember(info.events.SyllablePresented{t},voiced_syll)
        consonant{t} = 'voiced';
    else
        consonant{t} = 'unvoiced';
    end
    if strcmp(info.events.SyllablePresented{t},info.events.Response{t})
        info.events.Accr{t} = num2str(strcmp(info.events.SyllablePresented{t},info.events.Response{t}));
    elseif strcmp(info.events.Response{t},'')
        info.events.Accr{t} = 'noresp';
    else
        info.events.Accr{t} = '0';
    end
end
info.events = [info.events,cell2table(consonant)];

if strcmp(Sbj_Metadata.sbj_ID,'NS163') && strcmp(curr_block,'SyllablePerception_EStim')
    info.events = info.events(~ismember(1:64,[1,11,14,16,21,22,25,31,32,33,37,41,43,51,54,58,61,62,63,64]),:);
end
if strcmp(Sbj_Metadata.sbj_ID,'NS163') && strcmp(curr_block,'SyllablePerception')
    info.events = info.events(~ismember(1:92,1:46),:);
end

voiced_corr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'1'));
voiced_incorr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'0'));
voiced_noresp = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'noresp'));
voiced_total = voiced_corr +voiced_incorr+voiced_noresp;
unvoiced_corr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'1'));
unvoiced_incorr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'0'));
unvoiced_noresp = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'noresp'));
unvoiced_total = unvoiced_corr + unvoiced_incorr+unvoiced_noresp;

figure('Position',[100 100 1000 1000])
to_bar = 100*[voiced_corr/voiced_total voiced_incorr/voiced_total voiced_noresp/voiced_total;...
          unvoiced_corr/unvoiced_total unvoiced_incorr/unvoiced_total unvoiced_noresp/unvoiced_total];
hBar = bar(to_bar,'stacked');
set(gca,'xticklabels',{'Voiced','Unvoiced'},'FontSize',15)
grid on
legend('Correct response','Incorrect Response','No-response','Location','southoutside','NumColumns',3)
ctr=[];ydt=[];
for k1 = 1:size(to_bar,2)
    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, hBar(k1).XOffset');    % Note: ‘XOffset’ Is An Undocumented Feature, This Selects The ‘bar’ Centres
    ydt(k1,:) = hBar(k1).YData;                                     % Individual Bar Heights
    if k1==1
        text(ctr(k1,:), ydt(k1,:), sprintfc('%.0f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',18, 'Color','k')
    else
        text(ctr(k1,:), ydt(k1,:)+sum(ydt(1:k1-1,:),1), sprintfc('%.0f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',18, 'Color','k')
    end
end

% Observed data
n1 = voiced_corr; N1 = voiced_incorr+voiced_noresp+voiced_corr;
n2 = unvoiced_corr; N2 = unvoiced_incorr+unvoiced_noresp+unvoiced_corr;
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval] = crosstab(x1,x2)
title(['Comparing accuracy between voiced vs unvoiced, p-val: ' num2str(pval)])
ylim([0 100])

print(fullfile(Sbj_Metadata.results,['2ndhalf_' curr_block '_accuracy_onlyclean_percentage.jpg']),'-djpeg','-r300')

%% Comparing between blocks raw trial numbers
syllable_couples = {'pa','ba';'ta','da';'fa','va';'sa','za'};
unvoiced_syll = {'pa';'ta';'fa';'sa'};
voiced_syll = {'ba';'da';'va';'za'};
block=cell(1,2);
for i=1:2
    load(fullfile(Sbj_Metadata.iEEG_data,Sbj_Metadata.BlockLists{i},[Sbj_Metadata.BlockLists{i} '_info.mat']));
    info.events.Accr = num2cell(strcmp(info.events.SyllablePresented,info.events.Response));
    info.events.Accr(strcmp(info.events.Response,'')) = {'noresp'};
    
    consonant = cell(size(info.events,1),1);
    for t=1:size(info.events,1)
        if ismember(info.events.SyllablePresented{t},voiced_syll)
            consonant{t} = 'voiced';
        else
            consonant{t} = 'unvoiced';
        end
        if strcmp(info.events.SyllablePresented{t},info.events.Response{t})
            info.events.Accr{t} = num2str(strcmp(info.events.SyllablePresented{t},info.events.Response{t}));
        elseif strcmp(info.events.Response{t},'')
            info.events.Accr{t} = 'noresp';
        else
            info.events.Accr{t} = '0';
        end
    end
    info.events = [info.events,cell2table(consonant)];
    if strcmp(Sbj_Metadata.sbj_ID,'NS163') && strcmp(Sbj_Metadata.BlockLists{i},'SyllablePerception_EStim')
        info.events = info.events(~ismember(1:64,[1,11,14,16,21,22,25,31,32,33,37,41,43,51,54,58,61,62,63,64]),:);
    end
    if strcmp(Sbj_Metadata.sbj_ID,'NS163') && strcmp(Sbj_Metadata.BlockLists{i},'SyllablePerception')
        info.events = info.events(~ismember(1:92,1:46),:);
    end
    
    voiced_corr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'1'));
    voiced_incorr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'0'));
    voiced_noresp = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'noresp'));
    
    unvoiced_corr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'1'));
    unvoiced_incorr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'0'));
    unvoiced_noresp = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'noresp'));
    
    block{i}.voiced_corr = voiced_corr;
    block{i}.voiced_incorr = voiced_incorr;
    block{i}.voiced_noresp = voiced_noresp;
    block{i}.unvoiced_corr = unvoiced_corr;
    block{i}.unvoiced_incorr = unvoiced_incorr;
    block{i}.unvoiced_noresp = unvoiced_noresp;
end


figure('Position',[100 100 1200 1000])
to_bar = [block{1}.voiced_corr block{1}.voiced_incorr block{1}.voiced_noresp;...
          block{2}.voiced_corr block{2}.voiced_incorr block{2}.voiced_noresp;...
          block{1}.unvoiced_corr block{1}.unvoiced_incorr block{1}.unvoiced_noresp;...
          block{2}.unvoiced_corr block{2}.unvoiced_incorr block{2}.unvoiced_noresp];
hBar = bar(to_bar,'stacked');
set(gca,'xticklabels',{'Voiced-NoEStim','Voiced-EStim','Unvoiced-NoEStim','Unvoiced-EStim'},'FontSize',15)
grid on
legend('Correct response','Incorrect Response','No-response')                                                  % Return ‘bar’ Handle
ctr=[];ydt=[];
for k1 = 1:size(to_bar,2)
    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, hBar(k1).XOffset');    % Note: ‘XOffset’ Is An Undocumented Feature, This Selects The ‘bar’ Centres
    ydt(k1,:) = hBar(k1).YData;                                     % Individual Bar Heights
    if k1==1
        text(ctr(k1,:), ydt(k1,:), sprintfc('%.0f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',18, 'Color','k')
    else
        text(ctr(k1,:), ydt(k1,:)+sum(ydt(1:k1-1,:),1), sprintfc('%.0f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',18, 'Color','k')
    end
end

% Observed data
n1 = block{1}.voiced_corr; N1 = block{1}.voiced_corr+block{1}.voiced_incorr+block{1}.voiced_noresp;
n2 = block{2}.voiced_corr; N2 = block{2}.voiced_corr+block{1}.voiced_incorr+block{1}.voiced_noresp;
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl1,chi2stat1,pval1] = crosstab(x1,x2)

n1 = block{1}.unvoiced_corr; N1 = block{1}.unvoiced_corr+block{1}.unvoiced_incorr+block{1}.unvoiced_noresp;
n2 = block{2}.unvoiced_corr; N2 = block{2}.unvoiced_corr+block{1}.unvoiced_incorr+block{1}.unvoiced_noresp;
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl2,chi2stat2,pval2] = crosstab(x1,x2)
title({['Comparing accuracy of voiced between no-EStim vs EStim, p-val: ' num2str(pval1)];...
       ['Comparing accuracy of unvoiced between no-EStim vs EStim, p-val: ' num2str(pval2)]})

print(fullfile(Sbj_Metadata.results,['2ndhalf_' curr_block '_block_comp_onlyclean.jpg']),'-djpeg','-r300')
% print(fullfile(Sbj_Metadata.results,[curr_block '_block_comp.jpg']),'-djpeg','-r300')

%% Comparing between blocks with percentages, NOT raw trial numbers
syllable_couples = {'pa','ba';'ta','da';'fa','va';'sa','za'};
unvoiced_syll = {'pa';'ta';'fa';'sa'};
voiced_syll = {'ba';'da';'va';'za'};
block=cell(1,2);
for i=1:2
    load(fullfile(Sbj_Metadata.iEEG_data,Sbj_Metadata.BlockLists{i},[Sbj_Metadata.BlockLists{i} '_info.mat']));
    info.events.Accr = num2cell(strcmp(info.events.SyllablePresented,info.events.Response));
    info.events.Accr(strcmp(info.events.Response,'')) = {'noresp'};
    
    consonant = cell(size(info.events,1),1);
    for t=1:size(info.events,1)
        if ismember(info.events.SyllablePresented{t},voiced_syll)
            consonant{t} = 'voiced';
        else
            consonant{t} = 'unvoiced';
        end
        if strcmp(info.events.SyllablePresented{t},info.events.Response{t})
            info.events.Accr{t} = num2str(strcmp(info.events.SyllablePresented{t},info.events.Response{t}));
        elseif strcmp(info.events.Response{t},'')
            info.events.Accr{t} = 'noresp';
        else
            info.events.Accr{t} = '0';
        end
    end
    info.events = [info.events,cell2table(consonant)];
    % remove bad trials in e-stim block (the ones that extend into speaking or the ones that are not stimulating during the stim presentation)
    if strcmp(Sbj_Metadata.sbj_ID,'NS163') && strcmp(Sbj_Metadata.BlockLists{i},'SyllablePerception_EStim')
        info.events = info.events(~ismember(1:64,[1,11,14,16,21,22,25,31,32,33,37,41,43,51,54,58,61,62,63,64]),:);
    end
    if strcmp(Sbj_Metadata.sbj_ID,'NS163') && strcmp(Sbj_Metadata.BlockLists{i},'SyllablePerception')
        info.events = info.events(~ismember(1:92,1:46),:);
    end
    
    voiced_corr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'1'));
    voiced_incorr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'0'));
    voiced_noresp = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'noresp'));
    
    unvoiced_corr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'1'));
    unvoiced_incorr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'0'));
    unvoiced_noresp = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'noresp'));
    
    voiced_total = voiced_corr+voiced_incorr+voiced_noresp;
    unvoiced_total = unvoiced_corr+unvoiced_incorr+unvoiced_noresp;
    block{i}.voiced_corr = voiced_corr;
    block{i}.voiced_incorr = voiced_incorr;
    block{i}.voiced_noresp = voiced_noresp;
    block{i}.voiced_total = voiced_total;
    block{i}.unvoiced_corr = unvoiced_corr;
    block{i}.unvoiced_incorr = unvoiced_incorr;
    block{i}.unvoiced_noresp = unvoiced_noresp;
    block{i}.unvoiced_total = unvoiced_total;
end


figure('Position',[100 100 1200 1000])
to_bar = [block{1}.voiced_corr/block{1}.voiced_total block{1}.voiced_incorr/block{1}.voiced_total block{1}.voiced_noresp/block{1}.voiced_total;...
          block{2}.voiced_corr/block{2}.voiced_total block{2}.voiced_incorr/block{2}.voiced_total block{2}.voiced_noresp/block{2}.voiced_total;...
          block{1}.unvoiced_corr/block{1}.unvoiced_total block{1}.unvoiced_incorr/block{1}.unvoiced_total block{1}.unvoiced_noresp/block{1}.unvoiced_total;...
          block{2}.unvoiced_corr/block{2}.unvoiced_total block{2}.unvoiced_incorr/block{2}.unvoiced_total block{2}.unvoiced_noresp/block{2}.unvoiced_total]*100;
hBar = bar(to_bar,'stacked');
set(gca,'xticklabels',{'Voiced-NoEStim','Voiced-EStim','Unvoiced-NoEStim','Unvoiced-EStim'},'FontSize',15)
grid on
legend('Correct response','Incorrect Response','No-response','Location','southoutside','NumColumns',3)
ctr=[];ydt=[];
for k1 = 1:size(to_bar,2)
    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, hBar(k1).XOffset');    % Note: ‘XOffset’ Is An Undocumented Feature, This Selects The ‘bar’ Centres
    ydt(k1,:) = hBar(k1).YData;                                     % Individual Bar Heights
    if k1==1
        text(ctr(k1,:), ydt(k1,:), sprintfc('%.0f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',18, 'Color','k')
    else
        text(ctr(k1,:), ydt(k1,:)+sum(ydt(1:k1-1,:),1), sprintfc('%.0f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',18, 'Color','k')
    end
end

% Observed data
n1 = block{1}.voiced_corr; N1 = block{1}.voiced_corr+block{1}.voiced_incorr+block{1}.voiced_noresp;
n2 = block{2}.voiced_corr; N2 = block{2}.voiced_corr+block{1}.voiced_incorr+block{1}.voiced_noresp;
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl1,chi2stat1,pval1] = crosstab(x1,x2)

n1 = block{1}.unvoiced_corr; N1 = block{1}.unvoiced_corr+block{1}.unvoiced_incorr+block{1}.unvoiced_noresp;
n2 = block{2}.unvoiced_corr; N2 = block{2}.unvoiced_corr+block{1}.unvoiced_incorr+block{1}.unvoiced_noresp;
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl2,chi2stat2,pval2] = crosstab(x1,x2)
title({['Comparing accuracy of voiced between no-EStim vs EStim, p-val: ' num2str(pval1)];...
       ['Comparing accuracy of unvoiced between no-EStim vs EStim, p-val: ' num2str(pval2)]})

print(fullfile(Sbj_Metadata.results,['2ndhalf_' curr_block '_block_comp_onlyclean_percent.jpg']),'-djpeg','-r300')
% print(fullfile(Sbj_Metadata.results,[curr_block '_block_comp.jpg']),'-djpeg','-r300')
