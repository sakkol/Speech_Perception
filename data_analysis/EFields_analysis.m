function EFields_analysis(Sbj_Metadata,curr_block,aroundPeak)
vars=who; % around peak is not ready actually
if ~any(ismember(vars,'aroundPeak')) || strcmp(Sbj_Metadata.project_name,'Efields_Alone')
    aroundPeak = 0; 
end
clear vars

%% Load bipolar
fprintf('Loading bipolar data from: %s\n',fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_bp.mat']))
tmp = load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_ecog_bp.mat']));ecog_bp = tmp.ecog_bp;clear tmp
% Load events
load(fullfile(Sbj_Metadata.iEEG_data, curr_block,[curr_block '_info.mat']),'info')
events = info.events;

%% Generate trial structure
if strcmp(Sbj_Metadata.project_name,'Efields_Alone')
    trl = [];
    if ~aroundPeak
        trl(:,1) = events.eventsOnsetSec;
        trl(:,2) = events.eventsEndSec;
        trl(:,3) = 0;
        to_print_name = '';
    else % this is not ready for Efields_Alone
        %     [~,l] = max(events.Stimuli{noncontrol_idx(el)}(:,2));
        %     trl(:,1) = events.trial_onsets(noncontrol_idx(el)) + l/24000 -0.02; % 20ms around the peak e-stim
        %     trl(:,2) = events.trial_onsets(noncontrol_idx(el)) + l/24000 +0.02;
        %     trl(:,3) = 0;
        %     to_print_name = '_aroundPeak';
    end
elseif strcmp(Sbj_Metadata.project_name,'Speech_Perception')
    % Select events with stimulation
    noncontrol_idx = find(events.Cond_code ~= 1);
    
    % Generate trial structure
    trl = [];
    for i = 1:length(noncontrol_idx)
        if ~aroundPeak
            trl(i,1) = events.speech_onsets(noncontrol_idx(i)) + events.trial_details{noncontrol_idx(i)}.delay;
            trl(i,2) = trl(i,1) + events.word_info{noncontrol_idx(i)}.offset(end);
            trl(i,3) = 0;
            to_print_name = '';
        else
            [~,l] = max(events.Stimuli{noncontrol_idx(i)}(:,2));
            trl(i,1) = events.trial_onsets(noncontrol_idx(i)) + l/24000 -0.02; % 20ms around the peak e-stim
            trl(i,2) = events.trial_onsets(noncontrol_idx(i)) + l/24000 +0.02;
            trl(i,3) = 0;
            to_print_name = '_aroundPeak';
        end
    end
end

% Put events with different lengths
trl = floor(trl*ecog_bp.ftrip.fsample);

%% Remove pairs that include artefact channels (both artifact_patient and artifact_block from info.channelinfo)
el_id = zeros(length(ecog_bp.ftrip.label),1);
for el = 1:length(ecog_bp.ftrip.label)
    if ~any(ismember(strsplit(ecog_bp.ftrip.label{el},'-'),info.channelinfo.Label(info.channelinfo.artifact_patient==1 | info.channelinfo.artifact_block==1)))
        el_id(el,1) = 1;
    else
        el_id(el,1) = 0;
    end
end

% Epoch without bad channels and detrend
cfg         = [];
cfg.channel = ecog_bp.ftrip.label(logical(el_id));
cfg.detrend = 'yes';
tmp         = ft_preprocessing(cfg,ecog_bp.ftrip);

cfg         = [];
cfg.trl     = trl;
stim_epochs = ft_redefinetrial(cfg,tmp);
clear tmp

%% Find the locations of electrodes
elecInfoFname=fullfile(Sbj_Metadata.freesurfer,'elec_recon',[Sbj_Metadata.fsname '.electrodeNames']);
elecInfo=csv2Cell(elecInfoFname,' ',2);
elec_toplot_names = elecInfo(:,1);

% Load electrode coordinates from PIAL
elecCoordFname=fullfile(Sbj_Metadata.freesurfer,'elec_recon',[Sbj_Metadata.fsname '.LEPTO']);
elecCoordCsv=csv2Cell(elecCoordFname,' ',2);
fullElecCoord=zeros(size(elecCoordCsv));
% Convert coordinates from string to #
for a=1:size(elecCoordCsv,1),
    for b=1:3,
        fullElecCoord(a,b)=str2double(elecCoordCsv{a,b});
    end
end
% check if lengths are correct
if length(fullElecCoord) ~= length(elec_toplot_names)
    error('There is a problem with number of elec names and elec coords, please load and check!\n')
end

%% Find distances based on bipolar pairs
distances = zeros(length(stim_epochs.label),1);
labelss = [];labelss{length(stim_epochs.label),2} = ''; % labelss is the splitted version of bipolar namings
for el = 1:length(stim_epochs.label)
    if length(strsplit(stim_epochs.label{el},'-')) ~= 2
        error('There is a problem with label naming, please load and check!\n')
    else
        labelss(el,1:2) = strsplit(stim_epochs.label{el},'-');
    end
    for l = 1:2
        loc = zeros(2,3);
        if sum(contains(elec_toplot_names,labelss{el,l}))
            loc(l,1:3) = fullElecCoord(strcmpi(labelss{el,l},elec_toplot_names),:);
        else
            loc(l,1:3) = [0 0 0];
        end
    end
    distances(el,1) = norm(loc(2,:) -  loc(1,:));
end

%% Calculate electric fields
efields_mean = zeros(length(stim_epochs.label),1);
efields_sterr = zeros(length(stim_epochs.label),1);
for el = 1:length(stim_epochs.label) % Loop electrode pairs
    % collect max of the absolute voltage in each trial for every bipolar pair
    tmp = zeros(length(stim_epochs.trial),1);
    for t = 1:length(stim_epochs.trial) % loop trials
        tmp(t,1) = max(abs(stim_epochs.trial{t}(el,:)));
    end
    % Calculate average across trials
    efields_mean(el,1) = mean(tmp);
    efields_sterr(el,1) = stderr(tmp);
end

%% Put colors of stim elecs differently
stim_elecs = strsplit(ecog_bp.params.CurrBlockInfo.StimCh_anode{1},',');
stim_elecs = [stim_elecs strsplit(ecog_bp.params.CurrBlockInfo.StimCh_cathode{1},',')];
for el = 1:length(elec_toplot_names)
    if ~strcmp(elec_toplot_names{el},stim_elecs)
        elec_col(el,:) = [.8 0.8 .8];
    elseif any(strcmp(elec_toplot_names{el},strsplit(ecog_bp.params.CurrBlockInfo.StimCh_anode{1},',')))
        elec_col(el,:) = [1 0 0];
    elseif any(strcmp(elec_toplot_names{el},strsplit(ecog_bp.params.CurrBlockInfo.StimCh_cathode{1},',')))
        elec_col(el,:) = [0 0 1];
    else
        error(sprintf('Couldn''t assign elec color to this electrode: %s\n',elec_toplot_names{el}))
    end
end

%% Prepare pairs cell for plotPialSurf input
pairs=[];pairs{1,4}='';
efields_mean_toplot = 0;
for el = 1:length(stim_epochs.label) % Loop electrode pairs
    % distances with zero are nonexistent electrodes AND nonstimulated
    % electrodes
    if distances(el) ~= 0 && ~any(strcmp(labelss(el,1),stim_elecs)) && ~any(strcmp(labelss(el,2),stim_elecs))
        efields_mean_toplot(end+1,1) = efields_mean(el);
        pairs(end+1,1:2) = labelss(el,:);
        pairs{end,4} = elecInfo{strcmp(labelss{el,1},elec_toplot_names),3};
    end
end
pairs(1,:) = [];efields_mean_toplot(1) = [];


% bwr = load('bwr_cmap.mat');
% p=bwr.rgb_vals;
p=hot(256);
x=efields_mean_toplot; %data to be plotted
ran=range(x); %finding range of data
min_val=min(x);%finding maximum value of data
y=floor(((x-min_val)/ran)*length(p))+1;
y(y==257)=256;
% col=zeros(length(efields_mean_toplot),3);

for el=1:length(efields_mean_toplot)
  a=y(el);
%   col(i,:)=p(a,:);
  pairs{el,3}=p(a,:);
%   stem3(i,i,x(i),'Color',col(i,:))
%   hold on
end

%% Plot
cfg=[];
cfg.ignoreDepthElec = 'n';
cfg.elecCoord = 'LEPTO';
cfg.elecNames = elec_toplot_names;
cfg.elecColors = elec_col;
cfg.elecCbar = 'n';
cfg.pullOut = 0;
cfg.fsurfSubDir = erase(Sbj_Metadata.freesurfer,Sbj_Metadata.fsname);
if any(strcmpi(pairs(:,4),'R')) && any(strcmpi(pairs(:,4),'L'))
    cfg.view='omni';
elseif any(strcmpi(pairs(:,4),'L')) && ~any(strcmpi(pairs(:,4),'R'))
    cfg.view='lomni';
elseif ~any(strcmpi(pairs(:,4),'L')) && any(strcmpi(pairs(:,4),'R'))
    cfg.view='romni';
end
cfg.pairs=pairs;
cfg.showLabels='n';
if any(strcmp(elecInfo(:,2),'D')) && ~all(contains(elec_toplot_names(strcmp(elecInfo(:,2),'D')),'ref','IgnoreCase',1))
    cfg.opaqueness = 0.4;
else
    cfg.opaqueness = 0.9;
end
cfg.elecShape = 'sphere';
cfg.elecSize = .7;
cfg.lineWidth = 5;
if strcmp(Sbj_Metadata.project_name,'Efields_Alone')
    cfg.title = [Sbj_Metadata.sbj_ID '; Stimulation type: ' ecog_bp.params.CurrBlockInfo.pulseType{1}...
        '; Location: ' ecog_bp.params.CurrBlockInfo.Location{1}...
        '; Max stim amp: ' ecog_bp.params.CurrBlockInfo.Amp_mA{1} '; Anode(red):' ecog_bp.params.CurrBlockInfo.StimCh_anode{1} '; Cathode(blue):' ecog_bp.params.CurrBlockInfo.StimCh_cathode{1}];
elseif strcmp(Sbj_Metadata.project_name,'Speech_Perception')
    cfg.title = [Sbj_Metadata.sbj_ID '; Stimulation type: ' ecog_bp.params.CurrBlockInfo.Stim_type{1}...
        '; Location: ' ecog_bp.params.CurrBlockInfo.Location{1}...
        '; Max stim amp: ' ecog_bp.params.CurrBlockInfo.max_Amp{1} '; Anode(red):' ecog_bp.params.CurrBlockInfo.StimCh_anode{1} '; Cathode(blue):' ecog_bp.params.CurrBlockInfo.StimCh_cathode{1}];
end
plotPialSurf(Sbj_Metadata.fsname,cfg);

% Create and delete new axes to plot colorbar
ax = axes;
colormap(p);
cmaph = colorbar(ax);
cmaph.Ticks = linspace(0,1,8);
cmaph.TickLabels = num2cell(linspace(min(efields_mean_toplot), max(efields_mean_toplot),8));
cmaph.FontSize = 13;
cmaph.LineWidth = 1.5;
colorTitleHandle = get(cmaph,'Title');
set(colorTitleHandle ,'String','Electric field (V/m)','FontSize',15);
a=get(cmaph); %gets properties of colorbar
a = a.Position; %gets the positon and size of the color bar
set(cmaph,'Position',[a(1)+0.05 a(2) 0.03 0.8])% To change size
ax.Visible = 'off';
set(gcf,'Units','normalized');
set(gcf,'Position', [0 0 1 1]);

%% Save plot and data
to_print_folder = fullfile(Sbj_Metadata.results,'Efields',curr_block);
if ~exist(to_print_folder,'dir'),mkdir(to_print_folder),end
fprintf('\n\tSaving data and plot to: %s\n',to_print_folder)
save(fullfile(to_print_folder,[curr_block, '_efield_data' to_print_name '.mat']),'efields_mean','efields_sterr','pairs')
print(fullfile(to_print_folder,[curr_block, '_efield_plot' to_print_name '.jpg']),'-djpeg','-r300')

end