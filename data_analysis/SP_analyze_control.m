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

for t = 1:size(control_events,1)
    
    for w = 1:5
        % Separate words
        % Check if correct
        % Collect fft in a cell structure?
        if strcmp(control_events.word_info{t}.response{w},'1')
            % find closest timepoint
            cl_times = [nearest(control_wlt.time, control_events.word_info{t}.onset(w))...
                nearest(control_wlt.time,control_events.word_info{t}.offset(w))];
            corr_rspn_fft{end+1,1} = squeeze(control_wlt.fourierspctrm(t,:,:,cl_times(1):cl_times(2)));
            
        elseif strcmp(control_events.word_info{t}.response{w},'0')
            % find closest timepoint
            cl_times = [nearest(control_wlt.time, control_events.word_info{t}.onset(w))...
                nearest(control_wlt.time,control_events.word_info{t}.offset(w))];
            no_rspn_fft{end+1,1} = squeeze(control_wlt.fourierspctrm(t,:,:,cl_times(1):cl_times(2)));
            
        else
            % find closest timepoint
            cl_times = [nearest(control_wlt.time, control_events.word_info{t}.onset(w))...
                nearest(control_wlt.time,control_events.word_info{t}.offset(w))];
            wrng_rspn_fft{end+1,1} = squeeze(control_wlt.fourierspctrm(t,:,:,cl_times(1):cl_times(2)));
            
        end
        
        
    end
    
    
end



