function [control_events, control_ERP, control_wlt] = get_SP_controlwltERP(Sbj_Metadata,control_blocks,runagain)
% This function helps easily loading control conditions from within each
% block of run. 

%% Select blocks to import
vars=who;
if ~ismember(vars,'control_blocks')
    control_blocks = select_cont_blocks(Sbj_Metadata);
elseif isempty(control_blocks)
    control_blocks = select_cont_blocks(Sbj_Metadata);
end
if ~ismember(vars,'runagain')
    runagain = 0;
end
clear vars

%% bring in these blocks and combine only the control events
fprintf('These blocks are going to be used: %s\n',strjoin(control_blocks,', '))
% load info.mat file
load(fullfile(Sbj_Metadata.iEEG_data,Sbj_Metadata.BlockLists{1},[Sbj_Metadata.BlockLists{1} '_info.mat']),'info');

% check if this has already been run or need to be re-run
if ~exist(fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),[strjoin(control_blocks,'_') '_control_wltERP.mat']),'file') || runagain
    for b = 1:length(control_blocks)
        curr_block = control_blocks{b};
        fprintf('...loading:%s\n',curr_block)
        
        % Load iEEG
        load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_wlt.mat']))
        load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']))
        events = info.events;
        
        % Select only control events
        control_idx = events.Cond_code == 1;
        
        % Select ecog data
        cfg = [];
        cfg.trials = control_idx;
        [epoched_wlt] = ft_selectdata(cfg, epoched_wlt);
        curr_ERP = ft_selectdata(cfg, epoched_data);
        
        % Select events
        events = events(control_idx,:);
        
        % Append to overall list
        if b == 1
            control_wlt = epoched_wlt;
            control_events = events;
            control_ERP = curr_ERP;
        else
            cfg = [];
            cfg.parameter  = 'fourierspctrm';
            control_wlt = ft_appendfreq(cfg, control_wlt, epoched_wlt);
            cfg = [];
            control_ERP = ft_appendtimelock(cfg,control_ERP,curr_ERP);
            control_events = [control_events;events];
        end
        clear epoched_wlt events info curr_ERP control_idx epoched_data
    end
    
    control_wlt.cumtapcnt = ones([size(control_wlt.fourierspctrm,1),length(control_wlt.freq)]);
    
    fprintf('Saving to:\n\t->%s\n',fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),[strjoin(control_blocks,'_') '_control_wltERP.mat']))
    save(fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),[strjoin(control_blocks,'_') '_control_wltERP.mat']),'control_events','control_ERP','control_wlt','-v7.3')
else
    fprintf('Loading from:\n\t->%s\n',fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),[strjoin(control_blocks,'_') '_control_wltERP.mat']))
    load(fullfile(Sbj_Metadata.results, strjoin(control_blocks,'_'),[strjoin(control_blocks,'_') '_control_wltERP.mat']),'control_events','control_ERP','control_wlt')
    
end

end