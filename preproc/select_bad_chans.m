% select_bad_chans
% To be able to easily select the bad electrodes and a handful others from a nice table

soz_id = ismember(ecog.ftrip.label,ecog.szr_onset_chans);
spike_id = ismember(ecog.ftrip.label,ecog.spike_chans);
bad_id = ismember(ecog.ftrip.label,ecog.bad_chans);
vars=who;
if ~any(ismember(vars,'selected_chans'))
    select_id = false(length(ecog.ftrip.label),1);
    selected_chans={};
else
    select_id = ismember(ecog.ftrip.label,selected_chans);
end

d = table(ecog.ftrip.label,bad_id,spike_id,soz_id,select_id);
% newData = select_bad_chan_GUI(d);
% function data = MyDeleteFcn(t)
% data = get(t,'Data');
% end
h = figure('Name', 'Check appropriate boxes', 'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'none',...
    'Units', 'Normalized', 'Position', [0.4, 0.15, 0.2, 0.8]);


uit = uitable('Parent',h,...
    'Units', 'Normalized', 'Position', [0.1, 0.1, 0.85, 0.8],...
    'Data', table2cell(d),...
    'ColumnName',{'Label','bad chans','spikey','SOZ','Selected'},...
    'ColumnEditable',true,...
     'DeleteFcn','newData = MyDeleteFcn(gcbo);');


% % set(uit,'CellEditCallback','newData = get(h,''Data'');');
% % btn = uicontrol('Style', 'pushbutton', 'String', 'Update',...
% %         'Position', [420 480 100 40],...
% %         'Callback', 'newData = get(uit,''Data'');'); 
waitfor(h)
% % newData = get(uit,'Data');


% Vocalize the results
already_nonbad = ecog.ftrip.label(~ismember(ecog.ftrip.label,ecog.bad_chans));
new_bads = already_nonbad(ismember(already_nonbad,ecog.ftrip.label([newData{:,2}])));

already_nonsp = ecog.ftrip.label(~ismember(ecog.ftrip.label,ecog.spike_chans));
new_spikes = already_nonsp(ismember(already_nonsp,ecog.ftrip.label([newData{:,3}])));

already_nonSOZ = ecog.ftrip.label(~ismember(ecog.ftrip.label,ecog.szr_onset_chans));
new_SOZs = already_nonSOZ(ismember(already_nonSOZ,ecog.ftrip.label([newData{:,4}])));

already_selected = ecog.ftrip.label(~ismember(ecog.ftrip.label,selected_chans));
new_selected_chans = already_selected(ismember(already_selected,ecog.ftrip.label([newData{:,5}])));

fprintf('Newly assigned bad channels: %s\n',strjoin(new_bads,', '))
fprintf('Newly assigned spikey channels: %s\n',strjoin(new_spikes,', '))
fprintf('Newly assigned SOZ channels: %s\n',strjoin(new_SOZs,', '))
fprintf('Newly SELECTED channels: %s\n',strjoin(new_selected_chans,', '))

% Put it back
ecog.bad_chans = ecog.ftrip.label([newData{:,2}]);
ecog.spike_chans = ecog.ftrip.label([newData{:,3}]);
ecog.szr_onset_chans = ecog.ftrip.label([newData{:,4}]);
selected_chans = ecog.ftrip.label([newData{:,5}]);

clear new_bads new_spikes new_SOZs newData uit h soz_id bad_id spike_id select_id vars d already_nonbad already_nonsp already_nonSOZ already_selected new_selected_chans
