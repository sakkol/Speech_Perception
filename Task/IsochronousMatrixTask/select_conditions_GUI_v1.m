function [selections, select_block] = select_conditions_GUI(select_block_input,language,LvsR)
% new select_conditions with a GUI

if ~exist('language','var') || isempty(language)
    language = ' ';
end
if ~exist('LvsR','var') || isempty(LvsR)
    LvsR = ' ';
end

%% Define the initially displayed data
tmp = load('default_conditions.mat');
default_table = ...
    {false,language,2.4,LvsR,'1attention-1sentence','30','Sentence','iso','clean','4-word','Control no-stim condition',0;
    false,language,2.4,LvsR,'1attention-1sentence','30','Scrambled','iso','clean','4-word','Epicranial - short delay',50;
    false,language,2.4,LvsR,'1attention-1sentence','30','Sentence','a','clean','4-word','Epicranial - long delay',150;
    false,language,2.4,LvsR,'1attention-1sentence','30','Scrambled','a','clean','4-word','Control no-stim condition',0;
    false,language,2.4,LvsR,'1attention-1sentence','30','Sentence','iso','in-noise','4-word','Control no-stim condition',0;
    false,language,2.4,LvsR,'1attention-1sentence','30','Scrambled','iso','in-noise','4-word','Control no-stim condition',0;
    false,language,2.4,LvsR,'1attention-1sentence','30','Sentence','a','in-noise','4-word','Control no-stim condition',0;
    false,language,2.4,LvsR,'1attention-1sentence','30','Scrambled','a','in-noise','4-word','Control no-stim condition',0;
    false,language,2.4,LvsR,'1attention-1sentence','30','Sentence','iso','in-noise','3-/5-word','Control no-stim condition',0;
    false,language,2.4,LvsR,'1attention-1sentence','30','Scrambled','iso','in-noise','3-/5-word','Control no-stim condition',0;
    false,language,2.4,LvsR,'1attention-1sentence','30','Sentence','a','in-noise','3-/5-word','Control no-stim condition',0;
    false,language,2.4,LvsR,'1attention-1sentence','30','Scrambled','a','in-noise','3-/5-word','Control no-stim condition',0};
if ~exist('select_block_input','var') || isempty(select_block_input)
    fig_title = 'PLEASE SELECT AND ARRANGE THE CONDITIONS!';
    d = default_table;
elseif ismember(select_block_input,fieldnames(tmp))
    tmp.(select_block_input).Language(:) = {language};
    tmp.(select_block_input).Left_vs_Right(:) = {LvsR};
    d = table2cell([tmp.(select_block_input);default_table]);
    fig_title = 'YOU MAY PLAY WITH PARAMETERS AS YOU WISH!';
end

elecstim_cond_list = {  'Control no-stim condition',...
    'HG - short delay','HG - long delay',...
    'STG - short delay','STG - long delay',...
    'Epicranial - short delay','Epicranial - long delay'...
    'HG - sinewave - inphase','HG - sinewave - outphase',...
    'STG - sinewave - inphase','STG - sinewave - outphase',...
    'Epicranial - sinewave - inphase','Epicranial - sinewave - outphase',...
    'Epicranial - 100Hz - nonmodulated'};

% Column names and column format
columnname = {'Activate Conditions','Language','Frequency','Left vs Right','Number of sentences','trial no per block','Sentence vs Scrambled','Iso/A-chronous/Natural','Clean vs In-noise','Word per sentence','Electrical stim','E-stim delay'};
columnformat = {'logical',{'English','Spanish'},'numeric',{'L','R','both'},{'1attention-1sentence','3sentences','5sentences'},{'10','15','25','30'},...
    {'Sentence','Scrambled'},{'iso','a','natural'},{'clean','in-noise'},{'4-word','3-/5-word'},elecstim_cond_list,'numeric'}; %// Set the entries of the popup menu in a cell array. When the format is 'logical', the output in the table is a checkbox.

%% Create the figure and uitable
h = figure('Name', fig_title, 'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'none',...
    'Units', 'Normalized', 'Position', [0.05, 0.2, 0.9, 0.4]);

t = uitable(h,...
    'Units', 'Normalized', 'Position', [0.05, 0.1, 0.9, 0.7],...
    'Data', d,...
    'ColumnName', columnname,...
    'ColumnFormat', columnformat,...
    'ColumnEditable', [true true true true true true true true true true true true],... %// That's the important line. Entries set to true will allow you to create a popup menu for the whole column.
    'ColumnWidth',{'auto' 'auto' 'auto' 'auto' 180 'auto' 'auto' 'auto' 'auto' 'auto' 250 'auto'},...
    'RowName',{'Conditions'},...
    'DeleteFcn','CloseAndSave(gcbo)');

%% User Interfaces to change different options in the whole table
% ui to select from default options
c0 = uicontrol(h,'Style','popupmenu',...
    'Units', 'Normalized','Position', [0.73 0.925 0.13 0.05],...
    'String', {'','iso_mat_24','FR_3sent','passive_5sent','clear_table'},...
    'Callback', {@selection_block,t},...
    'FontSize',14);
cnon = uicontrol(h,'Style','edit','Enable','off',...
    'Units', 'Normalized','Position', [0.57 0.92 0.16 0.05],...
    'String', 'Choose from default options:',...
    'FontSize',14);

% ui to select Left vs Right ear
c2 = uicontrol(h,'Style','popupmenu',...
    'Units', 'Normalized','Position', [0.45 0.925 0.04 0.05],...
    'String', {' ','L','R','both'},...
    'Callback', {@selection_LR,t},...
    'FontSize',14);
cnon = uicontrol(h,'Style','edit','Enable','off',...
    'Units', 'Normalized','Position', [0.3 0.92 0.15 0.05],...
    'String', 'Change left vs right ear:',...
    'FontSize',14);

% ui to change language if needed
c1 = uicontrol(h,'Style','popupmenu',...
    'Units', 'Normalized','Position', [0.12 0.85 0.05 0.045],...
    'String', {'English','Spanish'},...
    'Callback', {@selection_lang,t});
cnon = uicontrol(h,'Style','edit','Enable','off',...
    'Units', 'Normalized','Position', [0.05 0.85 0.07 0.045],...
    'String', 'Change language:');

% ui to change frequency of word presentations
c2 = uicontrol(h,'Style','popupmenu',...
    'Units', 'Normalized','Position', [0.29 0.85 0.03 0.045],...
    'String', {2,2.4,2.5,2.8,3,3.2},...
    'Callback', {@change_freq,t});
cnon = uicontrol(h,'Style','edit','Enable','off',...
    'Units', 'Normalized','Position', [0.22 0.85 0.07 0.045],...
    'String', 'Change frequency:');

% ui to change number of trials
c2 = uicontrol(h,'Style','popupmenu',...
    'Units', 'Normalized','Position', [0.44 0.85 0.03 0.045],...
    'String', {'10','15','25','30'},...
    'Callback', {@selection_notrial,t});
cnon = uicontrol(h,'Style','edit','Enable','off',...
    'Units', 'Normalized','Position', [0.35 0.85 0.09 0.045],...
    'String', 'Change no of trials:');

uiwait

%% After table is
% general problem is that when figure closes, the data is also cleaned. To
% work around this problem, 'CloseAndSave(gcbo)' in uitable creates a
% temporary file to load it here. After loading tmp file, delete it.
load('tmp.mat','data_to_save');
delete tmp.mat
if exist('tmp2.mat','file')
    load('tmp2.mat','select_block');
    delete tmp2.mat
else
    select_block=select_block_input;
end

% selected ones are passed onto next step
selections = cell2table(data_to_save);
selections.Properties.VariableNames=replace(columnname,{'/',' ','-'},'_');
selections = selections(selections.Activate_Conditions==1,:);

% check the output if anything missing
if isempty(selections)
    error('No conditions were selected, change "Activate Conditions" column boxes for the ones you want to involve!')
end
if any(isempty(selections.Language)) || any(ismember(selections.Language,{' ',''}))
    error('Language is left blank, please correct it!')
end
if any(isempty(selections.Left_vs_Right)) || any(ismember(selections.Left_vs_Right,{' ',''}))
    error('Left_vs_Right is left blank, please correct it!')
end

end


%% Subfunctions of User Interfaces to change options in whole table
function selection_block(src,event,t)
tmp = load('default_conditions.mat');
srcstr = get(src,'String');
srcval = get(src,'Value');
select_block = srcstr{srcval};

switch select_block
    case 'iso_mat_24'
        selected_block = table2cell(tmp.iso_mat_24);
    case 'FR_3sent'
        selected_block = table2cell(tmp.FR_3sent);
    case 'passive_5sent'
        selected_block = table2cell(tmp.passive_5sent);
    case 'clear_table'
        selected_block = table2cell(tmp.clear_table);
    otherwise
        selected_block = cell(0,11);
end

if strcmp(select_block,'clear table')
    tdata = cell(0,11);
else
    tdata = get(t,'Data');
    selected_block(:,2) = tdata(1,2);
    selected_block(:,4) = tdata(1,4);
end

% save block selection for later use
save('tmp2.mat','select_block')

tdata = [selected_block;tdata];
set(t,'Data',tdata);
% disp(['Selection: ' srcstr{srcval}]);
end


function selection_lang(src,event,t)
srcstr = get(src,'String');
srcval = get(src,'Value');
tdata = get(t,'Data');
tdata(:,2) = srcstr(srcval);
set(t,'Data',tdata);
% disp(['Selection: ' srcstr{srcval}]);
end

function selection_LR(src,event,t)
srcstr = get(src,'String');
srcval = get(src,'Value');
tdata = get(t,'Data');
tdata(:,4) = srcstr(srcval);
set(t,'Data',tdata);
% disp(['Selection: ' srcstr{srcval}]);
end

function selection_notrial(src,event,t)
srcstr = get(src,'String');
srcval = get(src,'Value');
tdata = get(t,'Data');
tdata(:,6) = srcstr(srcval);
set(t,'Data',tdata);
% disp(['Selection: ' srcstr{srcval}]);
end

function change_freq(src,event,t)
srcstr = get(src,'String');
srcval = get(src,'Value');
tdata = get(t,'Data');
tdata(:,3) = {str2double(srcstr{srcval})};
set(t,'Data',tdata);
% disp(['Selection: ' srcstr{srcval}]);
end