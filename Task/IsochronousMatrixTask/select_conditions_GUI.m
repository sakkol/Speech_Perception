function [selections, select_block] = select_conditions_GUI(default_selection,language,LvsR)
% new select_conditions with a GUI

if ~exist('language','var') || isempty(language)
    language = ' ';
end
if ~exist('LvsR','var') || isempty(LvsR)
    LvsR = ' ';
end

% Define the initial displayed data
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
if ~exist('default_selection','var') || isempty(default_selection)
    fig_title = 'PLEASE SELECT AND ARRANGE THE CONDITIONS!';
    d = default_table;
elseif ismember(default_selection,fieldnames(tmp))
    d = table2cell([tmp.(default_selection);default_table]);
    fig_title = 'YOU MAY PLAY WITH PARAMETERS AS YOU WISH!';
else
    
end

elecstim_cond_list = {  'Control no-stim condition',...
                        'HG - short delay','HG - long delay',...
                        'STG - short delay','STG - long delay',...
                        'Epicranial - short delay','Epicranial - long delay'...
                        'HG - sinewave - inphase','HG - sinewave - outphase',...
                        'STG - sinewave - inphase','STG - sinewave - outphase',...
                        'Epicranial - sinewave - inphase','Epicranial - sinewave - outphase',...
                        'Epicranial - 100Hz - nonmodulated'};

h = figure('Name', fig_title, 'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'none',...
    'Units', 'Normalized', 'Position', [0.15, 0.3, 0.7, 0.35]);

% Column names and column format
columnname = {'Activate','Language','Frequency','Left vs Right','Number of sentences','trial no per block','Sentence vs Scrambled','Iso/A-chronous/Natural','Clean vs In-noise','Word per sentence','Electrical stim','E-stim delay'};
columnformat = {'logical',{'English','Spanish'},'numeric',{'L','R','both'},{'1attention-1sentence','3sentences','5sentences'},{'10','15','25','30'},...
    {'Sentence','Scrambled'},{'iso','a','natural'},{'clean','in-noise'},{'4-word','3-/5-word'},elecstim_cond_list,'numeric'}; %// Set the entries of the popup menu in a cell array. When the format is 'logical', the output in the table is a checkbox.

% Create the uitable
t = uitable(h,...
    'Units', 'Normalized', 'Position', [0.05, 0.1, 0.9, 0.65],...
    'Data', d,...
    'ColumnName', columnname,...
    'ColumnFormat', columnformat,...
    'ColumnEditable', [true true true true true true true true true true true true],... %// That's the important line. Entries set to true will allow you to create a popup menu for the whole column.
    'ColumnWidth',{'auto' 'auto' 'auto' 'auto' 180 'auto' 'auto' 'auto' 'auto' 'auto' 250 'auto'},...
    'RowName',{'Options'},...
    'DeleteFcn','CloseAndSave(gcbo)');

c0 = uicontrol(h,'Style','popupmenu',...
    'Units', 'Normalized','Position', [0.5 0.925 0.13 0.05],...
    'String', {'','2.4Hz-isochronous-matrix','3sentence-free-recall','5sentence-passive','clear table'},...
    'Callback', {@selection_block,t},...
    'FontSize',14);
cnon = uicontrol(h,'Style','edit','Enable','off',...
    'Units', 'Normalized','Position', [0.37 0.92 0.13 0.05],...
    'String', 'Select block option:',...
    'FontSize',14);

c1 = uicontrol(h,'Style','popupmenu',...
    'Units', 'Normalized','Position', [0.12 0.85 0.05 0.045],...
    'String', {'English','Spanish'},...
    'Callback', {@selection_lang,t});
cnon = uicontrol(h,'Style','edit','Enable','off',...
    'Units', 'Normalized','Position', [0.05 0.85 0.07 0.045],...
    'String', 'Change language:');
    
c2 = uicontrol(h,'Style','popupmenu',...
    'Units', 'Normalized','Position', [0.29 0.85 0.03 0.045],...
    'String', {2,2.4,2.5,2.8,3,3.2},...
    'Callback', {@change_freq,t});
cnon = uicontrol(h,'Style','edit','Enable','off',...
    'Units', 'Normalized','Position', [0.22 0.85 0.07 0.045],...
    'String', 'Change frequency:');

c2 = uicontrol(h,'Style','popupmenu',...
    'Units', 'Normalized','Position', [0.44 0.85 0.04 0.045],...
    'String', {' ','L','R','both'},...
    'Callback', {@selection_LR,t});
cnon = uicontrol(h,'Style','edit','Enable','off',...
    'Units', 'Normalized','Position', [0.35 0.85 0.09 0.045],...
    'String', 'Change left vs right ear:');
    
c2 = uicontrol(h,'Style','popupmenu',...
    'Units', 'Normalized','Position', [0.64 0.85 0.03 0.045],...
    'String', {'10','15','25','30'},...
    'Callback', {@selection_notrial,t});
cnon = uicontrol(h,'Style','edit','Enable','off',...
    'Units', 'Normalized','Position', [0.55 0.85 0.09 0.045],...
    'String', 'Change no of trials:');
    
uiwait

% general problem is that when figure closes, the data is also cleaned. To
% work around this problem, 'CloseAndSave(gcbo)' in uitable creates a
% temporary file to load it here. After loading tmp file, delete it.
load('tmp.mat','data_to_save');
delete tmp.mat
if exist('tmp2.mat','file')
    load('tmp2.mat','select_block');
    delete tmp2.mat
else
    select_block=[];
end

% selected ones are passed onto next step
selections = cell2table(data_to_save);
selections.Properties.VariableNames=replace(columnname,{'/',' ','-'},'_');
selections = selections(selections.Activate==1,:);

% check the output if anything missing
if isempty(selections)
    error('No conditions were selected, change Activate column boxes for the ones you want to involve!')
end
if any(isempty(selections.Language)) || any(strcmp(selections.Language,' '))
    error('Language is left blank, please correct it!')
end
if any(isempty(selections.Left_vs_Right)) || any(strcmp(selections.Left_vs_Right,' '))
    error('Left_vs_Right is left blank, please correct it!')
end

end


%% Subfunctions
function selection_block(src,event,t)
tmp = load('default_conditions.mat');
srcstr = get(src,'String');
srcval = get(src,'Value');
select_block = srcstr{srcval};

switch select_block
    case '2.4Hz-isochronous-matrix'
        selected_block = table2cell(tmp.iso_mat_24);
    case '3sentence-free-recall'
        selected_block = table2cell(tmp.FR_3sent);
    case '5sentence-passive'
        selected_block = table2cell(tmp.passive_5sent);
    case 'clear table'
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