% select_conditions (based on select_bad_chans)
% To be able to easily select the conditions from a nice table if needed to select manually

h = figure('Name', 'Select appropriate conditions', 'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'none',...
    'Units', 'Normalized', 'Position', [0.25, 0.4, 0.5, 0.35]);

% Column names and column format
columnname = {'Conditions','Language','Frequency','Left vs Right','Sentence vs Scrambled','Iso/A-chronous/Natural','Clean vs In-noise','Word per sentence','Electrical stim'};
columnformat = {'logical',{'English','Spanish'},'numeric',{'L','R','both'},{'Sentence','Scrambled'},{'iso','a','natural'},{'clean','in-noise'},{'4-word','3-/5-word'},'logical'}; %// Set the entries of the popup menu in a cell array. When the format is 'logical', the output in the table is a checkbox.

% Define the initial displayed data
d =    {true 'English',2.4,'L','Sentence','iso','clean','4-word',false;
        true 'English',2.4,'L','Scrambled','iso','clean','4-word',false;
        true 'English',2.4,'L','Sentence','a','clean','4-word',false;
        true 'English',2.4,'L','Scrambled','a','clean','4-word',false;
        false 'English',2.4,'L','Sentence','iso','in-noise','4-word',false;
        false 'English',2.4,'L','Scrambled','iso','in-noise','4-word',false;
        false 'English',2.4,'L','Sentence','a','in-noise','4-word',false;
        false 'English',2.4,'L','Scrambled','a','in-noise','4-word',false;
        false 'English',2.4,'L','Sentence','iso','in-noise','3-/5-word',false;
        false 'English',2.4,'L','Scrambled','iso','in-noise','3-/5-word',false;
        false 'English',2.4,'L','Sentence','a','in-noise','3-/5-word',false;
        false 'English',2.4,'L','Scrambled','a','in-noise','3-/5-word',false};

% Create the uitable
t = uitable('Parent',h,...
    'Units', 'Normalized', 'Position', [0.1, 0.1, 0.85, 0.8],...
    'Data', d,...
    'ColumnName', columnname,...
    'ColumnFormat', columnformat,...
    'ColumnEditable', [true true true true true true true true true],... %// That's the important line. Entries set to true will allow you to create a popup menu for the whole column.
    'RowName','Options',...
    'DeleteFcn','newData = MyDeleteFcn(gcbo);');

waitfor(h)

selections = cell2table(newData);
selections.Properties.VariableNames=replace(columnname,{'/',' ','-'},'_');
selections = selections(selections.Conditions==1,:);
clear h t d columnname columnformat newData