function selections = select_conditions_GUI
% new select_conditions with a GUI

elecstim_cond_list ={'Control condition','Control condition','Control condition',...
                'HG - short delay','HG - long delay',...
                'STG - short delay','STG - long delay',...
                'Epicranial - short delay','Epicranial - long delay'...
                'HG - sinewave - inphase','HG - sinewave - outphase',...
                'STG - sinewave - inphase','STG - sinewave - outphase',...
                'Epicranial - sinewave - inphase','Epicranial - sinewave - outphase',...
                'Epicranial - 100Hz - nonmodulated'};

h = figure('Name', 'Select appropriate conditions', 'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'none',...
    'Units', 'Normalized', 'Position', [0.25, 0.4, 0.5, 0.35]);

% Column names and column format
columnname = {'Conditions','Language','Frequency','Left vs Right','Sentence vs Scrambled','Iso/A-chronous/Natural','Clean vs In-noise','Word per sentence','Electrical stim','E-stim delay'};
columnformat = {'logical',{'English','Spanish'},'numeric',{'L','R','both'},{'Sentence','Scrambled'},{'iso','a','natural'},{'clean','in-noise'},{'4-word','3-/5-word'},elecstim_cond_list,'numeric'}; %// Set the entries of the popup menu in a cell array. When the format is 'logical', the output in the table is a checkbox.

% Define the initial displayed data
d =    {true 'English',2.4,'L','Sentence','iso','clean','4-word','Control condition',0;
        true 'English',2.4,'L','Scrambled','iso','clean','4-word','Epicranial - short delay',50;
        true 'English',2.4,'L','Sentence','a','clean','4-word','Epicranial - long delay',150;
        true 'English',2.4,'L','Scrambled','a','clean','4-word','Control condition',0;
        false 'English',2.4,'L','Sentence','iso','in-noise','4-word','Control condition',0;
        false 'English',2.4,'L','Scrambled','iso','in-noise','4-word','Control condition',0;
        false 'English',2.4,'L','Sentence','a','in-noise','4-word','Control condition',0;
        false 'English',2.4,'L','Scrambled','a','in-noise','4-word','Control condition',0;
        false 'English',2.4,'L','Sentence','iso','in-noise','3-/5-word','Control condition',0;
        false 'English',2.4,'L','Scrambled','iso','in-noise','3-/5-word','Control condition',0;
        false 'English',2.4,'L','Sentence','a','in-noise','3-/5-word','Control condition',0;
        false 'English',2.4,'L','Scrambled','a','in-noise','3-/5-word','Control condition',0};

% Create the uitable
t = uitable(h,...
    'Units', 'Normalized', 'Position', [0.1, 0.1, 0.85, 0.8],...
    'Data', d,...
    'ColumnName', columnname,...
    'ColumnFormat', columnformat,...
    'ColumnEditable', [true true true true true true true true true true],... %// That's the important line. Entries set to true will allow you to create a popup menu for the whole column.
    'RowName','Options',...
    'DeleteFcn','CloseAndSave(gcbo)');

uiwait

% general problem is that when figure closes, the data is also cleaned. To
% work around this problem, 'CloseAndSave(gcbo)' in uitable creates a
% temporary file to load it here. After loading tmp file, delete it.
load('tmp.mat','selections_table');
delete tmp.mat
data=selections_table;

% selected ones are passed onto next step
selections = cell2table(data);
selections.Properties.VariableNames=replace(columnname,{'/',' ','-'},'_');
selections = selections(selections.Conditions==1,:);

end