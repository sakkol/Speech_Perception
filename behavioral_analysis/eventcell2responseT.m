% function eventcell2responseT(Sbj_Metadata,curr_block)
% To quickly create the response table with empty cells. Optionally to
% write the table here in the function.

% load event_cell
load(fullfile(Sbj_Metadata.behavioral_root, curr_block, [curr_block '.mat']),'events_cell')

% create dummies
Order = [1:size(events_cell,1)]';
Sentence_Code = events_cell(:,1);
Sentence = events_cell(:,2);
Subject = cell(size(events_cell,1),1);
Verb = cell(size(events_cell,1),1);
Number = cell(size(events_cell,1),1);
Adjective = cell(size(events_cell,1),1);
Noun = cell(size(events_cell,1),1);
Acc_word_count = cell(size(events_cell,1),1);
Condition_Code = events_cell(:,3);
Condition_Name = events_cell(:,4);

newrestT = [array2table(Order),array2table(Sentence_Code),array2table(Sentence),...
    array2table(Subject),array2table(Verb),array2table(Number),array2table(Adjective),array2table(Noun),...
    array2table(Acc_word_count),array2table(Condition_Code),array2table(Condition_Name)];

%% Plot the table to be able to fill up
h = figure('Name', 'Check appropriate boxes', 'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'none',...
    'Units', 'Normalized', 'Position', [0.2, 0.15, 0.6, 0.5]);


uit = uitable('Parent',h,...
    'Units', 'Normalized', 'Position', [0.1, 0.1, 0.85, 0.8],...
    'Data', table2cell(newrestT),...
    'ColumnName',newrestT.Properties.VariableNames,...
    'ColumnEditable',true,...
     'DeleteFcn','newData = MyDeleteFcn(gcbo);');

waitfor(h)

%% Write the table if needed
% cell2table and display it
newDataTable = cell2table(newData,'VariableNames',newrestT.Properties.VariableNames)

answer = questdlg('Would you like to overwrite these to response table?', ...
    'Over-writing?', ...
	'Yes','No','Quit here!','Quit here!');
% Handle response
switch answer
    case 'Yes'
        disp([answer ' over-writing!'])
        
        writetable(newDataTable, fullfile(Sbj_Metadata.behavioral_root, curr_block, [curr_block '_response_table.xlsx']));
    case 'No'
        disp([answer ', ok, leaving here!'])
    case 'Quit here!'
        error('Why I''m here!!')
end

clear uit h events_cell Order Sentence_Code Sentence Subject Verb Number Adjective Noun Acc_word_count Condition_Code Condition_Name newDataTable newData answer newrestT
% end
