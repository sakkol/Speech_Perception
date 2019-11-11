function SP_beh_analysis(Sbj_Metadata,curr_block)
% In order to create behavioral data results (quicker)

% Load response table
response_table = readtable(fullfile(Sbj_Metadata.behavioral_root, curr_block, [curr_block '_response_table.xlsx']));

if iscell(response_table.Condition_Code)
    for i=1:size(response_table,1)
        coll(i,1) = str2double(cell2mat(response_table.Condition_Code(i)));
    end
    response_table.Condition_Code = coll;
end
% Read conditions:
all_cond_codes = sort(unique(response_table.Condition_Code));
for c=1:length(all_cond_codes)
    all_cond_names{c} = response_table.Condition_Name{response_table.Condition_Code==all_cond_codes(c)};
end

% collect correct response table for each condition
for c=1:length(all_cond_codes)
    cond_corrects_tables{c} = ismember(response_table{response_table.Condition_Code==all_cond_codes(c),4:8},'1');
    cond_corrects{c} = sum(sum(cond_corrects_tables{c}));
    subj_corr{c} = sum(ismember(response_table{response_table.Condition_Code==all_cond_codes(c),4},'1'));
    verb_corr{c} = sum(ismember(response_table{response_table.Condition_Code==all_cond_codes(c),5},'1'));
    num_corr{c} = sum(ismember(response_table{response_table.Condition_Code==all_cond_codes(c),6},'1'));
    adj_corr{c} = sum(ismember(response_table{response_table.Condition_Code==all_cond_codes(c),7},'1'));
    noun_corr{c} = sum(ismember(response_table{response_table.Condition_Code==all_cond_codes(c),8},'1'));
end

% Conditions in bar graph
figure('Position',[0 0 1000 1000])
bar([cond_corrects{:}]);
title('Correct responses for each condition (/100)')
set(gca,'xticklabels',all_cond_names(:),'FontSize',11)
ylim([0 100])
text(1:length([cond_corrects{:}]),[cond_corrects{:}],num2str([cond_corrects{:}]'),'vert','bottom','horiz','center','Fontsize',18);
print(fullfile(Sbj_Metadata.behavioral_root,curr_block,[curr_block '_conditions.jpg']),'-djpeg','-r300')

% Correct responses for each word
figure('Position',[100 100 1000 1000])
to_bar = [subj_corr{:};...
          verb_corr{:};...
          num_corr{:};...
          adj_corr{:};...
          noun_corr{:}];
hBar = bar(to_bar);

ylim([0 20])
set(gca,'xticklabels',{'Subject','Verb','Number','Adjective','Noun'},'FontSize',15)
title('Correct responses for each word (/20)')
grid on
legend(all_cond_names{:})                                                  % Return ‘bar’ Handle
for k1 = 1:size(to_bar,2)
    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, hBar(k1).XOffset');    % Note: ‘XOffset’ Is An Undocumented Feature, This Selects The ‘bar’ Centres
    ydt(k1,:) = hBar(k1).YData;                                     % Individual Bar Heights
    text(ctr(k1,:), ydt(k1,:), sprintfc('%.0f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'FontSize',10, 'Color','k')
end
print(fullfile(Sbj_Metadata.behavioral_root,curr_block,[curr_block '_words.jpg']),'-djpeg','-r300')

% Change as block goes
figure('Position',[200 200 1000 1000])
















end