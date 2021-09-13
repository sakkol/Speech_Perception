function SP_beh_analysis(Sbj_Metadata,curr_block)
% In order to create behavioral data results quickly

%% Administrative part
% Load response table
% response_table = readtable(fullfile(Sbj_Metadata.behavioral_root, curr_block, [curr_block '_response_table.xlsx']));
response_cell = readcell(fullfile(Sbj_Metadata.behavioral_root, curr_block, [curr_block '_response_table.xlsx']));

opts = detectImportOptions(fullfile(Sbj_Metadata.behavioral_root, curr_block, [curr_block '_response_table.xlsx']));
opts = setvartype(opts, response_cell(1,4:8), 'string');
response_table = readtable(fullfile(Sbj_Metadata.behavioral_root, curr_block, [curr_block '_response_table.xlsx']), opts);

% convert code column into double
if iscell(response_table.Condition_Code)
    for i=1:size(response_table,1)
        coll(i,1) = str2double(cell2mat(response_table.Condition_Code(i)));
    end
    response_table.Condition_Code = coll;
end

% convert sub/verb/numb/adj/noun column into string if needed
for c=4:8
    if isnumeric(response_table.(response_table.Properties.VariableNames{c}))
        response_table.(response_table.Properties.VariableNames{c}) = cellstr(num2str(response_table.(response_table.Properties.VariableNames{c})));
    end
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

select_cols = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560];

% Informations
CurrBlockInfo = makeBlockInfo(Sbj_Metadata,curr_block);
stim_type = strrep(CurrBlockInfo.Stim_type{:},'_',' ');

% Create folder to save
res_dir = fullfile(Sbj_Metadata.results,'behavioral');
if ~exist(res_dir,'dir')
    mkdir(res_dir)
end

%% Conditions in bar graph
figure('Position',[0 0 1000 1000])
bar([cond_corrects{:}]);
set(gca,'xticklabels',all_cond_names(:),'FontSize',11)
title({'Correct responses for each condition (/100)';['Stimulation type: ' stim_type]},'FontSize',15)
ylim([0 100])
grid on
text(1:length([cond_corrects{:}]),[cond_corrects{:}],num2str([cond_corrects{:}]'),'vert','bottom','horiz','center','Fontsize',25);
print(fullfile(res_dir,[curr_block '_conditions.jpg']),'-djpeg','-r300')

%% Correct responses for each word location
figure('Position',[100 100 1000 1000])
to_bar = [subj_corr{:};...
          verb_corr{:};...
          num_corr{:};...
          adj_corr{:};...
          noun_corr{:}];
hBar = bar(to_bar);

ylim([0 20])
set(gca,'xticklabels',{'Subject','Verb','Number','Adjective','Noun'},'FontSize',15)
title({'Correct responses for each word position (/20)';['Stimulation type: ' stim_type]})
grid on
legend(all_cond_names{:})                                                  % Return ‘bar’ Handle
for k1 = 1:size(to_bar,2)
    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, hBar(k1).XOffset');    % Note: ‘XOffset’ Is An Undocumented Feature, This Selects The ‘bar’ Centres
    ydt(k1,:) = hBar(k1).YData;                                     % Individual Bar Heights
    text(ctr(k1,:), ydt(k1,:), sprintfc('%.0f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'FontSize',18, 'Color','k')
end
print(fullfile(res_dir,[curr_block '_words.jpg']),'-djpeg','-r300')

%% Change as block goes
sorted_response_table = sortrows(response_table,'Order');

figure('Position',[200 200 1000 1000])
subplot(211)
plot(sorted_response_table.Acc_word_count,'k','LineWidth',3)
hold on
for c=1:length(all_cond_codes)
    curr_plot(:,2) = sorted_response_table.Acc_word_count(sorted_response_table.Condition_Code == all_cond_codes(c));
    curr_plot(:,1) = sorted_response_table.Order(sorted_response_table.Condition_Code == all_cond_codes(c));
    plot(curr_plot(:,1),curr_plot(:,2),'o','Color',select_cols(c,:),'LineWidth',4)
end
fit_res = fit(sorted_response_table.Order,sorted_response_table.Acc_word_count, 'poly1');
h=plot(fit_res);h.Color = [0 0 0];
h.LineWidth = .75;legend('hide');
xlim([0 length(sorted_response_table.Acc_word_count)+1])
yticks(0:1:5)
ylim([0 5.5])
xticks(0:5:length(sorted_response_table.Acc_word_count))
set(gca,'FontSize',13)
title('Correct responses as block continues','FontSize',15)
xlabel('');ylabel('');

subplot(212)
for c=1:length(all_cond_codes)
    curr_plot = sorted_response_table.Acc_word_count(sorted_response_table.Condition_Code == all_cond_codes(c));
    plot(curr_plot+(c-1)/15,'LineWidth',3,'Color',select_cols(c,:))
    hold on
    text(0.5,.005+c/15,[all_cond_names{c} ': ' num2str(cond_corrects{c}) '%'],'Color',select_cols(c,:),'FontSize',13,'FontWeight','bold','Units','normalized','HorizontalAlignment','center','BackgroundColor', 'white')
    
    fit_res = fit((1:20)',curr_plot, 'poly1');
    h=plot(fit_res);h.Color = select_cols(c,:);
    h.LineWidth = .75;legend('hide');
end
xlim([0 21])
xlabel('Trial #')
yticks(0:1:5)
ylim([0 5.5])
xticks(1:20)
set(gca,'FontSize',13)
title('Correct responses as block continues (separated for each condition)','FontSize',15)
xlabel('');ylabel('');
% legend(all_cond_names(:),'FontSize',13,'Location','best')
sgtitle(['Stimulation type: ' stim_type],'FontSize',15)
print(fullfile(res_dir,[curr_block '_change.jpg']),'-djpeg','-r300')


end