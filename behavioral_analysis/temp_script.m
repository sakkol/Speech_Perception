control_mean = mean(real1_response_table.Acc_word_count(strcmp(real1_response_table.Condition_Code,'1')));
control_std = std(real1_response_table.Acc_word_count(strcmp(real1_response_table.Condition_Code,'1')));

control_data = real1_response_table.Acc_word_count(strcmp(real1_response_table.Condition_Code,'1'));
STG_short_data = real1_response_table.Acc_word_count(strcmp(real1_response_table.Condition_Code,'4'));
STG_long_data = real1_response_table.Acc_word_count(strcmp(real1_response_table.Condition_Code,'5'));

figure
xlim([0.5 3.5])
hold on
boxplot([control_data; STG_short_data;STG_long_data],...
    [ones(length(control_data),1);2*ones(length(STG_short_data),1);3*ones(length(STG_long_data),1)],...
    'Widths',0.45,'Colors','b')

scatter(ones(length(control_data),1),control_data,'filled','jitter','on','MarkerFaceColor','k');
scatter(2*ones(length(STG_short_data),1),STG_short_data,'filled','jitter','on','MarkerFaceColor','k');
scatter(3*ones(length(STG_long_data),1),STG_long_data,'filled','jitter','on','MarkerFaceColor','k');

[~,contVSshort] = ttest2(control_data,STG_short_data);
[~,contVSlong] = ttest2(control_data,STG_long_data);
[~,shortVSlong] = ttest2(STG_short_data,STG_long_data);