control_mean = mean(real1_response_table.Acc_word_count(strcmp(real1_response_table.Condition_Code,'1')));
control_std = std(real1_response_table.Acc_word_count(strcmp(real1_response_table.Condition_Code,'1')));

STG_short_mean = mean(real1_response_table.Acc_word_count(strcmp(real1_response_table.Condition_Code,'4')));
STG_short_std = std(real1_response_table.Acc_word_count(strcmp(real1_response_table.Condition_Code,'4')));

STG_long_mean = mean(real1_response_table.Acc_word_count(strcmp(real1_response_table.Condition_Code,'5')));
STG_long_std = std(real1_response_table.Acc_word_count(strcmp(real1_response_table.Condition_Code,'5')));

control_data = real1_response_table.Acc_word_count(strcmp(real1_response_table.Condition_Code,'1'));
STG_short_data = real1_response_table.Acc_word_count(strcmp(real1_response_table.Condition_Code,'4'));
STG_long_data = real1_response_table.Acc_word_count(strcmp(real1_response_table.Condition_Code,'5'));

control_corr = ismember(real1_response_table{strcmp(real1_response_table.Condition_Code,'1'),4:8},'1');
STG_short_corr = ismember(real1_response_table{strcmp(real1_response_table.Condition_Code,'4'),4:8},'1');
STG_long_corr = ismember(real1_response_table{strcmp(real1_response_table.Condition_Code,'5'),4:8},'1');

bar([sum(sum(control_corr)),sum(sum(STG_short_corr)),sum(sum(STG_long_corr))])

figure
xlim([0.5 3.5])
hold on
boxplot([control_data; STG_short_data;STG_long_data],...
    [ones(length(control_data),1);2*ones(length(STG_short_data),1);3*ones(length(STG_long_data),1)],...
    'Widths',0.45,'Colors','b')

scatter(ones(length(control_data),1),control_data,'filled','jitter','on','MarkerFaceColor','k');
scatter(2*ones(length(STG_short_data),1),STG_short_data,'filled','jitter','on','MarkerFaceColor','k');
scatter(3*ones(length(STG_long_data),1),STG_long_data,'filled','jitter','on','MarkerFaceColor','k');

[~,contVSshort] = ttest2(control_data,STG_short_data)
[~,contVSlong] = ttest2(control_data,STG_long_data)
[~,shortVSlong] = ttest2(STG_short_data,STG_long_data)

cont_table = [sum(sum(control_corr)), sum(sum(STG_long_corr));...
    100-sum(sum(control_corr)),100-sum(sum(STG_long_corr))];
[h,p,stats] = fishertest(cont_table)

%% Real2
control_mean = mean(real2_response_table.Acc_word_count(real2_response_table.Condition_Code==1));
control_std = std(real2_response_table.Acc_word_count(real2_response_table.Condition_Code==1));

epic_short_mean = mean(real2_response_table.Acc_word_count(real2_response_table.Condition_Code==6));
epic_short_std = std(real2_response_table.Acc_word_count(real2_response_table.Condition_Code==6));

epic_long_mean = mean(real2_response_table.Acc_word_count(real2_response_table.Condition_Code==7));
epic_long_std = std(real2_response_table.Acc_word_count(real2_response_table.Condition_Code==7));


control_data = real2_response_table.Acc_word_count(real2_response_table.Condition_Code==1);
epic_short_data = real2_response_table.Acc_word_count(real2_response_table.Condition_Code==6);
epic_long_data = real2_response_table.Acc_word_count(real2_response_table.Condition_Code==7);

figure
xlim([0.5 3.5])
hold on
boxplot([control_data; epic_short_data;epic_long_data],...
    [ones(length(control_data),1);2*ones(length(epic_short_data),1);3*ones(length(epic_long_data),1)],...
    'Widths',0.45,'Colors','b')

scatter(ones(length(control_data),1),control_data,'filled','jitter','on','MarkerFaceColor','k');
scatter(2*ones(length(epic_short_data),1),epic_short_data,'filled','jitter','on','MarkerFaceColor','k');
scatter(3*ones(length(epic_long_data),1),epic_long_data,'filled','jitter','on','MarkerFaceColor','k');

[~,contVSshort] = ttest2(control_data,epic_short_data)
[~,contVSlong] = ttest2(control_data,epic_long_data)
[~,shortVSlong] = ttest2(epic_short_data,epic_long_data)

control_corr = ismember(real2_response_table{real2_response_table.Condition_Code==1,4:8},'1');
epic_short_corr = ismember(real2_response_table{real2_response_table.Condition_Code==6,4:8},'1');
epic_long_corr = ismember(real2_response_table{real2_response_table.Condition_Code==7,4:8},'1');

bar([sum(sum(control_corr)),sum(sum(epic_short_corr)),sum(sum(epic_long_corr))])

wilx_contVSshort = ranksum(control_corr(:),epic_short_corr(:))
wilx_contVSlong = ranksum(control_corr(:),epic_long_corr(:))
wilx_shortVSlong = ranksum(epic_short_corr(:),epic_long_corr(:))

cont_table = [sum(sum(control_corr)), sum(sum(epic_long_corr));...
    100-sum(sum(control_corr)),100-sum(sum(epic_long_corr))];
[h,p,stats] = fishertest(cont_table)

%% Bootstrapping

bootsam_cont = bootstrp(1000,@mean, control_corr(:));
bootsam_epic_short = bootstrp(1000,@mean, epic_short_corr(:));
bootsam_epic_long = bootstrp(1000,@mean, epic_long_corr(:));

[x,boot_contVSshort,~,tablestat] = ttest2(bootsam_cont,bootsam_epic_short)
[x,boot_contVSlong,~,tablestat] = ttest2(bootsam_cont,bootsam_epic_long)

nReps = 100000;
perm = zeros(nReps,1);
for i=1:nReps
    %shuffle the lsat scores and recalculate the correlation
    perm(i) = ttest2([control_corr(:);epic_long_corr(:)],shuffle([ones(length(control_corr(:)),1);2*ones(length(epic_long_corr(:)),1)]));
end


p = sum(perm>sampStat)/nReps;

[pval, t_orig, ~, ~, ~]=mult_comp_perm_t2(control_corr(:),epic_long_corr(:),10000,0);

