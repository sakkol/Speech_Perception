data_root = '/media/sakkol/HDD1/HBML/';
project_name = 'IsochronousListening';

AllBlockInfo = readtable(fullfile(data_root,'PROJECTS_DATA',project_name,[project_name '_BlockInfo.xlsx']));
subjects = unique(AllBlockInfo.sbj_ID);
subjects = subjects(~cellfun(@isempty,subjects));

plot_col=[];
all_cond_accr=[];
for s = 1:length(subjects)
    
    sbj_ID = subjects{s};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    
    % Get params directly from BlockList excel sheet
    curr_block = AllBlockInfo.BlockList(strcmp(AllBlockInfo.sbj_ID,sbj_ID) & AllBlockInfo.preproc_FU==1);
    if isempty(curr_block)
        subjects{s}='';
        continue
    else
        curr_block = curr_block{1};
        plot_col(s,:) = [1.000000000000000   0.411764705882353   0.160784313725490];
    end
    fprintf('Running for %s\n',sbj_ID)
%     StimSide{end+1,1} = upper(AllBlockInfo.StimSide{strcmp(AllBlockInfo.sbj_ID,sbj_ID) & strcmp(AllBlockInfo.BlockList,curr_block)});
%     StimMethod{end+1,1} = lower(AllBlockInfo.StimLoc{strcmp(AllBlockInfo.sbj_ID,sbj_ID) & strcmp(AllBlockInfo.BlockList,curr_block)});
%     Lang{end+1,1} = lower(AllBlockInfo.Language{strcmp(AllBlockInfo.sbj_ID,sbj_ID) & strcmp(AllBlockInfo.BlockList,curr_block)});
    
    % Load info table
    load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']));
    
    % Read conditions:
    events = info.events;
    [event_types] = IL_event_types(events);
    
    cond_name=[];cond_accr=[];
loop1={'iso','a'};
loop2={4,3,5};
loop3={'sentence','scrambled'};
for l1=1:2
    for l2=1:3
        for l3=1:2
            cond_name{end+1,1} = strjoin({loop3{l3}, loop1{l1}, num2str(loop2{l2})},'-');
            cond_accr(end+1,1) = length(IL_get_eventsOI(events, loop3{l3}, loop1{l1}, loop2{l2}, 0)); % only the accurate trials
            cond_accr(end,2)   = length(IL_get_eventsOI(events, loop3{l3}, loop1{l1}, loop2{l2}, 1)); % only the accurate trials
        end
    end
end

all_cond_accr{end+1,1} = cond_accr;

end

subjects = subjects(~cellfun(@isempty,subjects));
plot_col(plot_col(:,1)==0,:) = [];

%% plot the accuracy of iso vs a combined - V1

isoVSa{1} = contains(cond_name,'-iso-');
isoVSa{2} = contains(cond_name,'-a-');
isoVSa_accr = zeros(length(all_cond_accr),2);

figure('Position',[0 0 800 1000])
for s = 1:length(all_cond_accr)
    for b = 1:2
        isoVSa_accr(s,b) = sum(all_cond_accr{s}(isoVSa{b},2));
    end
    plot(1:2,isoVSa_accr(s,:),'-s','Color',plot_col(s,:),'LineWidth',4,'MarkerFaceColor',plot_col(s,:),'MarkerSize',14)
    hold on
end

xlim([.5 2.5])
ylim([28 40])
yticks([30:5:40])
set(gca,'FontSize',24,'FontWeight','bold')
set(gca,'XTick',1:2,'xticklabels',{'',''})
% set(gca,'xticklabels',{'No-stimulation','50ms\\Delay','200ms\\newlineDelay'})
grid on
box off


savedir = '/media/sakkol/HDD1/HBML/PROJECTS_DATA/IsochronousListening/Collected_Results/beh_results';
print(fullfile(savedir,'beh_results_v1.png'),'-dpng','-r300')



%% plot the accuracy of iso vs a combined - V2

isoVSa{1} = contains(cond_name,'-iso-');
isoVSa{2} = contains(cond_name,'-a-');
isoVSa_accr = zeros(length(all_cond_accr),2);

figure('Position',[0 0 800 1000])
for s = 1:length(all_cond_accr)
    for b = 1:2
        isoVSa_accr(s,b) = sum(all_cond_accr{s}(isoVSa{b},2));
    end
    plot(1:2,isoVSa_accr(s,:)/.4,'-s','Color',plot_col(s,:),'LineWidth',4,'MarkerFaceColor',plot_col(s,:),'MarkerSize',14)
    hold on
end

xlim([.8 2.2])
ylim([70 100])
yticks([70:10:100])
set(gca,'FontSize',26,'FontWeight','bold')
set(gca,'XTick',1:2,'xticklabels',{'Isochr.','Achr.'})
ylabel('Accuracy (%)')
% set(gca,'xticklabels',{'No-stimulation','50ms\\Delay','200ms\\newlineDelay'})
grid on
box off


savedir = '/media/sakkol/HDD1/HBML/PROJECTS_DATA/IsochronousListening/Collected_Results/beh_results';
print(fullfile(savedir,'beh_results_v2.png'),'-dpng','-r300')
