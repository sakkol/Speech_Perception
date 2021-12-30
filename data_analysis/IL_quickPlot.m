function IL_quickPlot(Sbj_Metadata,curr_block)




load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_mtmconvol_pow.mat']),'epoched_wlt');
load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_info.mat']),'info')
events=info.events;
save_folder = fullfile(Sbj_Metadata.results,'Quick_plot_mtmconvol_pow');
if ~exist(save_folder,'dir'),mkdir(save_folder),end

% from fourierspectrum to powerspectrum
% cfg = [];
% cfg.output='abs';
% cfg.keeptrials = 'yes';
% epoched_wlt=ft_freqdescriptives(cfg,epoched_wlt);
% epoched_wlt.fourierspctrm=abs(epoched_wlt.fourierspctrm);
% epoched_wlt = renamefields(epoched_wlt,'fourierspctrm','powspctrm');
% epoched_wlt.dimord = 'rpt_chan_freq_time';

% for HFA
cfg=[];
cfg.frequency=[70 200];
cfg.avgoverfreq='no';
epoched_hfa = ft_selectdata(cfg, epoched_wlt);
cfg              = [];
cfg.baseline     = [-.5 -.05];
cfg.baselinetype = 'relative';
cfg.parameter    = 'powspctrm';
[epoched_hfa]    = ft_freqbaseline(cfg, epoched_hfa);
cfg=[];
cfg.frequency=[70 200];
cfg.avgoverfreq='yes';
epoched_hfa = ft_selectdata(cfg, epoched_hfa);
epoched_hfa = master_smoothData(epoched_hfa);

% for 
load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_mtmconvol_pow.mat']),'epoched_data');
epoched_data = master_smoothData(epoched_data);

isofrequency = events.cond_info{1}.Frequency;

% get the trial no for each condition
iso_4_sentence = IL_event_nos(events, 'all', 'Sentence', 'iso', '4-word');
iso_4_scrambled = IL_event_nos(events, 'all', 'Scrambled', 'iso', '4-word');
a_4_sentence = IL_event_nos(events, 'all', 'Sentence', 'a', '4-word');
a_4_scrambled = IL_event_nos(events, 'all', 'Scrambled', 'a', '4-word');
iso_35_sentence = IL_event_nos(events, 'all', 'Sentence', 'iso', '3-/5-word');
iso_35_scrambled = IL_event_nos(events, 'all', 'Scrambled', 'iso', '3-/5-word');
a_35_sentence = IL_event_nos(events, 'all', 'Sentence', 'a', '3-/5-word');
a_35_scrambled = IL_event_nos(events, 'all', 'Scrambled', 'a', '3-/5-word');
% 3 or 5 word sentences
nnnnnw=[];
for i=1:height(events)
    nnnnn=fieldnames(events.cfgs{i}.part2);
    nnnnnw(i,1)=str2double(nnnnn{end}(end));
end
nnnnnw3=find(nnnnnw==3);nnnnnw5=find(nnnnnw==5);
iso_3_sentence=intersect(iso_35_sentence,nnnnnw3);
iso_5_sentence=intersect(iso_35_sentence,nnnnnw5);
a_3_sentence=intersect(a_35_sentence,nnnnnw3);
a_5_sentence=intersect(a_35_sentence,nnnnnw5);
iso_3_scrambled=intersect(iso_35_scrambled,nnnnnw3);
iso_5_scrambled=intersect(iso_35_scrambled,nnnnnw5);
a_3_scrambled=intersect(a_35_scrambled,nnnnnw3);
a_5_scrambled=intersect(a_35_scrambled,nnnnnw5);


% 3-4-5 -- iso-a -- sentence-scrambled -- wlt-hfa-erp
loop1={'iso','a'};
loop2={'4','3','5'};
loop3={'sentence','scrambled'};
loop4={'wlt','hfa','data'};
sampstr=[];
datastr = [];datastr.epoched_wlt=epoched_wlt;datastr.epoched_hfa=epoched_hfa;datastr.epoched_data=epoched_data;
for l1=1:2
    for l2=1:3
        for l3=1:2
            for l4=1:3
                cfg                = [];
                cfg.latency        = [-.1 (str2double(loop2{l2})*5/isofrequency)+.1];
                eval(['cfg.trials         = ' loop1{l1} '_' loop2{l2} '_' loop3{l3} ';'])
                sampstr.([loop1{l1} '_' loop2{l2} '_' loop3{l3} '_' loop4{l4}]) = ft_selectdata(cfg,datastr.(['epoched_' loop4{l4}]));
            end
        end
    end
end

%% Now plot
for el = 1:length(epoched_wlt.label)
    figure('Units','normalized','Position', [0 0  1 1]);
    curr_label = epoched_wlt.label{el};

    for l4 = 1:3 % for rows
        for l2 = 1:3 % for columns
            subplot(3,3,subplotno(3,l4,l2))
            %% subplot for PSD
            if l4==1
                for_avg=squeeze(nanmean(sampstr.(['iso_' loop2{l2} '_sentence_' loop4{l4}]).powspctrm(:,el,:,:),4));
                hold on;
                if isnan(for_avg)
                    for_avg=nan([1,length(epoched_wlt.freq)]);
                elseif isvector(for_avg)
                    for_avg=for_avg';
                end
%                 s1=shadedErrorBar(epoched_wlt.freq,mean(for_avg,1),stderr(for_avg),'lineprops','r');
                s1=plot(epoched_wlt.freq,nanmean(for_avg,1),'r');

                for_avg=squeeze(nanmean(sampstr.(['iso_' loop2{l2} '_scrambled_' loop4{l4}]).powspctrm(:,el,:,:),4));
                hold on;
                if isnan(for_avg)
                    for_avg=nan([1,length(epoched_wlt.freq)]);
                elseif isvector(for_avg)
                    for_avg=for_avg';
                end
%                 s2=shadedErrorBar(epoched_wlt.freq,mean(for_avg,1),stderr(for_avg),'lineprops','b');
                s2=plot(epoched_wlt.freq,nanmean(for_avg,1),'b');

                for_avg=squeeze(nanmean(sampstr.(['a_' loop2{l2} '_sentence_' loop4{l4}]).powspctrm(:,el,:,:),4));
                hold on;
                if isnan(for_avg)
                    for_avg=nan([1,length(epoched_wlt.freq)]);
                elseif isvector(for_avg)
                    for_avg=for_avg';
                end
%                 s3=shadedErrorBar(epoched_wlt.freq,mean(for_avg,1),stderr(for_avg),'lineprops','g');
                s3=plot(epoched_wlt.freq,nanmean(for_avg,1),'g');

                for_avg=squeeze(nanmean(sampstr.(['a_' loop2{l2} '_scrambled_' loop4{l4}]).powspctrm(:,el,:,:),4));
                hold on;
                if isnan(for_avg)
                    for_avg=nan([1,length(epoched_wlt.freq)]);
                elseif isvector(for_avg)
                    for_avg=for_avg';
                end
%                 s4=shadedErrorBar(epoched_wlt.freq,mean(for_avg,1),stderr(for_avg),'lineprops','m');
                s4=plot(epoched_wlt.freq,nanmean(for_avg,1),'m');

                % additional stuff: frequency lines
                plot([isofrequency isofrequency], ylim,'k','HandleVisibility','off')
                if l2==1
                    plot([isofrequency/2 isofrequency/2], ylim,'k','HandleVisibility','off')
                    plot([isofrequency/4 isofrequency/4], ylim,'k','HandleVisibility','off')
                elseif l2==2
                    plot([isofrequency/3 isofrequency/3], ylim,'k','HandleVisibility','off')
                elseif l2==3
                    plot([isofrequency/5 isofrequency/5], ylim,'k','HandleVisibility','off')
                end

                %% Subplot for HFA
            elseif l4==2 % hfa
                for_avg=squeeze(nanmean(sampstr.(['iso_' loop2{l2} '_sentence_' loop4{l4}]).powspctrm(:,el,:,:),1));
                for_time = sampstr.(['iso_' loop2{l2} '_sentence_' loop4{l4}]);
                hold on;
                if isnan(for_avg)
                    for_avg=nan([1,length(for_time.time)]);
                elseif isvector(for_avg)
                    for_avg=for_avg';
                end
%                 s1=shadedErrorBar(epoched_wlt.freq,mean(for_avg,1),stderr(for_avg),'lineprops','r');
                s1=plot(for_time.time,nanmean(for_avg,1),'r');

                for_avg=squeeze(nanmean(sampstr.(['iso_' loop2{l2} '_scrambled_' loop4{l4}]).powspctrm(:,el,:,:),1));
                hold on;
                if isnan(for_avg)
                    for_avg=nan([1,length(for_time.time)]);
                elseif isvector(for_avg)
                    for_avg=for_avg';
                end
%                 s2=shadedErrorBar(epoched_wlt.freq,mean(for_avg,1),stderr(for_avg),'lineprops','b');
                s2=plot(for_time.time,nanmean(for_avg,1),'b');

                for_avg=squeeze(nanmean(sampstr.(['a_' loop2{l2} '_sentence_' loop4{l4}]).powspctrm(:,el,:,:),1));
                hold on;
                if isnan(for_avg)
                    for_avg=nan([1,length(for_time.time)]);
                elseif isvector(for_avg)
                    for_avg=for_avg';
                end
%                 s3=shadedErrorBar(epoched_wlt.freq,mean(for_avg,1),stderr(for_avg),'lineprops','g');
                s3=plot(for_time.time,nanmean(for_avg,1),'g');

                for_avg=squeeze(nanmean(sampstr.(['a_' loop2{l2} '_scrambled_' loop4{l4}]).powspctrm(:,el,:,:),1));
                hold on;
                if isnan(for_avg)
                    for_avg=nan([1,length(for_time.time)]);
                elseif isvector(for_avg)
                    for_avg=for_avg';
                end
%                 s4=shadedErrorBar(epoched_wlt.freq,mean(for_avg,1),stderr(for_avg),'lineprops','m');
                s4=plot(for_time.time,nanmean(for_avg,1),'m');

                % lines for sentence onsets
                plot([0 0], ylim,'k','HandleVisibility','off')
                for i=1:5 %% (((str2double(loop2{l2})*5+1)/isofrequency)/isofrequency)+.1])
                    plot([str2double(loop2{l2})*(i-1)/isofrequency str2double(loop2{l2})*(i-1)/isofrequency], [-.5 5],'k','HandleVisibility','off')
                end
            %% Subplot for ERP:
            else
                for_avg=squeeze(nanmean(sampstr.(['iso_' loop2{l2} '_sentence_' loop4{l4}]).trial(:,el,:),1));
                for_time = sampstr.(['iso_' loop2{l2} '_sentence_' loop4{l4}]);
                hold on;
                if isnan(for_avg)
                    for_avg=nan([1,length(for_time.time)]);
                elseif isvector(for_avg)
                    for_avg=for_avg';
                end
%                 s1=shadedErrorBar(sample_data.time,mean(for_avg,1),stderr(for_avg),'lineprops','r');
                s1=plot(for_time.time,nanmean(for_avg,1),'r');

                for_avg=squeeze(nanmean(sampstr.(['iso_' loop2{l2} '_scrambled_' loop4{l4}]).trial(:,el,:),1));
                hold on;
                if isnan(for_avg)
                    for_avg=nan([1,length(for_time.time)]);
                elseif isvector(for_avg)
                    for_avg=for_avg';
                end
%                 s2=shadedErrorBar(sample_data.time,mean(for_avg,1),stderr(for_avg),'lineprops','b');
                s2=plot(for_time.time,nanmean(for_avg,1),'b');

                for_avg=squeeze(nanmean(sampstr.(['a_' loop2{l2} '_sentence_' loop4{l4}]).trial(:,el,:),1));
                hold on;
                if isnan(for_avg)
                    for_avg=nan([1,length(for_time.time)]);
                elseif isvector(for_avg)
                    for_avg=for_avg';
                end
%                 s3=shadedErrorBar(sample_data.time,mean(for_avg,1),stderr(for_avg),'lineprops','g');
                s3=plot(for_time.time,nanmean(for_avg,1),'g');

                for_avg=squeeze(nanmean(sampstr.(['a_' loop2{l2} '_scrambled_' loop4{l4}]).trial(:,el,:),1));
                hold on;
                if isnan(for_avg)
                    for_avg=nan([1,length(for_time.time)]);
                elseif isvector(for_avg)
                    for_avg=for_avg';
                end
%                 s4=shadedErrorBar(sample_data.time,mean(for_avg,1),stderr(for_avg),'lineprops','m');
                s4=plot(for_time.time,nanmean(for_avg,1),'m');

                % lines for sentence onsets
                plot([0 0], ylim,'k','HandleVisibility','off')
                for i=1:5 %% (((str2double(loop2{l2})*5+1)/isofrequency)/isofrequency)+.1])
                    plot([str2double(loop2{l2})*(i-1)/isofrequency str2double(loop2{l2})*(i-1)/isofrequency], ylim,'k','HandleVisibility','off')
                end
            end

            %% Additional stuff: legend, title, xlim
%             legend([s1.mainLine,s2.mainLine,s3.mainLine,s4.mainLine],{'Iso-Sentence','Iso-Scrambled','Achr-Sentence','Achr-Scrambled'})
            legend([s1,s2,s3,s4],{'Iso-Sentence','Iso-Scrambled','Achr-Sentence','Achr-Scrambled'})
            
            %info: loop2={'4','3','5'};loop4={'wlt','hfa','data'};
            if l4==1
                xlim([epoched_wlt.freq(1) 6])
                title('PSD')
                if l2==1
                    title(['4-word sentences; word frequency ' num2str(isofrequency) 'Hz'])
                elseif l2==2
                    title(['3-word sentences; word frequency ' num2str(isofrequency) 'Hz'])
                elseif l2==3
                    title(['5-word sentences; word frequency ' num2str(isofrequency) 'Hz'])
                end
                ylabel('Power')
                xlabel('Frequency (Hz)')
            elseif l4==2
                title('HFA')
                xlim([-.1 ((str2double(loop2{l2}))*5/isofrequency)+.1])
                ylim([.5 5])
                ylabel('Power (relative)')
                xlabel('Time (s)')
            elseif l4==3
                title('ERP')
                xlim([-.1 ((str2double(loop2{l2}))*5/isofrequency)+.1])
                ylabel('Voltage (\muV)')
                xlabel('Time (s)')
            end
        end
    end
    
    xxx=makeBlockInfo(Sbj_Metadata,curr_block);
    ax = axes;
    ax.Visible = 'off';
    text(gca,.5,1.07,[curr_label ' in ' info.channelinfo.Desikan_Killiany{el} '; ' xxx.Stim_type{1} ' Isochronous Listening Task; block: ' curr_block],...
        'Units','normalized','FontSize',18,'FontWeight','bold','HorizontalAlignment','center')
%     sgtitle([curr_label ' in ' info.channelinfo.Desikan_Killiany{el} '; Isochronous Listening Task; block: ' curr_block])
    
    % Save the figure
    fprintf('\t-Saving electrode #%d %s, out of %d\n',el,curr_label,length(epoched_wlt.label))
    print(fullfile(save_folder,[curr_label , '_' curr_block '_Qplot_mtmconvol_pow.jpg']),'-djpeg','-r300')
    close all
    
end
