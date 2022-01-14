function IL_power_peaks(Sbj_Metadata,curr_block,elec,ERP_HFA_TF)

if ~exist('ERP_HFA_TF','var') || isempty(ERP_HFA_TF)
    ERP_HFA_TF = 'elec_TF';
end

load(fullfile(Sbj_Metadata.iEEG_data, curr_block, ERP_HFA_TF, [elec '_' curr_block '_mtmconvol_fourier.mat']),'epoched_wlt');
load(fullfile(Sbj_Metadata.iEEG_data, curr_block, [curr_block '_info.mat']),'info')
events=info.events;
save_folder = fullfile(Sbj_Metadata.results,'power_peaks');
if ~exist(save_folder,'dir'),mkdir(save_folder),end

% from fourierspectrum to powerspectrum
cfg = [];
cfg.output='abs';
cfg.keeptrials = 'yes';
epoched_wlt=ft_freqdescriptives(cfg,epoched_wlt);
epoched_wlt.dimord = 'rpt_chan_freq_time';

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

% load ERP
load(fullfile(Sbj_Metadata.iEEG_data, curr_block, ERP_HFA_TF, [elec '_' curr_block '_mtmconvol_fourier.mat']),'epoched_data');
epoched_data = master_smoothData(epoched_data);

isofrequency = events.cond_info{1}.Frequency;

% 3-4-5 -- iso-a -- sentence-scrambled -- wlt-hfa-erp
loop1={'iso','a'};
loop2={4,3,5};
loop3={'sentence','scrambled'};
loop4={'wlt','hfa','data'};
sampstr=[];
datastr = [];datastr.epoched_wlt=epoched_wlt;datastr.epoched_hfa=epoched_hfa;datastr.epoched_data=epoched_data;
for l1=1:2
    for l2=1:3
        for l3=1:2
            for l4=1:3
                cfg                = [];
                cfg.latency        = [-.1 (loop2{l2}*5/isofrequency)+.1];
                cfg.trials         = IL_get_eventsOI(events, loop3{l3}, loop1{l1}, loop2{l2}, [1 0]); % only the accurate trials
                sampstr.([loop1{l1} '_' num2str(loop2{l2}) '_' loop3{l3} '_' loop4{l4}]) = ft_selectdata(cfg,datastr.(['epoched_' loop4{l4}]));
            end
        end
    end
end

ranksumresults=[];
%% Now plot
for el = 1:length(epoched_wlt.label)
    figure('Units','normalized','Position', [0 0  1 1]);
    curr_label = epoched_wlt.label{el};

    for l4 = 1:3 % for rows
        for l2 = 1:3 % for columns
            subplot(3,3,subplotno(3,l4,l2))
            %% subplot for PSD
            if l4==1
                % first iso sentence
                [s1, for_avg] = IL_powerpeak_plot(sampstr,['iso_' num2str(loop2{l2}) '_sentence_' loop4{l4}],'r');
                [~,ranksumresults.(['iso_' num2str(loop2{l2}) '_sentence_24Hz'])] = IL_powerpeak_stat(for_avg,epoched_wlt.freq,s1,isofrequency);
                if l2==1
                    [~,ranksumresults.(['iso_' num2str(loop2{l2}) '_sentence_12Hz'])] = IL_powerpeak_stat(for_avg,epoched_wlt.freq,s1,isofrequency/2);
                    [~,ranksumresults.(['iso_' num2str(loop2{l2}) '_sentence_06Hz'])] = IL_powerpeak_stat(for_avg,epoched_wlt.freq,s1,isofrequency/4);
                elseif l2==2
                    [~,ranksumresults.(['iso_' num2str(loop2{l2}) '_sentence_08Hz'])] = IL_powerpeak_stat(for_avg,epoched_wlt.freq,s1,isofrequency/3);
                elseif l2==3
                    [~,ranksumresults.(['iso_' num2str(loop2{l2}) '_sentence_048Hz'])] = IL_powerpeak_stat(for_avg,epoched_wlt.freq,s1,isofrequency/5);
                end
                
                % second iso scr
                [s2, for_avg] = IL_powerpeak_plot(sampstr,['iso_' num2str(loop2{l2}) '_scrambled_' loop4{l4}],'b');
                [~,ranksumresults.(['iso_' num2str(loop2{l2}) '_scrambled_24Hz'])] = IL_powerpeak_stat(for_avg,epoched_wlt.freq,s2,isofrequency);
                if l2==1
                    [~,ranksumresults.(['iso_' num2str(loop2{l2}) '_scrambled_12Hz'])] = IL_powerpeak_stat(for_avg,epoched_wlt.freq,s2,isofrequency/2);
                    [~,ranksumresults.(['iso_' num2str(loop2{l2}) '_scrambled_06Hz'])] = IL_powerpeak_stat(for_avg,epoched_wlt.freq,s2,isofrequency/4);
                elseif l2==2
                    [~,ranksumresults.(['iso_' num2str(loop2{l2}) '_scrambled_08Hz'])] = IL_powerpeak_stat(for_avg,epoched_wlt.freq,s2,isofrequency/3);
                elseif l2==3
                    [~,ranksumresults.(['iso_' num2str(loop2{l2}) '_scrambled_048Hz'])] = IL_powerpeak_stat(for_avg,epoched_wlt.freq,s2,isofrequency/5);
                end
                
                % third a sent
                [s3, for_avg] = IL_powerpeak_plot(sampstr,['a_' num2str(loop2{l2}) '_sentence_' loop4{l4}],'g');
                [~,ranksumresults.(['a_' num2str(loop2{l2}) '_sentence_24Hz'])] = IL_powerpeak_stat(for_avg,epoched_wlt.freq,s3,isofrequency);
                if l2==1
                    [~,ranksumresults.(['a_' num2str(loop2{l2}) '_sentence_12Hz'])] = IL_powerpeak_stat(for_avg,epoched_wlt.freq,s3,isofrequency/2);
                    [~,ranksumresults.(['a_' num2str(loop2{l2}) '_sentence_06Hz'])] = IL_powerpeak_stat(for_avg,epoched_wlt.freq,s3,isofrequency/4);
                elseif l2==2
                    [~,ranksumresults.(['a_' num2str(loop2{l2}) '_sentence_08Hz'])] = IL_powerpeak_stat(for_avg,epoched_wlt.freq,s3,isofrequency/3);
                elseif l2==3
                    [~,ranksumresults.(['a_' num2str(loop2{l2}) '_sentence_048Hz'])] = IL_powerpeak_stat(for_avg,epoched_wlt.freq,s3,isofrequency/5);
                end
                
                % fourth a scr
                [s4, for_avg] = IL_powerpeak_plot(sampstr,['a_' num2str(loop2{l2}) '_scrambled_' loop4{l4}],'m');
                [~,ranksumresults.(['a_' num2str(loop2{l2}) '_scrambled_24Hz'])] = IL_powerpeak_stat(for_avg,epoched_wlt.freq,s4,isofrequency);
                if l2==1
                    [~,ranksumresults.(['a_' num2str(loop2{l2}) '_scrambled_12Hz'])] = IL_powerpeak_stat(for_avg,epoched_wlt.freq,s4,isofrequency/2);
                    [~,ranksumresults.(['a_' num2str(loop2{l2}) '_scrambled_06Hz'])] = IL_powerpeak_stat(for_avg,epoched_wlt.freq,s4,isofrequency/4);
                elseif l2==2
                    [~,ranksumresults.(['a_' num2str(loop2{l2}) '_scrambled_08Hz'])] = IL_powerpeak_stat(for_avg,epoched_wlt.freq,s4,isofrequency/3);
                elseif l2==3
                    [~,ranksumresults.(['a_' num2str(loop2{l2}) '_scrambled_048Hz'])] = IL_powerpeak_stat(for_avg,epoched_wlt.freq,s4,isofrequency/5);
                end

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
                IL_powerpeak_plot(sampstr,['iso_' num2str(loop2{l2}) '_sentence_' loop4{l4}],'r');
                IL_powerpeak_plot(sampstr,['iso_' num2str(loop2{l2}) '_scrambled_' loop4{l4}],'b');
                IL_powerpeak_plot(sampstr,['a_' num2str(loop2{l2}) '_sentence_' loop4{l4}],'g');
                IL_powerpeak_plot(sampstr,['a_' num2str(loop2{l2}) '_scrambled_' loop4{l4}],'m');

                % lines for sentence onsets
                plot([0 0], ylim,'k','HandleVisibility','off')
                for i=1:5 %% (((loop2{l2}*5+1)/isofrequency)/isofrequency)+.1])
                    plot([loop2{l2}*(i-1)/isofrequency loop2{l2}*(i-1)/isofrequency], [-.5 5],'k','HandleVisibility','off')
                end
            %% Subplot for ERP:
            else
                IL_powerpeak_plot(sampstr,['iso_' num2str(loop2{l2}) '_sentence_' loop4{l4}],'r');
                IL_powerpeak_plot(sampstr,['iso_' num2str(loop2{l2}) '_scrambled_' loop4{l4}],'b');
                IL_powerpeak_plot(sampstr,['a_' num2str(loop2{l2}) '_sentence_' loop4{l4}],'g');
                IL_powerpeak_plot(sampstr,['a_' num2str(loop2{l2}) '_scrambled_' loop4{l4}],'m');

                % lines for sentence onsets
                plot([0 0], ylim,'k','HandleVisibility','off')
                for i=1:5 %% (((loop2{l2}*5+1)/isofrequency)/isofrequency)+.1])
                    plot([loop2{l2}*(i-1)/isofrequency loop2{l2}*(i-1)/isofrequency], ylim,'k','HandleVisibility','off')
                end
            end

            %% Additional stuff: legend, title, xlim
%             legend([s1.mainLine,s2.mainLine,s3.mainLine,s4.mainLine],{'Iso-Sentence','Iso-Scrambled','Achr-Sentence','Achr-Scrambled'})
            
            %info: loop2={'4','3','5'};loop4={'wlt','hfa','data'};
            if l4==1
                xlim([epoched_wlt.freq(1) 5])
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
                if l2==3
                    legend([s1,s2,s3,s4],{'Iso-Sentence','Iso-Scrambled','Achr-Sentence','Achr-Scrambled'})
                end
            elseif l4==2
                title('HFA')
                xlim([-.1 ((loop2{l2})*5/isofrequency)+.1])
                ylim([.5 5])
                ylabel('Power (relative)')
                xlabel('Time (s)')
            elseif l4==3
                title('ERP')
                xlim([-.1 ((loop2{l2})*5/isofrequency)+.1])
                ylabel('Voltage (\muV)')
                xlabel('Time (s)')
            end
        end
    end
    
    xxx=makeBlockInfo(Sbj_Metadata,curr_block);
    ax = axes;
    ax.Visible = 'off';
    text(gca,.5,1.07,[curr_label ' in ' info.channelinfo.Desikan_Killiany{strcmp(info.channelinfo.Label,elec)} '; ' xxx.Stim_type{1} ' task; block: ' curr_block],...
        'Units','normalized','FontSize',18,'FontWeight','bold','HorizontalAlignment','center','Interpreter','none')
%     sgtitle([curr_label ' in ' info.channelinfo.Desikan_Killiany{el} '; Isochronous Listening Task; block: ' curr_block])
    
    % Save the figure
    fprintf('\t-Saving electrode #%d %s, out of %d\n',el,curr_label,length(epoched_wlt.label))
    print(fullfile(save_folder,[curr_label , '_' curr_block '_' ERP_HFA_TF '_powerstats.jpg']),'-djpeg','-r300')
    save(fullfile(save_folder,[curr_label , '_' curr_block '_' ERP_HFA_TF '_powerstats.mat']),'ranksumresults')
    close all
    
end
