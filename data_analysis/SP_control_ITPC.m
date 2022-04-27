function SP_control_ITPC(Sbj_Metadata)

% Select blocks to import
control_blocks = select_cont_blocks(Sbj_Metadata);
save_dir = fullfile(Sbj_Metadata.results, [strjoin(control_blocks,'_') '_v2']);
fprintf('Loading ''fouri_of_words'' from:\n-->%s\n',fullfile(save_dir, [strjoin(control_blocks,'_') '_ctrl_word_fouri.mat']))
load(fullfile(save_dir, [strjoin(control_blocks,'_') '_ctrl_word_fouri.mat']));
load(fullfile(Sbj_Metadata.iEEG_data,Sbj_Metadata.BlockLists{1},[Sbj_Metadata.BlockLists{1} '_info.mat']),'info');

freq = fouri_of_words.freq_band_dtls{1};
time = linspace(-0.05,0.7,76);

save_dir = fullfile(Sbj_Metadata.results, [strjoin(control_blocks,'_') '_v2'],'PICS');
if ~exist(save_dir,'dir'),mkdir(save_dir),end

bwr = load('bwr_cmap.mat');

itpc = [];
for cond = 1:3
    % get data
    if cond == 1
        curr_fouri_all = fouri_of_words.corr_rspn_fouri_word{1};
    elseif cond == 2
        curr_fouri_all = fouri_of_words.no_rspn_fouri_word{1};
    else
        curr_fouri_all = fouri_of_words.wrng_rspn_fouri_word{1};
    end
    
    % compute inter-trial phase coherence (itpc) for each conditions
    tmp      = curr_fouri_all./abs(curr_fouri_all);    % divide by amplitude
    tmp      = sum(tmp,1);                            % sum angles across trials
    tmp      = abs(tmp)/size(curr_fouri_all,1);       % take the absolute value and normalize
    itpc(cond,:,:,:) = squeeze(tmp);                          % remove the first singleton dimension
    
    
end

% plot ITPC of each electrode
for el = 1:size(itpc,2)
    figure('Units','normalized','Position', [0 0  .6 .3]);
    for cond = 1:3
        
        if cond == 1
            totitle = 'Correct responses';
        elseif cond == 2
            totitle = 'No responses';
        else
            totitle = 'Wrong responses';
        end
        
        subplot(1, 3, cond);
        imagesc(time, freq, squeeze(itpc(cond,el,:,:)));
        axis xy
        title(totitle);
        set(gca, 'FontSize',13,'FontWeight','bold');
        if cond==1,ylabel('Frequency (Hz)');end
        if cond==2,xlabel('Time (s)');end
        caxis([0 0.5])
    end
    colormap(bwr.rgb_vals)
    colorbar
    sgtitle(['Elec: ' info.channelinfo.Label{el} ' - inter-trial phase coherence (word onset locked)'], 'FontSize',15,'FontWeight','bold')
    print('-r300','-djpeg',fullfile(save_dir,[info.channelinfo.Label{el} '_ITPC_word.jpg']))
    close all
end


end