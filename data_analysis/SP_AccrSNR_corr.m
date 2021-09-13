function [sbj_block] = SP_AccrSNR_corr(sbj_block)
%% to correlate the Accuracy and SNR
data_root = '/media/sakkol/HDD1/HBML/';
project_name = 'Speech_Perception';
if ~exist('sbj_block','var') || isempty(sbj_block)
    
    AllBlockInfo = readtable(fullfile(data_root,'PROJECTS_DATA',project_name,[project_name '_BlockInfo.xlsx']));
    sbj_block = AllBlockInfo(contains(AllBlockInfo.conditions_code,'1')&AllBlockInfo.preproc_FU==1,1:2);
    sbj_block = table2cell(sbj_block);
end

%% gather the values
Accr = zeros(height(sbj_block),1);
SNR = zeros(height(sbj_block),1);
for s = 1:height(sbj_block)
    sbj_ID = sbj_block{s,1};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    
    curr_block = sbj_block{s,2};
    % load one of the info files
    load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']),'info');
    
    Accr(s) = sum(info.events.Acc_word_count(info.events.Cond_code==1));
    try
        SNR(s) = info.events.trial_details{1}.SNR;
    catch
        info.events.trial_details = struct2table(info.events.trial_details);
        SNR(s) = info.events.trial_details.SNR(1);
    end
end

%% plot
figure('Position',[0,0,700,700])
scatter(SNR,Accr)
l=lsline;l.LineWidth = 2;l.Color = 'k';
xlabel('SNR')
ylabel('Accuracy (%)')
set(gca, 'FontSize',13,'FontWeight','bold');
box on
[rho,pval] = corr(SNR,Accr);
text(0.04,0.96,['r=' sprintf('%.3f',rho) ';p=' sprintf('%.3f',pval)], 'Units','normalized','FontSize',13,'FontWeight','bold')

end