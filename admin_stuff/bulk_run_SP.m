%% Select 'subjects'
data_root = '/media/sakkol/HDD1/HBML/';
project_name = 'Speech_Perception';

AllBlockInfo = readtable(fullfile(data_root,'PROJECTS_DATA',project_name,[project_name '_BlockInfo.xlsx']));
subjects = unique(AllBlockInfo.sbj_ID);

[indx,~] = listdlg('ListString',subjects);


%% Behavioral
% create subjects in the first section
for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    
    for c = 1:size(Sbj_Metadata.BlockLists,1)
        
        curr_block = Sbj_Metadata.BlockLists{c};
        fprintf('Running for Subject:%s, Block:%s of %d\n',sbj_ID,curr_block,size(Sbj_Metadata.BlockLists,1))
        
        SP_beh_analysis(Sbj_Metadata,curr_block)
        close all
    end
end

%% Bulk run word evoked ERP-HFA plots
for s = 1:length(indx)
    sbj_ID = subjects{indx(s)};
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    SP_word_evoked(Sbj_Metadata)
end
