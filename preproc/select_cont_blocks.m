function [control_blocks] = select_cont_blocks(Sbj_Metadata)
% Select blocks to import
blocklistsheet = [char(Sbj_Metadata.project_root) filesep char(Sbj_Metadata.project_name) '_BlockInfo.xlsx']; % "F:\HBML\PROJECTS_DATA\CL_Train\CL_Train_BlockLists.xlsx";
blocklistall = readtable(blocklistsheet);
control_blocks = blocklistall.BlockList(strcmpi(blocklistall.sbj_ID,Sbj_Metadata.sbj_ID) & contains(blocklistall.conditions_code,'1'));
end