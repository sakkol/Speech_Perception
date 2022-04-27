function CloseAndSave(t)
data_to_save = get(t,'Data');
save('tmp.mat','data_to_save')
% for reference and some more info see: https://www.mathworks.com/matlabcentral/answers/509960-i-can-t-retrieve-the-user-input-information-in-a-table-when-the-user-closes-the-table
end