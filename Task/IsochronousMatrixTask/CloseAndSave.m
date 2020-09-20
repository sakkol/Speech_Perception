function CloseAndSave(t)
selections_table = get(t,'Data');
save('tmp.mat','selections_table')
% for reference and some more info see: https://www.mathworks.com/matlabcentral/answers/509960-i-can-t-retrieve-the-user-input-information-in-a-table-when-the-user-closes-the-table
end