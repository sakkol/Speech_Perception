%% Load table
sentence_table = readtable('C:\Users\user\Desktop\PROJECT MANAGEMENT\PhD\TASK\Matrix_Sentence_Table.xlsx');

%% create many many sentences and write them to individual files
save_dir = 'C:\Users\user\Desktop\PROJECT MANAGEMENT\PhD\TASK\Matrix_Sentences\All_Texts';
i=1;
collect_fname{100000,2}='';
fprintf('Iteration = 000000');
for a=1:10
    for b=1:10
        for c=1:10
            for d=1:10
                for e=1:10
                    fprintf('\b\b\b\b\b\b%6d',i)
                    
                    curr_sent = [sentence_table.Name{a} ' ' sentence_table.Verb{b} ' ' sentence_table.Number{c} ...
                        ' ' sentence_table.Adjective{d} ' ' sentence_table.Noun{e} '.'];
                    curr_sent_fname = strrep(curr_sent,' ','_');
                    curr_sent_fname = erase(curr_sent_fname,'.');
                    collect_fname{i,1} = [save_dir filesep curr_sent_fname '.txt'];
                    collect_fname{i,2} = [curr_sent_fname '.txt'];
                    fid = fopen([save_dir filesep curr_sent_fname '.txt'],'w+');
                    fwrite(fid,curr_sent);
                    fclose(fid);
                    i=i+1;
                end
            end
        end
    end
end
fprintf('\n')

%% Randomly select 400
selected_save_dir = 'C:\Users\user\Desktop\PROJECT MANAGEMENT\PhD\TASK\Matrix_Sentences\Selected_Texts';

randomnos = randi([1 100000],1,400);

for r=1:400
    copyfile(collect_fname{randomnos(r),1}, [selected_save_dir filesep collect_fname{randomnos(r),1}]);
end