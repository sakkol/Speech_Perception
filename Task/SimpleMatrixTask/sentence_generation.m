%% Load table
sentence_table = readtable('/home/sakkol/Documents/Speech_Perception_stim/Matrix_Sentence_Table.xlsx');

%% create many many sentences and write them to individual files
save_dir = '/home/sakkol/Documents/Speech_Perception_stim/All_Texts';
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

%% Randomly select 404
selected_save_dir = 'C:\Users\user\Desktop\PROJECT MANAGEMENT\PhD\TASK\Matrix_Sentences\02nd_Generation\Selected_Texts_2ndGen';

randomnos = randi([1 100000],1,404);

for r=1:810
    copyfile(collect_fname{randomnos(r),1}, [selected_save_dir filesep collect_fname{randomnos(r),2}]);
end

%% Randomly select 4000 sentences directly from All_Texts
all_text_dir = dir('/home/sakkol/Documents/Speech_Perception_stim/All_Texts/*.txt');
selected_save_dir = '/home/sakkol/Documents/Speech_Perception_stim/4th_Generation/Selected_Texts_4thGen';

randomnos = randi([1 100000],1,2000);

for i=1:2000
    to_copy = fullfile(all_text_dir(randomnos(i)).folder,all_text_dir(randomnos(i)).name);
    to_save = fullfile(selected_save_dir,all_text_dir(randomnos(i)).name);
    copyfile(to_copy,to_save)
end

%% Spanish version:
% Load table
sentence_table = readtable('/home/sakkol/Documents/Spanish_Matrix_Sentence/spanish_matrix_words.xlsx');

% create many many sentences and write them to individual files
save_dir = '/home/sakkol/Documents/Spanish_Matrix_Sentence/Version_1/All_Texts';
i=1;
collect_fname{100000,2}='';
fprintf('Iteration = 000000');
for a=1:10
    for b=1:10
        for c=1:10
            for d=1:10
                for e=1:10
                    fprintf('\b\b\b\b\b\b%6d',i)
                    
                    curr_sent = [sentence_table.Name{a} ' ' sentence_table.Verb{b} ' ' sentence_table.Numeral{c} ...
                        ' ' sentence_table.Object{d} ' ' sentence_table.Adjective{e} '.'];
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

%% Randomly select 10000 sentences directly from All_Texts
all_text_dir = dir('/home/sakkol/Documents/Spanish_Matrix_Sentence/Version_2/All_Texts/*.txt');
selected_save_dir = '/home/sakkol/Documents/Spanish_Matrix_Sentence/Version_2/Selected_20000';

randomnos = randperm(100000,20000);

fprintf('Iteration = 00000');
for i=1:20000
    fprintf('\b\b\b\b\b%5d',i)
    to_copy = fullfile(all_text_dir(randomnos(i)).folder,all_text_dir(randomnos(i)).name);
    to_save = fullfile(selected_save_dir,all_text_dir(randomnos(i)).name);
    copyfile(to_copy,to_save)
end
fprintf('\n')