%% For English:
% import excel to read sentences
sentences = readtable('/home/sakkol/Documents/BACKUPS/Speech_Perception_BACKUP(2019.09.11)/SpeechPerception_Task(backup-11.09.2019)/SpeechPerception_Task/main_stim_loc/Sentence_to_use.xlsx');
sentences = sentences.Sentence;

% where texts are:
txtdir = '/home/sakkol/Documents/BACKUPS/Speech_Perception_BACKUP(2019.09.11)/Speech_Perception_stim/4th_Generation/Selected_Texts_4thGen';

wavdir = '/home/sakkol/Documents/BACKUPS/Speech_Perception_BACKUP(2019.09.11)/SpeechPerception_Task(backup-11.09.2019)/SpeechPerception_Task/main_stim_loc/Sentences_Rate0.9';
newtxtdir = '/home/sakkol/Documents/Forced_Alignment/FORCE/Slow-English/txt';
newwavdir = '/home/sakkol/Documents/Forced_Alignment/FORCE/Slow-English/wav';

wavdir = '/home/sakkol/Documents/BACKUPS/Speech_Perception_BACKUP(2019.09.11)/SpeechPerception_Task(backup-11.09.2019)/SpeechPerception_Task/main_stim_loc/Sentences_Rate1.3';
newtxtdir = '/home/sakkol/Documents/Forced_Alignment/FORCE/Fast-English/txt';
newwavdir = '/home/sakkol/Documents/Forced_Alignment/FORCE/Fast-English/wav';



% Loop for each sentence
for s = 1:length(sentences)

% find and copy sentence txt -no dash
txtname = [replace(sentences{s},' ','_') '.txt'];
textitself = fileread(fullfile(txtdir,txtname));
textitself = replace(textitself,'eigth','eight');

fid = fopen(fullfile(newtxtdir,txtname),'w');
fprintf(fid,'%s',upper(textitself));
fclose(fid);

% find and copy sentence wav -no dash
wavname = [replace(sentences{s},' ','_') '-M.wav'];
wavname_new = erase(wavname,'-M');
copyfile(fullfile(wavdir,wavname),fullfile(newwavdir,wavname_new))

end

%% For Spanish:
% import excel to read sentences
sentences = readtable('/home/sakkol/Documents/BACKUPS/Spanish_SpeechPerception_Task_BACKUP_31.10.2019/main_stim_loc/sentence_to_use.xlsx');
sentences = sentences.Sentence;

% where texts are:
txtdir = '/home/sakkol/Documents/BACKUPS/Spanish_SpeechPerception_Task_BACKUP_31.10.2019/Preparations/Version_4/All_Texts';

wavdir = '/home/sakkol/Documents/BACKUPS/Spanish_SpeechPerception_Task_BACKUP_31.10.2019/main_stim_loc/Sentences_Rate0.9';
newtxtdir = '/home/sakkol/Documents/Forced_Alignment/FORCE/Slow-Spanish/txt';
newwavdir = '/home/sakkol/Documents/Forced_Alignment/FORCE/Slow-Spanish/wav';

wavdir = '/home/sakkol/Documents/BACKUPS/Spanish_SpeechPerception_Task_BACKUP_31.10.2019/main_stim_loc/Sentences_Rate1.3';
newtxtdir = '/home/sakkol/Documents/Forced_Alignment/FORCE/Fast-Spanish/txt';
newwavdir = '/home/sakkol/Documents/Forced_Alignment/FORCE/Fast-Spanish/wav';



% Loop for each sentence
for s = 1:length(sentences)

% find and copy sentence txt -no dash
txtname = [replace(sentences{s},' ','_') '.txt'];
fID = fopen(fullfile(txtdir,txtname), 'r','n','latin-1');
textitself = fread(fID);
fclose(fID);
textitself = native2unicode(textitself, 'latin-1')';

fid = fopen(fullfile(newtxtdir,txtname),'w');
fprintf(fid,'%s',upper(textitself));
fclose(fid);

% find and copy sentence wav -no dash
wavname = [replace(sentences{s},' ','_') '-F.wav'];
wavname_new = erase(wavname,'-F');
copyfile(fullfile(wavdir,wavname),fullfile(newwavdir,wavname_new))

end


%% Now it's time to combine all the output (WIP part)
examplesent = '/home/sakkol/Documents/Forced_Alignment/FORCE/Fast-English/TextGrid/Alan_kept_two_pretty_sofas (copy).TextGrid';

fid = fopen(examplesent);
clear A
A{1,1} = fgetl(fid);
n=1;
while ischar(A{n})
    A{n+1,1} = fgetl(fid);
    n=n+1;
end
fclose(fid);

innttierlines = find(strcmp(A,'"IntervalTier"'));

% phoneme infos
pho_lines = A(innttierlines(1)+5:innttierlines(2)-1);
% word infos
w_lines = A(innttierlines(2)+5:end-1);

phoneme_info = table;i=1;
clear phoneme onset offset
for p = 3:3:length(pho_lines)
    if strcmpi(pho_lines{p},'"sp"'),continue,end
    phoneme{i,1} = erase(pho_lines{p},'"');
    onset(i,1) = str2double(pho_lines{p-2});
    offset(i,1) = str2double(pho_lines{p-1});
    i=i+1;
end
phoneme_info.phoneme = phoneme;
phoneme_info.onset = onset;
phoneme_info.offset = offset;

word_info = table;i=1;
clear word onset offset
for w = 3:3:length(w_lines)
    if strcmpi(w_lines{w},'"sp"'),continue,end
    word{i,1} = erase(w_lines{w},'"');
    onset(i,1) = str2double(w_lines{w-2});
    offset(i,1) = str2double(w_lines{w-1});
    i=i+1;
end
word_info.word = word;
word_info.onset = onset;
word_info.offset = offset;

%% Let's do it
FORCE_dir = '/home/sakkol/Documents/Forced_Alignment/FORCE';
lang = {'English','Spanish'};
Rate = {'Fast','Slow'};
onset_table = struct;

for l = 1% There is only English currently, Spanish doesn't have forced alignment
    for r = 1:2
        tgrids = dir(fullfile(FORCE_dir,[Rate{r} '-' lang{l}],'TextGrid','*.TextGrid'));
        
        for t = 1:length(tgrids)
            curr_sent = erase(tgrids(t).name,'.TextGrid');
            
            fid = fopen(fullfile(tgrids(t).folder,tgrids(t).name));
            clear A
            A{1,1} = fgetl(fid);
            n=1;
            while ischar(A{n})
                A{n+1,1} = fgetl(fid);
                n=n+1;
            end
            fclose(fid);
            
            innttierlines = find(strcmp(A,'"IntervalTier"'));
            
            % phoneme lines
            pho_lines = A(innttierlines(1)+5:innttierlines(2)-1);
            % word lines
            w_lines = A(innttierlines(2)+5:end-1);
            
            % Prepare the table containing info for phonemes and words
            phoneme_info = table;i=1;
            clear phoneme onset offset
            for p = 3:3:length(pho_lines)
                if strcmpi(pho_lines{p},'"sp"'),continue,end
                phoneme{i,1} = erase(pho_lines{p},'"');
                onset(i,1) = str2double(pho_lines{p-2});
                offset(i,1) = str2double(pho_lines{p-1});
                i=i+1;
            end
            phoneme_info.phoneme = phoneme;
            phoneme_info.onset = onset;
            phoneme_info.offset = offset;
            
            word_info = table;i=1;
            clear word onset offset
            for w = 3:3:length(w_lines)
                if strcmpi(w_lines{w},'"sp"'),continue,end
                word{i,1} = erase(w_lines{w},'"');
                onset(i,1) = str2double(w_lines{w-2});
                offset(i,1) = str2double(w_lines{w-1});
                i=i+1;
            end
            word_info.word = word;
            word_info.onset = onset;
            word_info.offset = offset;
            
            % Gather all info
            onset_table(end+1).lang = lang{l};
            onset_table(end).rate = Rate{l};
            onset_table(end).sentence = curr_sent;
            onset_table(end).phoneme_info = phoneme_info;
            onset_table(end).word_info = word_info;
        end
    end
end

onset_table(1)=[];
onset_table = struct2table(onset_table);
save(fullfile(FORCE_dir,'English_onset_info.mat'), 'onset_table')

%% Now add information about the syllables

wsp_table = readtable('/home/sakkol/Documents/Forced_Alignment/word_phoneme_syllable.xlsx');
onset_table = onset_table_orig;
onset_table = table2struct(onset_table);

for i=1:size(onset_table,1)
    curr_sentence = onset_table(i).sentence;
    sent_words = split(curr_sentence,' ');
    pids=1;
    syllable_info = struct;all_info = struct;
    for w=1:length(sent_words)
        
        curr_wsp = wsp_table(strcmpi(sent_words{w},wsp_table.word),:);
        curr_word_sylls = split(curr_wsp.syllables,',');
        curr_word_phos = split(curr_wsp.phoneme,',');
        curr_word_syll_ids = str2double(split(curr_wsp.syllable_ids,','));
        curr_pho_info = onset_table(i).phoneme_info(pids:pids+length(curr_word_syll_ids)-1,:);
        pids = pids+length(curr_word_syll_ids);
        
        for s=1:length(curr_word_sylls)
            syllable_info(end+1).syllable = curr_word_sylls{s};
            syllable_info(end).syllable_id = s;
            syllable_info(end).onset = curr_pho_info.onset(min(find(curr_word_syll_ids==s)));
            syllable_info(end).offset = curr_pho_info.offset(max(find(curr_word_syll_ids==s)));
            
            syllable_info(end).word = curr_wsp.word;
            syllable_info(end).word_id = w;
            syllable_info(end).phonemes = curr_word_phos(curr_word_syll_ids==s);
            syllable_info(end).phoneme_info = curr_pho_info(curr_word_syll_ids==s,:);
            
            % all_info
            for p = 1:size(syllable_info(end).phonemes,1)
                all_info(end+1).word = curr_wsp.word;
                all_info(end).word_id = w;
                all_info(end).word_onset = onset_table(i).word_info.onset(w);
                all_info(end).word_offset = onset_table(i).word_info.offset(w);
                
                all_info(end).syllable = curr_word_sylls{s};
                all_info(end).syllable_id = s;
                all_info(end).onset = curr_pho_info.onset(min(find(curr_word_syll_ids==s)));
                all_info(end).offset = curr_pho_info.offset(max(find(curr_word_syll_ids==s)));
                all_info(end).phos_of_syll = curr_word_phos(curr_word_syll_ids==s);
                all_info(end).phos_info_of_syll = curr_pho_info(curr_word_syll_ids==s,:);
                
                all_info(end).phoneme = syllable_info(end).phonemes(p);
                all_info(end).phoneme_id = p;
                all_info(end).phoneme_onset = syllable_info(end).phoneme_info.onset(p);
                all_info(end).phoneme_offset = syllable_info(end).phoneme_info.offset(p);
                
            end
            
            
        end
    end
    syllable_info(1)=[];
    syllable_info = struct2table(syllable_info);
    onset_table(i).syllable_info = syllable_info;
    
    all_info(1)=[];
    all_info = struct2table(all_info);
    onset_table(i).all_info = all_info;
    
end

all_info_table = struct2table(onset_table);
save(fullfile(FORCE_dir,'English_all_info.mat'), 'all_info_table')