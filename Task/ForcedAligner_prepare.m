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