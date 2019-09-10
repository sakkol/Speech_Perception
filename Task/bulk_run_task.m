%% Create many many stimuli :)

cd /home/sakkol/Documents/Speech_Perception_stim/02nd_Generation/GoogleTTS_M_1.7
speech_files = dir('*.wav');
stim_save_dir = '/home/sakkol/Documents/Speech_Perception_stim/02nd_Generation/Created_Stim_1.7';

cfg=[];
cfg.SNR = -4;
cfg.prespeech.part1.length= 0.5;
cfg.prespeech.part1.noise = 'pink';
cfg.prespeech.part2.noise = 'pink';
cfg.prespeech.part2.signal = 'C:\Users\user\Desktop\PROJECT MANAGEMENT\PhD\TASK\Matrix_Sentences\03rd_Generation\Pre-stim-Attention-comma-M-Rate0.9.wav';

cfg.speech.noise = 'pink';
cfg.postspeech.part1.length=2;
cfg.postspeech.part1.noise = 'pink';
cfg.LvsR = 'L';

for s=1:length(speech_files)
    cfg.speech.file = [speech_files(s).folder filesep speech_files(s).name];
    
    curr_sentence = erase(speech_files(s).name,'.wav');
    cfg.stim_save_filename = [stim_save_dir filesep curr_sentence 'LpatientRtdt.wav'];
    cfg.envelope_save_filename = [stim_save_dir filesep curr_sentence '.mat'];
    cfg.plot_save_filename = [stim_save_dir filesep curr_sentence '.jpg'];
    
    [stimulus,envelope]=stimuli_creator(cfg);
    
end


%% Create with v2 convention
main_stim_loc = 'C:\Users\user\Desktop\PROJECT MANAGEMENT\PhD\TASK\Matrix_Sentences\03rd_Generation';
speech_rate = 0.9;

cfg=[];
cfg.SNR = -4;
cfg.prespeech.part1.length= 0.5;
cfg.prespeech.part1.noise = 'pink';

cfg.prespeech.part2.noise = 'pink';
cfg.prespeech.part2.signal = find_sentence('Pre-stim-Attention-comma-Rate0.9',main_stim_loc,0.9);

cfg.speech.noise = 'pink';
cfg.speech.file = find_sentence(curr_sentence,main_stim_loc,speech_rate);

cfg.postspeech.part1.length=2;
cfg.postspeech.part1.noise = 'pink';

cfg.LvsR = 'L';
cfg.delay = 1;

%     cfg.stim_save_filename = [stim_save_dir filesep curr_sentence 'LpatientRtdt.wav'];
%     cfg.envelope_save_filename = [stim_save_dir filesep curr_sentence '.mat'];
    cfg.plot_save_filename = ['C:\Users\user\Desktop\PROJECT MANAGEMENT\PhD\TASK\Matrix_Sentences\03rd_Generation\Sentences_Rate0.9' filesep curr_sentence '.jpg'];

[stimulus,envelope]=stim_creatorv2(cfg);
    
%% Create stim for thresholding

main_stim_loc = '/home/sakkol/Documents/Speech_Perception_stim/4th_Generation';
speech_rate = 1.3;

all_sentences = readtable('/home/sakkol/Documents/Speech_Perception_stim/4th_Generation/Sentence_to_use.xlsx','Sheet','Threshold');

threshold_sentences = all_sentences.Sentence;

stim_save_dir = fullfile(main_stim_loc,'threshold_sounds_fast');

for i=21:30
cfg=[];
cfg.SNR = 2;
cfg.prespeech.part1.length= 0.5;
cfg.prespeech.part1.noise = 'pink';

cfg.prespeech.part2.noise = 'pink';
cfg.prespeech.part2.signal = find_sentence('Pre-stim-Attention-comma',main_stim_loc,speech_rate);

cfg.speech.noise = 'pink';
cfg.speech.file = find_sentence(threshold_sentences{i},main_stim_loc,speech_rate);

cfg.postspeech.part1.length=2;
cfg.postspeech.part1.noise = 'pink';

cfg.LvsR = 'both';
% cfg.delay = 1;

    cfg.stim_save_filename = [stim_save_dir filesep num2str(cfg.SNR) 'SNR_' threshold_sentences{i} '_bilateral.wav'];
%     cfg.envelope_save_filename = [stim_save_dir filesep curr_sentence '.mat'];
    cfg.plot_save_filename = [stim_save_dir filesep num2str(cfg.SNR) 'SNR_' threshold_sentences{i} '_bilateral.jpg'];

[stimulus,envelope]=stimuli_creator(cfg);
end

