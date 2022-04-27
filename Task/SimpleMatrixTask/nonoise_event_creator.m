function events_cell = nonoise_event_creator(main_stim_loc,slowVSfast)
% This function is a derivative of event_creator.m function, specifically
% to be used if sentences are going to be presented without any noise.
% July 2020, Serdar Akkol, HBML.

%% INPUTS
% select speech rate and ear to present speech
if nargin > 1 && exist('slowVSfast','var')
    if slowVSfast == 1
        speech_rate_input = 'Slow (x0.9 = 2.95Hz)';
    elseif slowVSfast == 2
        speech_rate_input = 'Fast (x1.3 = 4.4Hz)';
    end
else
    speech_rate_input = questdlg('Slow is default!', ...
        'Please select speech rate:', ...
        'Slow (x0.9 = 2.95Hz)','Fast (x1.3 = 4.4Hz)','Slow (x0.9 = 2.95Hz)');
end
% set speech rate
switch speech_rate_input
    case 'Slow (x0.9 = 2.95Hz)'
        speech_rate = '0.9';slowVSfast=1;
    case 'Fast (x1.3 = 4.4Hz)'
        speech_rate = '1.3';slowVSfast=2;
    case ''
        error('Quiting, no input for speech rate!')
end

% select ear to present speech
LvsR_input = questdlg('You need to select one!', ...
    'Which ear to present speech?', ...
    'Present to Left','Present to Right','Quit','Quit');
if strcmp(LvsR_input,'Quit')
    error('You wanted to quit!')
elseif strcmp(LvsR_input,'Present to Left')
    LvsR = 'L';
elseif strcmp(LvsR_input,'Present to Right')
    LvsR = 'R';
end
clear LvsR_input

%%
%% Randomization section
n_of_each_cond = 20;    % if you want to increase the number of trials per condition
all_sentence_list = readtable(fullfile(main_stim_loc, 'Sentence_to_use.xlsx'));
cond_indx=1;
cond_list={'No Noise Stimuli'};

% random selection from main sentence repository
if contains(main_stim_loc,'Spanish')
    random_no = randperm(350);
else
    random_no = randperm(367);
end
random_no = random_no(1:n_of_each_cond*length(cond_indx));
random_sentences = all_sentence_list.Sentence(random_no);
random_codes = all_sentence_list.Code(random_no);

% randomize conditions in block
pool  = repelem(cond_indx, n_of_each_cond);
index = randperm(numel(pool), n_of_each_cond*length(cond_indx));
random_conds_code = pool(index)';
random_conds_name = cond_list(random_conds_code)';

% correct 2 and 3 to 1, because they are also Control condition (names are correct)
random_conds_code(ismember(random_conds_code,[2 3])) = 1;

%% Now create stimuli and put it into event_cell
% common parameters
cfg = [];

switch speech_rate_input
    case 'Slow (x0.9 = 2.95Hz)'
        cfg.frequency = 2.95; % if needed in sinewave conditions
    case 'Fast (x1.3 = 4.4Hz)'
        cfg.frequency = 4.4;
    case ''
        error('Quiting, no input!')
end

cfg.prespeech.part1.length= 0.5;
cfg.prespeech.part1.noise = 'silence';

cfg.prespeech.part2.noise = 'silence';

if contains(main_stim_loc,'Spanish')
    cfg.prespeech.part2.signal = find_sentence('Por_favor,_presta_atención_y_recuerda_esta_oración',main_stim_loc,speech_rate);
else
    cfg.prespeech.part2.signal = find_sentence('Pre-stim-Attention-comma',main_stim_loc,speech_rate);
end

cfg.speech.noise = 'silence';

cfg.postspeech.part1.length=1;
cfg.postspeech.part1.noise = 'silence';

cfg.LvsR = LvsR;
cfg.nostim = 1;

% Prepare events_cell
events_cell=cell(length(random_sentences),6);

for cond_i = 1:length(random_sentences)
    
    events_cell{cond_i,1} = random_codes(cond_i);
    events_cell(cond_i,2) = random_sentences(cond_i);
    events_cell{cond_i,3} = random_conds_code(cond_i);
    events_cell(cond_i,4) = random_conds_name(cond_i);
    
    % sentence gets ready
    cfg.speech.file = find_sentence(events_cell{cond_i,2},main_stim_loc,speech_rate);
    
    % speech rate: need to do this in each loop so that if there is 100Hz,
    % it will be renewed in the beginning of the loop

    
    [stimulus,~]=stim_creatorv2(cfg);
    events_cell{cond_i,5} = stimulus;
    events_cell{cond_i,6} = cfg;
    
end

end