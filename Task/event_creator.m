function events_cell = event_creator(main_stim_loc)
% Input is only directory of where Sentences_Rate are. (eg. /home/sakkol/Documents/Speech_Perception_stim/4th_Generation)
% events_cell, everything ready output:
% first column is Code of sentence; second is Sentence itself; third and
% fourth are Condition code and name; Fifth is stimulus which will be given
% and Fifth is cfg input to stim_creatorv2 (saving this to be on the safer side.)

%% INPUTS
% select speech rate
speech_rate_input = questdlg('Slow is default!', ...
	'Please select speech rate:', ...
	'Slow (x0.9 = 2.95Hz)','Fast (x1.3 = 4.4Hz)','Slow (x0.9 = 2.95Hz)');

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

% SNR input
SNR = inputdlg('Please write SNR in dB!');
SNR = str2double(SNR{1});

% Condition list to select from, for this block
cond_list = {'Control condition','HG - short delay','HG - long delay',...
    'STG - short delay','STG - long delay',...
    'Epicranial - short delay','Epicranial - long delay'...
    'HG - sinewave - inphase','HG - sinewave - outphase',...
    'STG - sinewave - inphase','STG - sinewave - outphase'};

[cond_indx,~] = listdlg('ListString',cond_list);

% get delays for each condition which has delay
if any(cond_indx~=1)
    all_delays = inputdlg(cond_list(cond_indx(cond_indx~=1)),'Please inout delays in miliseconds',[1 55]);
end
if isempty(all_delays),error('Quiting, no input!'),end

% Get delays
use_cond_list = cond_list(cond_indx(cond_indx~=1));
for i=1:length(cond_list(cond_indx~=1))
    condname_erase = erase(use_cond_list{i},{' ','-'});
    eval([condname_erase ' = ' all_delays{i} '/1000;'])
end

%% Randomization section
n_of_each_cond = 20;    % if you want to increase the number of trials per condition
all_sentence_list = readtable([main_stim_loc filesep 'Sentence_to_use.xlsx']);

% random selection from main sentence repository
random_no = randi([1 367],n_of_each_cond*length(cond_indx),1);
random_sentences = all_sentence_list.Sentence(random_no);
random_codes = all_sentence_list.Code(random_no);

% randomize conditions in block
pool  = repelem(cond_indx, n_of_each_cond);
index = randperm(numel(pool), n_of_each_cond*length(cond_indx));
random_conds_code = pool(index)';
random_conds_name = cond_list(random_conds_code)';

%% Now create stimuli and put it into event_cell
% common parameters
cfg = [];
% speech rate
switch speech_rate_input
    case 'Slow (x0.9 = 2.95Hz)'
        cfg.frequency = 2.95; % if needed in sinewave conditions
        speech_rate = '0.9';
    case 'Fast (x1.3 = 4.4Hz)'
        cfg.frequency = 4.4;
        speech_rate = '1.3';
    case ''
        error('Quiting, no input!')
end

cfg.SNR = SNR;
cfg.prespeech.part1.length= 0.5;
cfg.prespeech.part1.noise = 'pink';

cfg.prespeech.part2.noise = 'pink';
cfg.prespeech.part2.signal = find_sentence('Pre-stim-Attention-comma',main_stim_loc,0.9);

cfg.speech.noise = 'pink';

cfg.postspeech.part1.length=2;
cfg.postspeech.part1.noise = 'pink';

cfg.LvsR = LvsR;

% Prepare events_cell
events_cell=cell(length(random_sentences),6);

for cond_i = 1:length(random_sentences)
    
    events_cell{cond_i,1} = random_codes(cond_i);
    events_cell(cond_i,2) = random_sentences(cond_i);
    events_cell{cond_i,3} = random_conds_code(cond_i);
    events_cell(cond_i,4) = random_conds_name(cond_i);
    
    % sentence gets ready
    cfg.speech.file = find_sentence(events_cell{cond_i,2},main_stim_loc,speech_rate);
    
    % clean cfg from previous loop
    if isfield(cfg,'delay'),cfg=rmfield(cfg,'delay');end
    if isfield(cfg,'nostim'),cfg=rmfield(cfg,'nostim');end
    if isfield(cfg,'phase'),cfg=rmfield(cfg,'phase');end
    
    % create cfg: input into stim_creatorv2
    if strcmp(random_conds_name{cond_i},'Control condition')
        cfg.nostim = 1;

    elseif strcmp(random_conds_name{cond_i},'HG - short delay')
        cfg.delay = HGshortdelay;
        
    elseif strcmp(random_conds_name{cond_i},'HG - long delay')
        cfg.delay = HGlongdelay;
        
    elseif strcmp(random_conds_name{cond_i},'STG - short delay')
        cfg.delay = STGshortdelay;
        
    elseif strcmp(random_conds_name{cond_i},'STG - long delay')
        cfg.delay = STGlongdelay;
        
    elseif strcmp(random_conds_name{cond_i},'Epicranial - short delay')
        cfg.delay = Epicranialshortdelay;
        
    elseif strcmp(random_conds_name{cond_i},'Epicranial - long delay')
        cfg.delay = Epicraniallongdelay;
        
    elseif strcmp(random_conds_name{cond_i},'HG - sinewave - inphase')
        cfg.delay = 'attention';
        cfg.phase = 'in';
        
    elseif strcmp(random_conds_name{cond_i},'HG - sinewave - outphase')
        cfg.delay = 'attention';
        cfg.phase = 'out';
        
    elseif strcmp(random_conds_name{cond_i},'STG - sinewave - inphase')
        cfg.delay = 'attention';
        cfg.phase = 'in';
        
    elseif strcmp(random_conds_name{cond_i},'STG - sinewave - outphase')
        cfg.delay = 'attention';
        cfg.phase = 'out';
    
    else
        error('There is a mismatch between cond_list and nested if names')
    end
    
    [stimulus,~]=stim_creatorv2(cfg);
    events_cell{cond_i,5} = stimulus;
    events_cell{cond_i,6} = cfg;
    
end


end