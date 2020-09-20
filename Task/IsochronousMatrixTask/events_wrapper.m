function [events_cell] = events_wrapper(SNR_input,default_selection)

if ~exist('default_selection','var') || isempty(default_selection)
    select_conditions
else
    selections = default_conds(default_selection);
end













% SNR input
if nargin > 2 && exist('SNR_input','var') && length(SNR_input) == 1
    SNR = SNR_input;
elseif nargin > 2 && exist('SNR_input','var') && length(SNR_input) == 3 % run psignifit here
    [SNR] = estimate_threshold(slowVSfast,SNR_input);
elseif ~exist('SNR_input','var') || isempty(SNR_input)
    SNR = estimate_threshold(slowVSfast);
end


%% create the trial order and loop to create trials

%% Randomization section
if strcmp(cfg.language,'English')
    load('EnglishWordsInfo.mat','WordsInfo','avg_sig_pow','avg_noise_pow','words_table')
else
    load('SpanishWordsInfo.mat','WordsInfo','avg_sig_pow','avg_noise_pow','words_table')
end
all_sentence_list = readtable(fullfile(main_stim_loc, 'Sentence_to_use.xlsx'));
n_of_each_cond = 25;    % if you want to increase the number of trials per condition

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


%% switch the possible conditions



