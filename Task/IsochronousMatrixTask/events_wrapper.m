function [events_cell] = events_wrapper(SNR_input,default_selection)

if ~exist('default_selection','var') || isempty(default_selection)
    select_conditions
end



if any(selections.Electrical_stim)
    % Condition list to select from, for this block
    cond_list ={'Control condition','Control condition','Control condition',...
                'HG - short delay','HG - long delay',...
                'STG - short delay','STG - long delay',...
                'Epicranial - short delay','Epicranial - long delay'...
                'HG - sinewave - inphase','HG - sinewave - outphase',...
                'STG - sinewave - inphase','STG - sinewave - outphase',...
                'Epicranial - sinewave - inphase','Epicranial - sinewave - outphase',...
                'Epicranial - 100Hz - nonmodulated'};
    
    [cond_indx,~] = listdlg('ListString',cond_list);
    if sum(cond_indx)==0,error('Quiting, no input!'),end
    
    % get delays for each condition which needs delay to be specified
    if any(~ismember(cond_indx,[1 2 3]))
        all_delays = inputdlg(cond_list(cond_indx(~ismember(cond_indx,[1 2 3]))),'Please input delays in miliseconds',[1 55]);
    else
        all_delays=0;
    end
    if isempty(all_delays),error('Quiting, no input!'),end
    
    % Get delays
    use_cond_list = cond_list(cond_indx(~ismember(cond_indx,[1 2 3])));
    for i=1:length(cond_list(~ismember(cond_indx,[1 2 3])))
        condname_erase = erase(use_cond_list{i},{' ','-'});
        eval([condname_erase ' = ' all_delays{i} '/1000;'])
    end
end

%% create the trial order and loop to create trials

%% Randomization section
n_of_each_cond = 25;    % if you want to increase the number of trials per condition
if strcmp(cfg.language,'English')
    load('EnglishWordsInfo.mat','WordsInfo','avg_sig_pow','avg_noise_pow','words_table')
else
    load('SpanishWordsInfo.mat','WordsInfo','avg_sig_pow','avg_noise_pow','words_table')
end
all_sentence_list = readtable(fullfile(main_stim_loc, 'Sentence_to_use.xlsx'));

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



