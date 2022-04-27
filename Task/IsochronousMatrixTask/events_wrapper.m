function [events_table] = events_wrapper(SNR_input,language,LvsR,select_block_input)
% This is the wrapper for trial_creator. With simple inputs and needs, it
% returns everything you'd need in the block: stimuli, cfgs (to be able to
% use later), cond_info (to be able to use later)

%% Preparations
if ~exist('language','var') || isempty(language)
    language = [];
end
if ~exist('LvsR','var') || isempty(LvsR)
    LvsR = [];
end

% select the block and arrange several parameters
if ~exist('select_block_input','var') || isempty(select_block_input)
    [selections, select_block] = select_conditions_GUI([],language,LvsR);
else
    [selections, select_block] = select_conditions_GUI(select_block_input,language,LvsR);
end
language = selections.Language{1};
LvsR = selections.Left_vs_Right{1};
% SNR input
if strcmp(select_block,'iso_mat_24')
    if exist('SNR_input','var') && length(SNR_input) == 1
        SNR = SNR_input;
    elseif exist('SNR_input','var') && length(SNR_input) == 3 % run psignifit here
        [SNR] = estimate_threshold(3,SNR_input);
    elseif ~exist('SNR_input','var') || isempty(SNR_input)
        SNR = estimate_threshold(3); % 3 is the condition for 40-word list
    end
else
    SNR=1;
end

%% create the trial order and loop to create trials
if strcmp(language,'English')
    load('EnglishWordsInfo.mat','WordsInfo','avg_sig_pow','avg_noise_pow','words_table')
else
    load('SpanishWordsInfo.mat','WordsInfo','avg_sig_pow','avg_noise_pow','words_table')
end

% randomize word list
rng('shuffle')
all4wordsent = allcomb(words_table.Name,words_table.Verb,[words_table.Number;words_table.Adjective],words_table.Noun);
all4wordscramb = [all4wordsent(:,[1,3,2,4]);...
                  all4wordsent(:,[3,2,1,4]);...
                  all4wordsent(:,[4,1,2,3]);...
                  all4wordsent(:,[4,1,3,2]);...
                  all4wordsent(:,[2,1,4,3]);...
                  all4wordsent(:,[4,3,2,1])];
all5wordsent = allcomb(words_table.Name,words_table.Verb,words_table.Number,words_table.Adjective,words_table.Noun);
all5wordscramb = [all5wordsent(:,[1,3,4,5,2]);...
                  all5wordsent(:,[2,1,3,5,4]);...
                  all5wordsent(:,[4,3,2,1,5]);...
                  all5wordsent(:,[4,2,1,5,3]);...
                  all5wordsent(:,[3,4,2,1,5]);...
                  all5wordsent(:,[5,3,2,1,4])];
all3wordsent = all5wordsent(:,[1,2,5]);
all3wordscramb = [all5wordsent(:,[1,3,4]);...
                  all5wordsent(:,[2,1,3]);...
                  all5wordsent(:,[2,1,4]);...
                  all5wordsent(:,[5,4,3]);...
                  all5wordsent(:,[5,1,3]);...
                  all5wordsent(:,[4,2,3])];
% shuffle orders
all4wordsent = all4wordsent(randperm(length(all4wordsent)),:);
all4wordscramb = all4wordscramb(randperm(length(all4wordscramb)),:);
all5wordsent = all5wordsent(randperm(length(all5wordsent)),:);
all5wordscramb = all5wordscramb(randperm(length(all5wordscramb)),:);
all3wordsent = all3wordsent(randperm(length(all3wordsent)),:);
all3wordscramb = all3wordscramb(randperm(length(all3wordscramb)),:);

%% first create then shuffle the trials
w3sent=1;w4sent=1;w5sent=1;
w3scram=1;w4scram=1;w5scram=1;
mod35=1;

events_table=table;
events_table.trials(1:sum(str2double(selections.trial_no_per_block)))={''};
events_table.cfgs(1:sum(str2double(selections.trial_no_per_block)))={''};
events_table.cond_info(1:sum(str2double(selections.trial_no_per_block)))={''};
s=1;

fprintf('\n\t\t CREATING ALL EVENTS, THIS TAKES FEW SECONDS!\n\n')

for c = 1:size(selections,1)
    
    for t = 1:str2double(selections.trial_no_per_block{c})
        
        trial_words = {};
        
        cfg=[];
        cfg.language=selections.Language{c};
        if strcmp(selections.Clean_vs_In_noise{c},'clean')
            cfg.noise = 'silence';
        elseif strcmp(selections.Clean_vs_In_noise{c},'in-noise')
            cfg.noise='pink';
        end
        cfg.SNR=SNR;
        cfg.LvsR=LvsR;
        cfg.part1.length=0.5;
        
        if strcmp(selections.Number_of_sentences{c},'1attention-1sentence')
            cfg.part2.chronicity='iso';
            cfg.part2.frequency=selections.Frequency(c);
            cfg.part2.word1={'now'};
            cfg.part2.word2={'catch'};
            cfg.part2.word3={'these'};
            cfg.part2.word4={'words'};
            cfg.part2.estim=0;
            
            trial_words = {'now','catch','these','words'};
            
            if strcmp(selections.Word_per_sentence{c},'4-word') && strcmp(selections.Sentence_vs_Scrambled{c},'Sentence')
                curr_sent = all4wordsent(w4sent,:);w4sent=w4sent+1;
            elseif strcmp(selections.Word_per_sentence{c},'4-word') && strcmp(selections.Sentence_vs_Scrambled{c},'Scrambled')
                curr_sent = all4wordscramb(w4scram,:);w4scram=w4scram+1;
            elseif strcmp(selections.Word_per_sentence{c},'3-/5-word') && strcmp(selections.Sentence_vs_Scrambled{c},'Scrambled')
                if mod(mod35,2) % change from 3 to 5 word in each loop, so that there will be equal no of them
                    curr_sent = all3wordscramb(w3scram,:);w3scram=w3scram+1;
                else
                    curr_sent = all5wordscramb(w5scram,:);w5scram=w5scram+1;
                end
                mod35=mod35+1;
            elseif strcmp(selections.Word_per_sentence{c},'3-/5-word') && strcmp(selections.Sentence_vs_Scrambled{c},'Sentence')
                if mod(mod35,2) % change from 3 to 5 word in each loop, so that there will be equal no of them
                    curr_sent = all3wordsent(w3sent,:);w3sent=w3sent+1;
                else
                    curr_sent = all5wordsent(w5sent,:);w5sent=w5sent+1;
                end
                mod35=mod35+1;
            elseif strcmp(selections.Word_per_sentence{c},'3-word') && strcmp(selections.Sentence_vs_Scrambled{c},'Scrambled')
                curr_sent = all3wordscramb(w3scram,:);w3scram=w3scram+1;
            elseif strcmp(selections.Word_per_sentence{c},'3-word') && strcmp(selections.Sentence_vs_Scrambled{c},'Sentence')
                curr_sent = all3wordsent(w3sent,:);w3sent=w3sent+1;
            elseif strcmp(selections.Word_per_sentence{c},'5-word') && strcmp(selections.Sentence_vs_Scrambled{c},'Scrambled')
                curr_sent = all5wordscramb(w5scram,:);w5scram=w5scram+1;
            elseif strcmp(selections.Word_per_sentence{c},'5-word') && strcmp(selections.Sentence_vs_Scrambled{c},'Sentence')
                curr_sent = all5wordsent(w5sent,:);w5sent=w5sent+1;
            end
            
            for w=1:length(curr_sent)
                cfg.part3.(['word' num2str(w)]) = curr_sent(w);
            end
            trial_words = [trial_words,curr_sent];
            
            cfg.part3.chronicity=selections.Iso_A_chronous_Natural{c};
            cfg.part3.frequency=selections.Frequency(c);
            cfg.part3.estim=~strcmp(selections.Electrical_stim{c},'Control no-stim condition');
            cfg.part3.delay = (selections.E_stim_delay(c))/1000;
            
            cfg.part4 = cfg.part1;
            cfg.part4.length = 2;
            
        elseif strcmp(selections.Number_of_sentences{c},'3sentences') || strcmp(selections.Number_of_sentences{c},'5sentences')
            pss = str2double(selections.Number_of_sentences{c}(1))+1;
            
            % loop sentences (equals to "parts" field)
            p=2;
            while p<pss+1
                
                % select words from randomly generated list
                if strcmp(selections.Word_per_sentence{c},'4-word') && strcmp(selections.Sentence_vs_Scrambled{c},'Sentence')
                    curr_sent = all4wordsent(w4sent,:);w4sent=w4sent+1;
                elseif strcmp(selections.Word_per_sentence{c},'4-word') && strcmp(selections.Sentence_vs_Scrambled{c},'Scrambled')
                    curr_sent = all4wordscramb(w4scram,:);w4scram=w4scram+1;
                elseif strcmp(selections.Word_per_sentence{c},'3-/5-word') && strcmp(selections.Sentence_vs_Scrambled{c},'Scrambled')
                    if mod(mod35,2)
                        curr_sent = all3wordscramb(w3scram,:);w3scram=w3scram+1;
                    else
                        curr_sent = all5wordscramb(w5scram,:);w5scram=w5scram+1;
                    end
                    mod35=mod35+1;
                elseif strcmp(selections.Word_per_sentence{c},'3-/5-word') && strcmp(selections.Sentence_vs_Scrambled{c},'Sentence')
                    if mod(mod35,2)
                        curr_sent = all3wordsent(w3sent,:);w3sent=w3sent+1;
                    else
                        curr_sent = all5wordsent(w5sent,:);w5sent=w5sent+1;
                    end
                    mod35=mod35+1;
                elseif strcmp(selections.Word_per_sentence{c},'3-word') && strcmp(selections.Sentence_vs_Scrambled{c},'Scrambled')
                    curr_sent = all3wordscramb(w3scram,:);w3scram=w3scram+1;
                elseif strcmp(selections.Word_per_sentence{c},'3-word') && strcmp(selections.Sentence_vs_Scrambled{c},'Sentence')
                    curr_sent = all3wordsent(w3sent,:);w3sent=w3sent+1;
                elseif strcmp(selections.Word_per_sentence{c},'5-word') && strcmp(selections.Sentence_vs_Scrambled{c},'Scrambled')
                    curr_sent = all5wordscramb(w5scram,:);w5scram=w5scram+1;
                elseif strcmp(selections.Word_per_sentence{c},'5-word') && strcmp(selections.Sentence_vs_Scrambled{c},'Sentence')
                    curr_sent = all5wordsent(w5sent,:);w5sent=w5sent+1;
                end
                
                % if there are words used in this trial already, change it
                if strcmp(select_block,'FreeRecall_3sent')
                    if any(ismember(trial_words,curr_sent))
                        continue
                    end
                end
                
                % select general parameters
                cfg.(['part' num2str(p)]).chronicity=selections.Iso_A_chronous_Natural;
                cfg.(['part' num2str(p)]).frequency=selections.Frequency;
                for w=1:length(curr_sent)
                    cfg.(['part' num2str(p)]).(['word' num2str(w)]) = curr_sent(w);
                end
                trial_words = [trial_words,curr_sent];
                cfg.(['part' num2str(p)]).chronicity=selections.Iso_A_chronous_Natural{c};
                cfg.(['part' num2str(p)]).frequency=selections.Frequency(c);
                
                p = p+1;
            end
            
            cfg.(['part' num2str(p+1)]) = cfg.part1;
            
        end
        
        [events_table.trials{s}, events_table.cfgs{s}] = trial_creator(cfg);
        events_table.cond_info{s} = selections(c,:);
        events_table.trial_words{s} = trial_words;
        s=s+1;
        
    end
end


%% Randomization of trials
events_table = events_table(randperm(size(events_table,1)),:);

end
