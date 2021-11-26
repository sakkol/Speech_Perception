function [trial_codes] = IL_event_nos(events, Number_of_sentences, Sentence_vs_Scrambled, Iso_A_chronous_Natural, Word_per_sentence)
% to easily extract which trials you want
%% Preparation
if ~exist('Number_of_sentences','var') || strcmp(Number_of_sentences,'all')
    Number_of_sentences={'5sentences'};
end
if ~exist('Sentence_vs_Scrambled','var') || strcmp(Sentence_vs_Scrambled,'all')
    Sentence_vs_Scrambled = {'Sentence' ,'Scrambled'} ;
end
if ~exist('Iso_A_chronous_Natural','var') || strcmp(Iso_A_chronous_Natural,'all')
    Iso_A_chronous_Natural = {'iso','a'};
end
if ~exist('Word_per_sentence','var') || strcmp(Word_per_sentence,'all')
    Word_per_sentence = {'4-word','3-/5-word'};
end

%% Find them
trial_codes = [];
for t=1:size(events,1)
    if any(strcmp(events.cond_info{t}.Number_of_sentences,Number_of_sentences)) && ...
            any(strcmp(events.cond_info{t}.Sentence_vs_Scrambled,Sentence_vs_Scrambled)) && ...
            any(strcmp(events.cond_info{t}.Iso_A_chronous_Natural,Iso_A_chronous_Natural)) && ...
            any(strcmp(events.cond_info{t}.Word_per_sentence,Word_per_sentence))
        trial_codes(1,end+1) = t;
    end
end

end