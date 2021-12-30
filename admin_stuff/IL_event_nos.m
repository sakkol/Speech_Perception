function [trial_codes] = IL_event_nos(events, Number_of_sentences, Sentence_vs_Scrambled, Iso_A_chronous_Natural, Word_per_sentence)
% low level function for first grouping of events

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