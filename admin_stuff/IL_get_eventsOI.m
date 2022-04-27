function [eventsOI] = IL_get_eventsOI(events, Sentence_vs_Scrambled, Iso_A_chronous_Natural, Word_per_sentence,accuracy)
% a different way to extract events: that different way is IL_event_types.m

[event_types] = IL_event_types(events);

eventsOI = find(strcmp(event_types.sent_scr,Sentence_vs_Scrambled) & ...
    strcmp(event_types.iso_a,Iso_A_chronous_Natural) & ...
    event_types.word_count==Word_per_sentence & ...
    ismember(events.accuracy,accuracy))';
end