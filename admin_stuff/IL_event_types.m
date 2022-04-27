function [event_types] = IL_event_types(events)
% a different way to extract events: that different way is IL_get_eventsOI.m

% get the trial no for each condition
iso_4_sentence = IL_event_nos(events, '5sentences', 'Sentence', 'iso', '4-word');
iso_4_scrambled = IL_event_nos(events, '5sentences', 'Scrambled', 'iso', '4-word');
a_4_sentence = IL_event_nos(events, '5sentences', 'Sentence', 'a', '4-word');
a_4_scrambled = IL_event_nos(events, '5sentences', 'Scrambled', 'a', '4-word');
iso_35_sentence = IL_event_nos(events, '5sentences', 'Sentence', 'iso', '3-/5-word');
iso_35_scrambled = IL_event_nos(events, '5sentences', 'Scrambled', 'iso', '3-/5-word');
a_35_sentence = IL_event_nos(events, '5sentences', 'Sentence', 'a', '3-/5-word');
a_35_scrambled = IL_event_nos(events, '5sentences', 'Scrambled', 'a', '3-/5-word');
% 3 or 5 word sentences
nnnnnw=[];
for i=1:height(events)
    nnnnn=fieldnames(events.cfgs{i}.part2);
    nnnnnw(i,1)=str2double(nnnnn{end}(end));
end
nnnnnw3=find(nnnnnw==3);nnnnnw5=find(nnnnnw==5);
iso_3_sentence=intersect(iso_35_sentence,nnnnnw3);
iso_5_sentence=intersect(iso_35_sentence,nnnnnw5);
a_3_sentence=intersect(a_35_sentence,nnnnnw3);
a_5_sentence=intersect(a_35_sentence,nnnnnw5);
iso_3_scrambled=intersect(iso_35_scrambled,nnnnnw3);
iso_5_scrambled=intersect(iso_35_scrambled,nnnnnw5);
a_3_scrambled=intersect(a_35_scrambled,nnnnnw3);
a_5_scrambled=intersect(a_35_scrambled,nnnnnw5);

%% put them into a table

mycell=cell(height(events),3);
mycell(iso_4_sentence,:) = repmat([{4},{'iso'},{'sentence'}],length(iso_4_sentence),1);
mycell(iso_4_scrambled,:) = repmat([{4},{'iso'},{'scrambled'}],length(iso_4_scrambled),1);
mycell(a_4_sentence,:) = repmat([{4},{'a'},{'sentence'}],length(a_4_sentence),1);
mycell(a_4_scrambled,:) = repmat([{4},{'a'},{'scrambled'}],length(a_4_scrambled),1);

mycell(iso_3_sentence,:) = repmat([{3},{'iso'},{'sentence'}],length(iso_3_sentence),1);
mycell(iso_3_scrambled,:) = repmat([{3},{'iso'},{'scrambled'}],length(iso_3_scrambled),1);
mycell(a_3_sentence,:) = repmat([{3},{'a'},{'sentence'}],length(a_3_sentence),1);
mycell(a_3_scrambled,:) = repmat([{3},{'a'},{'scrambled'}],length(a_3_scrambled),1);

mycell(iso_5_sentence,:) = repmat([{5},{'iso'},{'sentence'}],length(iso_5_sentence),1);
mycell(iso_5_scrambled,:) = repmat([{5},{'iso'},{'scrambled'}],length(iso_5_scrambled),1);
mycell(a_5_sentence,:) = repmat([{5},{'a'},{'sentence'}],length(a_5_sentence),1);
mycell(a_5_scrambled,:) = repmat([{5},{'a'},{'scrambled'}],length(a_5_scrambled),1);

event_types = cell2table(mycell);
event_types.Properties.VariableNames = {'word_count','iso_a','sent_scr'};

end