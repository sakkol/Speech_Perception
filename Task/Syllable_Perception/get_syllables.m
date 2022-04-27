cd /home/sakkol/Documents/TASKS/Syllable_Perception/Sounds_copy

wav60s = dir('*60.wav');

select_sylls = {'pa','ba','ta','da','fa','va','sa','za'};
select_voices = {'f106','f119','m114','m111'};

sel_wav60s=cell(0,0);
for i=1:length(wav60s)
    if any(contains(wav60s(i).name,select_sylls)) && any(contains(wav60s(i).name,select_voices))
        sel_wav60s(end+1,1:6) = table2cell(struct2table(wav60s(i)));
        splitted = strsplit(wav60s(i).name,'_');
        sel_wav60s{end,7} = splitted{3};
        
        [sel_wav60s{end,8},sel_wav60s{end,9}] = audioread(fullfile(sel_wav60s{end,2},sel_wav60s{end,1}));
        sel_wav60s{end,8} = resample(sel_wav60s{end,8},44100,sel_wav60s{end,9});
        sel_wav60s{end,9} = 44100;
    end
end
sel_wav60s(:,2:6)=[];
syll_table = cell2table(sel_wav60s);
syll_table.Properties.VariableNames={'filename','StimName','Stim','SampleRate'};

save('/home/sakkol/Documents/Codes_git/Speech_Perception/Task/Syllable_Perception/syll_table.mat','syll_table')
%% arrange trials



