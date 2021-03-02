%% plot the syllables
figure('Units','normalized','Position', [0 0  1 1]);
for s = 1:size(syll_table,1)
    subplot(4,8,s)
    plot((1:length(syll_table.Stim{s}))/44100,syll_table.Stim{s})
%     plot(syll_table.Stim{s})
    title([num2str(s) ' - ' syll_table.filename{s}],'Interpreter','none')
end

%% trim from beginning and the end for a few of them

syll_table.Stim{11} = syll_table.Stim{11}(9848:28383);
syll_table.Stim{12} = syll_table.Stim{12}(12400:24670);
syll_table.Stim{13} = syll_table.Stim{13}(16720:33580);
syll_table.Stim{14} = syll_table.Stim{14}(24280:end);
syll_table.Stim{15} = syll_table.Stim{15}(12490:29520);
syll_table.Stim{16} = syll_table.Stim{16}(14100:30950);

%% maximize the amplitude
for s = 1:size(syll_table,1)
    syll_table.Stim{s} = syll_table.Stim{s}/(max(abs(syll_table.Stim{s})));
end

%% play the sounds
for s = 1:size(syll_table,1)
    s
    playblocking(audioplayer(syll_table.Stim{s},44100))
    pause(0.1)
end

%% select only 2 speakers (1 F, 1 M)
% syll_table([9:16,25:32],:) = [];

figure('Units','normalized','Position', [0 0  1 1]);
for s = 1:size(syll_table,1)
    subplot(2,8,s)
    plot((1:length(syll_table.Stim{s}))/44100,syll_table.Stim{s})
%     plot(syll_table.Stim{s})
    title([num2str(s) ' - ' syll_table.filename{s}],'Interpreter','none')
end

%% save
save syll_table.mat syll_table