gusWin= gausswin(75)/sum(gausswin(75));
isofrequency = 2.4;
manySentence = 5;
word_per_sentence = {3,4,5};

for w = 1:3
    no_words = word_per_sentence{w}*manySentence;
    latency        = [-.1 (no_words/isofrequency)+.1];
    
    % preparing the features
    timeLine = latency(1):0.01:latency(2);
    word_onsets = ((1:1:no_words)-1)/isofrequency;
    words_timeline = zeros(size(timeLine));words_timeline(nearest_multiple(timeLine,word_onsets)) = 1;
    sentence_onsets = ((1:word_per_sentence{w}:no_words)-1)/isofrequency;
    sentence_timeline = zeros(size(timeLine));sentence_timeline(nearest_multiple(timeLine,sentence_onsets)) = 1;
    
    
    % run gaussian convolution
    words_timeline = convn(words_timeline,gusWin,'same');
    sentence_timeline = convn(sentence_timeline,gusWin,'same');
    
    % plot the results
    figure('Position',[0 0 1400 500])
    plot(timeLine,words_timeline,'b')
    hold on
    plot(timeLine,sentence_timeline - 0.03,'g')
    
    
    if word_per_sentence{w}==4
        phrase_onsets = ((1:2:no_words)-1)/isofrequency;
        phrase_timeline = zeros(size(timeLine));phrase_timeline(nearest_multiple(timeLine,phrase_onsets)) = 1;
        phrase_timeline = convn(phrase_timeline,gusWin,'same');
        plot(timeLine,phrase_timeline - 0.06,'r')
        legend({'Words','Sentences','Phrases'})
    else
        legend({'Words','Sentences'})
    end
    yticks([]);
    xlim(latency)
    title(['Feature structure for ' num2str(word_per_sentence{w}) '-word sentences'])
    print(fullfile('/media/sakkol/HDD1/HBML/PROJECTS_DATA/word2vec',['features_' num2str(word_per_sentence{w}) 'wordsentence.png']),'-dpng')
end
