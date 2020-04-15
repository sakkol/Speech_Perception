% cd /home/sakkol/Documents/Speech_Perception_stim/03rd_Generation/GoogleTTS_M_Rate1

files = dir('*.wav');
amp=[];
amp{length(files),1}='';
prominent_peaks=[];
prominent_peaks{length(files),1}='';
pow_avg=[];

fprintf('Iteration = 00000');
for i = 1:length(files)
    % load sentences
    fprintf('\b\b\b\b\b%5d',i)
    filename = files(i).name;
     
    [y,Fs] = audioread(filename);
     
     
    %compute envelope
 
    % 1. band pass filtering between 100 Hz and 4000 Hz 
    filter_range = [100 4000];
    [b,a] = butter(2, filter_range/(Fs/2), 'bandpass');
    yfilt = filter(b,a,y);
 
     
    % 2. Compute Hilbert transform and take the absolute value to get the
    % amplitude of the envelope
    amp{i,1} = filename;
    amp{i,2} = abs(hilbert(yfilt));     
     
    % compute envelope's power spectrum
%     L=length(amp{i,2});
    L = 48000;  % in order to make every pow in same length
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    Y = fft(amp{i,2},NFFT)/L;
    freqf = Fs/2*linspace(0,1,NFFT/2+1);
    pow = abs(Y(1:NFFT/2+1)).^2;
    
    % find most prominent peaks
    [pks,locs] = findpeaks(pow(freqf>1.7&freqf<20)/max(pow((freqf>1.7&freqf<20))),freqf(freqf>1.7&freqf<20),'SortStr','descend');
    
    prominent_peaks{i,1} = filename;
    prominent_peaks{i,2} = locs(1);
    
    % Plot amplitude spectrum for each stimuli: optional
    figure('Position', [50 50  1000 1000]);
    hold on
    subplot(2,1,1)
    plot(freqf,pow/max(pow((freqf>0.7))),'b','linewidth',3)
    set(gca, 'XTick',[1, 3, 5, 8],'fontname','arial');
    set(gca,'YTick',[0,0.1,0.2,.5,1],'fontname','arial');
    xlabel('Frequency (Hz)','fontname','arial')
    ylabel('Normalized power','fontname','arial')
    xlim([1 8])
    ylim([0,1])
    title(sprintf('Envelope power of: %s',filename))
    
    subplot(2,1,2)
    
    findpeaks(pow/max(pow((freqf>0.7))),freqf)
    text(locs+.02,pks,num2str(locs',3));
    xlim([1 8])
    ylim([0,1.1])
    title(['Peaks in current sentence; speech is filtered in ' num2str(filter_range) 'Hz'])
    set(gca, 'fontname','arial');
    
    if ~exist('Envelopes','dir'),mkdir Envelopes;end
    print('-r300','-djpeg',['Envelopes' filesep sprintf('Envelope_%s',erase(filename,'.wav'))])
    close all
    
    pow_avg(i,:) = pow;
end
fprintf('\n')

save('prominent_peaks_each_stim.mat','prominent_peaks')

%% To get overall: average of envelope amplitudes
figure('Position', [50 50  900 1000]);
hold on 
subplot(2,1,1)
plot(freqf,mean(pow_avg)/max(mean(pow_avg(:,(freqf>0.7)))),'r','linewidth',3) 
set(gca, 'XTick',[1, 3, 5, 8],'fontsize',13,'fontname','arial'); 
set(gca,'YTick',[0,.5,1],'fontsize',13,'fontname','arial'); 
xlabel('Frequency (Hz)','fontsize',13,'fontname','arial')
ylabel('Normalized power','fontsize',13,'fontname','arial')
xlim([1 8])
ylim([0,1.1])
title('Average of envelope power of all stim','fontsize',15,'fontname','arial')

subplot(2,1,2)
[pks,locs] = findpeaks(mean(pow_avg(:,freqf>0.7&freqf<20),1)/max(mean(pow_avg(:,freqf>0.7&freqf<20),1)),freqf(freqf>0.7&freqf<20));
findpeaks(mean(pow_avg(:,freqf>0.7&freqf<20),1)/max(mean(pow_avg(:,freqf>0.7&freqf<20),1)),freqf(freqf>0.7&freqf<20))
text(locs+.02,pks,num2str(locs',3));
xlim([1 8])
ylim([0,1.1])
title(['Peaks in current sentence; speech is filtered in ' num2str(filter_range) 'Hz'],'fontsize',15,'fontname','arial')

if ~exist('Envelopes','dir'),mkdir Envelopes;end
print('-r300','-djpeg',['Envelopes' filesep 'Mean_envelope_power'])

%% Find sentences with prominent peaks in +/-0.3Hz range of attention sentence

p=0;
for i=1:length(prominent_peaks)
    if prominent_peaks{i,2}>3 && prominent_peaks{i,2}<3.6
        p=p+1;
        in_range_sentences(p,:) = prominent_peaks(i,:);
    end
end

p=0;
for i=1:length(prominent_peaks)
    if prominent_peaks{i,2}>4.1 && prominent_peaks{i,2}<4.7
        p=p+1;
        in_range_sentences(p,:) = prominent_peaks(i,:);
    end
end

save('in_range_sentences(0.3).mat','in_range_sentences')

inrange09 = load('in_range_sentences(0.3).mat')
inrange13 = load('in_range_sentences(0.3).mat')

%% Average of in_range_sentences

pow_avg_in_range=[];
p=0;
for i = 1:length(files)
    if prominent_peaks{i,2}>3 && prominent_peaks{i,2}<3.6
        p=p+1;
        pow_avg_in_range(p,:) = pow_avg(i,:);
    end
end
p=0;
for i=1:length(prominent_peaks)
    if prominent_peaks{i,2}>4.1 && prominent_peaks{i,2}<4.7
        p=p+1;
        pow_avg_in_range(p,:) = pow_avg(i,:);
    end
end

figure('Position', [50 50  900 1000]);
hold on 
subplot(2,1,1)
plot(freqf,mean(pow_avg_in_range)/max(mean(pow_avg_in_range(:,(freqf>0.7)))),'r','linewidth',3) 
set(gca, 'XTick',[1, 3, 5, 8],'fontsize',13,'fontname','arial'); 
set(gca,'YTick',[0,.5,1],'fontsize',13,'fontname','arial'); 
xlabel('Frequency (Hz)','fontsize',13,'fontname','arial')
ylabel('Normalized power','fontsize',13,'fontname','arial')
xlim([1 8])
ylim([0,1.1])
title('Average of envelope power of in range stim (+/-.3)','fontsize',15,'fontname','arial')

subplot(2,1,2)
[pks,locs] = findpeaks(mean(pow_avg_in_range(:,freqf>0.7&freqf<20),1)/max(mean(pow_avg_in_range(:,freqf>0.7&freqf<20),1)),freqf(freqf>0.7&freqf<20));
findpeaks(mean(pow_avg_in_range(:,freqf>0.7&freqf<20),1)/max(mean(pow_avg_in_range(:,freqf>0.7&freqf<20),1)),freqf(freqf>0.7&freqf<20))
text(locs+.02,pks,num2str(locs',3));
xlim([1 8])
ylim([0,1.1])
title(['Peaks in average; speech is filtered in ' num2str(filter_range) 'Hz'],'fontsize',15,'fontname','arial')

if ~exist('Envelopes','dir'),mkdir Envelopes;end
print('-r300','-djpeg',['Envelopes' filesep 'Mean_envelope_power_in_range'])


%% After grouping
% find common sentences from each speech rates
for i=1:length(inrange09.in_range_sentences)
    inrange09.in_range_sentences{i,3}=ismember(inrange09.in_range_sentences{i,1},inrange13.in_range_sentences(:,1));
end

n=1;
for i=1:length(inrange09.in_range_sentences)
    if inrange09.in_range_sentences{i,3} == 1
        copyfile(inrange09.in_range_sentences{i,1},'/home/sakkol/Documents/Spanish_Matrix_Sentence/Preparations/Version_4/Sentences_Rate0.9')
        fprintf('copied%d\n',i)
%         temp = erase(inrange09.in_range_sentences{i,1},'-M.wav');
%         sentence_to_use_list{n,1} = replace(temp,'_',' ');
        sentence_to_use_list{n,1} = inrange09.in_range_sentences{i,1};
        n=n+1;
    end
end

%% Seelcting for adaptation and thresholding

commons = readtable('/home/sakkol/Documents/Spanish_Matrix_Sentence/Version_2/common_sentences.xlsx');
comm_noThresh = readtable('/home/sakkol/Documents/Spanish_Matrix_Sentence/Version_2/common_sentences_noThresh.xlsx');

% Thresholding: 30 sentences, need to be as much different as possible
ufff=1;not_gogogo=1;fprintf('Iteration = 000000');
while not_gogogo && ufff <900000
    fprintf('\b\b\b\b\b\b%6d',ufff)
    chosen_rand = randperm(711,30);
    
    [ii,jj,kk]=unique(comm_noThresh.Name(chosen_rand,:));
    Name_f=histc(kk,1:numel(jj)); % Frequencies corresponding to ii
    [ii,jj,kk]=unique(comm_noThresh.Verb(chosen_rand,:));
    Verb_f=histc(kk,1:numel(jj)); % Frequencies corresponding to ii
    [ii,jj,kk]=unique(comm_noThresh.Adjective(chosen_rand,:));
    Adjective_f=histc(kk,1:numel(jj)); % Frequencies corresponding to ii
    [ii,jj,kk]=unique(comm_noThresh.Numeral(chosen_rand,:));
    Number_f=histc(kk,1:numel(jj)); % Frequencies corresponding to ii
    [ii,jj,kk]=unique(comm_noThresh.Object(chosen_rand,:));
    Noun_f=histc(kk,1:numel(jj)); % Frequencies corresponding to ii
    
    
    if sum(Name_f==3) == 7 && ...
            sum(Noun_f==3) == 7 && ...
            sum(Verb_f==3) == 7 && ...
            sum(Adjective_f==3) == 7 && ...
            sum(Number_f==3) == 7 %&& ...
%             size(unique(commons(chosen_rand,3:4)),1) > 15 && ...
%             size(unique(commons(chosen_rand,5:6)),1) > 15 && ...
%             size(unique(commons(chosen_rand,6:7)),1) > 15
            
            not_gogogo = 0;
        threshold_sentences = comm_noThresh(chosen_rand,:)
    end
    ufff=ufff+1;
end
fprintf('\n')

% 
word_list = readtable('/home/sakkol/Documents/Spanish_Matrix_Sentence/spanish_matrix_words.xlsx');
ufff=1;not_gogogo=1;fprintf('Iteration = 0000000');
% selc_sents={''};p=1;
while not_gogogo && ufff <9000000
    fprintf('\b\b\b\b\b\b\b%7d',ufff)
%     rand5 = randperm(10,5);
    c_name = word_list.Name(randi(10,1));
    c_verb = word_list.Verb(randi(10,1));
    c_numeral = word_list.Numeral(randi(10,1));
    c_obj = word_list.Object(randi(10,1));
    c_adj = word_list.Adjective(randi(10,1));
    c_sentence = strjoin([c_name(:)',c_verb(:)',c_numeral(:)',c_obj(:)',c_adj(:)'],' ');
    if ismember(c_sentence,comm_noThresh.Sentence) && ...
            ~ismember(c_sentence,selc_sents) && ...
            sum(contains(selc_sents,c_name))<4 && ...
            sum(contains(selc_sents,c_verb))<4 && ...
            sum(contains(selc_sents,c_numeral))<4 && ...
            sum(contains(selc_sents,c_obj))<4 && ...
            sum(contains(selc_sents,c_adj))<4
        
        selc_sents{p} = c_sentence;
        p=p+1;
        if p>30
            not_gogogo=0;
            selc_sents
        end
        
    end
    
    ufff=ufff+1;
end
fprintf('\n')

% I'm angry
word_list = readtable('/home/sakkol/Documents/Spanish_Matrix_Sentence/spanish_matrix_words.xlsx');
fprintf('Iteration = 0000000');
selc_sents_angry={''};p=1;
ufff=1;
for i1=1:10
    for i2=1:10
        for i3=1:10
            for i4=1:10
                for i5=1:10
                    fprintf('\b\b\b\b\b\b\b%7d',ufff)
                    ufff=ufff+1;
                    c_name = word_list.Name(i1);
                    c_verb = word_list.Verb(i2);
                    c_numeral = word_list.Numeral(i3);
                    c_obj = word_list.Object(i4);
                    c_adj = word_list.Adjective(i5);
                    c_sentence = strjoin([c_name(:)',c_verb(:)',c_numeral(:)',c_obj(:)',c_adj(:)'],' ');
                    if ismember(c_sentence,comm_noThresh.Sentence) && ...
                            ~ismember(c_sentence,selc_sents_angry) && ...
                            sum(contains(selc_sents_angry,c_name))<4 && ...
                            sum(contains(selc_sents_angry,c_verb))<4 && ...
                            sum(contains(selc_sents_angry,c_numeral))<3 && ...
                            sum(contains(selc_sents_angry,c_obj))<4 && ...
                            sum(contains(selc_sents_angry,c_adj))<3
                        
                        selc_sents_angry{p} = c_sentence;
                        p=p+1;
                        if p>30
                            
                            selc_sents
                            break;
                        end
                    end
                end
            end
        end
    end
    
end
fprintf('\n')

% Adaptation: 30 sentences, no need to be very different from each other
ufff=1;not_gogogo=1;fprintf('Iteration = 000000');
while not_gogogo && ufff <10000
    fprintf('\b\b\b\b\b\b%6d',ufff)
    chosen_rand = randperm(713,30);
    if length(unique(commons.Name(chosen_rand,:))) == 10 && ...
            length(unique(commons.Verb(chosen_rand,:))) == 10 && ...
            length(unique(commons.Number(chosen_rand,:))) == 10 && ...
            length(unique(commons.Adjective(chosen_rand,:))) == 10 && ...
            length(unique(commons.Noun(chosen_rand,:))) == 10 && ...
            ~any(ismember(commons.Sentence(chosen_rand,:),threshold_sentences.Sentence))
        not_gogogo = 0;
        adaptation_sentences = commons(chosen_rand,:)
    end
    ufff=ufff+1;
end
fprintf('\n')

% remove these adaptation and threshold from main list
n=1;
for i=1:size(commons,1)
    if ~ismember(commons.Sentence{i},threshold_sentences.Sentence) && ~ismember(commons.Sentence{i},adaptation_sentences.Sentence)
        sentence_to_use(n,:) = commons(i,:);n=n+1;
    end
end

% copy 