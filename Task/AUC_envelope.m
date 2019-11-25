cd /home/sakkol/Documents/BACKUPS/Speech_Perception_BACKUP(2019.09.11)/SpeechPerception_Task(backup-11.09.2019)/SpeechPerception_Task/main_stim_loc/Sentences_Rate0.9

files = dir('*.wav');
amp=[];
amp{length(files),3}='';
prominent_peaks=[];
prominent_peaks{length(files),1}='';
% pow_avg=[];

fprintf('Iteration = 00000');
for i = 1:length(files)
    % load sentences
    fprintf('\b\b\b\b\b%5d',i)
    filename = files(i).name;
     
    [y,Fs] = audioread(filename);
    % Remove trailing zeros
    y = y(1:find(y,1,'last'));

    %compute envelope
    % 1. band pass filtering between 100 Hz and 4000 Hz 
    filter_range = [100 4000];
    [b,a] = butter(2, filter_range/(Fs/2), 'bandpass');
    yfilt = filter(b,a,y);
    
	
    % 2. Compute Hilbert transform and take the absolute value to get the
    % amplitude of the envelope
    amp{i,1} = filename;
    amp{i,2} = smoothdata(abs(hilbert(yfilt)),'gaussian',250);
    amp{i,3} = trapz(smoothdata(abs(hilbert(yfilt)),'gaussian',250));
    amp{i,4} = length(amp{i,2});
    
%     % compute envelope's power spectrum
%     L=length(amp{i,2});
%     L = 48000;  % in order to make every pow in same length
%     NFFT = 2^nextpow2(L); % Next power of 2 from length of y
%     Y = fft(amp{i,2},NFFT)/L;
%     freqf = Fs/2*linspace(0,1,NFFT/2+1);
%     pow = abs(Y(1:NFFT/2+1)).^2;
%     
%     % find most prominent peaks
%     [pks,locs] = findpeaks(pow(freqf>1.7&freqf<20)/max(pow((freqf>1.7&freqf<20))),freqf(freqf>1.7&freqf<20),'SortStr','descend');
%     
%     prominent_peaks{i,1} = filename;
%     prominent_peaks{i,2} = locs(1);
%     
%     % Plot amplitude spectrum for each stimuli: optional
%     figure('Position', [50 50  1000 1000]);
%     hold on
%     subplot(2,1,1)
%     findpeaks(pow/max(pow((freqf>0.7))),freqf)
%     text(locs+.02,pks,num2str(locs',3));
%     xlim([1 8])
%     ylim([0,1.1])
%     title(['Peaks in current sentence; speech is filtered in ' num2str(filter_range) 'Hz'])
%     set(gca, 'fontname','arial');
%     sgtitle(sprintf('Envelope power of: %s',filename))
%     
%     subplot(2,1,2)
%     
%     plot(1/Fs:1/Fs:length(amp{i,2})/Fs,smoothdata(amp{i,2},'gaussian',250))
%     xlim([1/Fs length(amp{i,2})/Fs])
%     
%     
%     if ~exist('Envelopes','dir'),mkdir Envelopes;end
%     print('-r300','-djpeg',['Envelopes' filesep sprintf('Envelope_%s',erase(filename,'.wav'))])
%     close all
% 
%     pow_avg(i,:) = pow;
end
fprintf('\n')

%%
% This is the multiplication index that we need to multiply by the max we
% aim
((sum([amp{:,3}])/sum([amp{:,4}]))/max(vertcat(amp{:,2})))