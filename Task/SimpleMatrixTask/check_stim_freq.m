files = dir('*.wav');
amp=[];
amp{length(files),1}='';
pow_avg=[];
% pow_avg{length(files),1}=[];
for i = 1:length(files)
    % load sentences
 
    filename = files(i).name;
     
    [y,Fs] = audioread(filename);
     
     
    %compute envelope
 
    % 1. band pass filtering between 100 Hz and 4000 Hz 
    filter_range = [1 12000];
    [b,a] = butter(2, filter_range/(Fs/2), 'bandpass');
    yfilt = filter(b,a,y);
 
     
    % 2. Compute Hilbert transform and take the absolute value to get the
    % amplitude of the envelope
    amp{i,1} = filename;
    amp{i,2} = abs(hilbert(yfilt));     
     
    % compute envelope's power spectrum
%     L=length(amp{i,2});
    L = 53978;  % in order to make every pow in same length
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    Y = fft(amp{i,2},NFFT)/L;
    freqf = Fs/2*linspace(0,1,NFFT/2+1);
    pow = abs(Y(1:NFFT/2+1)).^2;
    
    % Plot amplitude spectrum for each stimuli
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
    [pks,locs] = findpeaks(pow(freqf>0.7&freqf<20)/max(pow((freqf>0.7&freqf<20))),freqf(freqf>0.7&freqf<20));
    findpeaks(pow/max(pow((freqf>0.7))),freqf)
    text(locs+.02,pks,num2str(locs',3));
    xlim([1 8])
    ylim([0,1.1])
    title(['Peaks in current sentence; speech is filtered in ' num2str(filter_range) 'Hz'])
    set(gca, 'fontname','arial');
    print('-r300','-djpeg',sprintf('Envelope_%s',erase(filename,'.wav')))
    close all
    
    pow_avg(i,:) = pow;
end


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
title('Average of envelope power of all stim Rate:1.7','fontsize',15,'fontname','arial')

subplot(2,1,2)
[pks,locs] = findpeaks(mean(pow_avg(:,freqf>0.7&freqf<20),1)/max(mean(pow_avg(:,freqf>0.7&freqf<20),1)),freqf(freqf>0.7&freqf<20));
findpeaks(mean(pow_avg(:,freqf>0.7&freqf<20),1)/max(mean(pow_avg(:,freqf>0.7&freqf<20),1)),freqf(freqf>0.7&freqf<20))
text(locs+.02,pks,num2str(locs',3));
xlim([1 8])
ylim([0,1.1])
title(['Peaks in current sentence; speech is filtered in ' num2str(filter_range) 'Hz'],'fontsize',15,'fontname','arial')

print('-r300','-djpeg','Mean_envelope_power_Rate1.7')
%save('envelopes_of_stims.mat','amp')