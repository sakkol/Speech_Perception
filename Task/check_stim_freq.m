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
    [b,a] = butter(2, [100 4000]/(Fs/2), 'bandpass');
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
%     figure(1)
%     subplot(5,1,i)
%     hold on
%     plot(freqf,pow,'r') 
%     pow_avg{i} = pow;
    pow_avg(i,:) = pow;
end

figure(2)
hold on 
for i=1:length(files)
    lengths_all(i) = length(amp{i,2});
end
plot(freqf,mean(pow_avg)/max(mean(pow_avg(:,(freqf>0.7)))),'r','linewidth',4) 
set(gca, 'XTick',[1, 3, 5, 8],'fontsize',20,'fontname','arial'); 
set(gca,'YTick',[0,1],'fontsize',20,'fontname','arial'); 
xlabel('Frequency (Hz)','fontsize',22,'fontname','arial')
ylabel('Normalized power','fontsize',22,'fontname','arial')
xlim([1 8])
ylim([0,1.1])

%save('envelopes_of_stims.mat','amp')