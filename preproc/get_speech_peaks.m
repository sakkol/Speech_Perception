function [peakRate, peakEnv, amp_env, deriv_amp_env] = get_speech_peaks(speech,Fs,ifPlot)
% This function calculates broadband speech amplitude envelope (amp_env),
% half-wave rectified first derivative of amplitude envelope
% (deriv_amp_env), peaks in the envelope (peakEnv), peaks in the rate
% changes of envelope (peakRate) and also plots speech spectrogram along
% with envelope and its first derivative.
% Spectrogram calculation is based on Cochlear Output function as
% implemented in NSL toolbox with critical filters (CF) = 440*2.^((-31:97)/24)).
% peakRate and peakEnv calculations are based on Oganian and Chang, 2019
% (Science Advances).
%
% Serdar Akkol, Human Brain Mapping Lab, Feinstein Institutes for Medical
% Research
% June 2020

%% Calculations
if nargin == 2, ifPlot=1;end

% remove zeros in the end
speech = speech(1:find(speech,1,'last'));

% get amplitude envelope and first derivative of amplitude envelope (using cochlear function scale)
[amp_env,deriv_amp_env,spectrogram] = get_speech_2features(speech,Fs,'cochlear');

% Remove frequency bands that exceeds half of sampling rate
CF = 440 * 2 .^ ((-31:97)/24);
CF(CF>Fs/2) = [];

% half wave rectification of derivative
deriv_amp_env(deriv_amp_env<0) = 0;

% find the peaks in both signals to get peakEnv and peak Rate locations
[~,peakEnv]=findpeaks(amp_env*(1/max(abs(speech))),'MinPeakDistance',400,'MinPeakProminence',0.0001);
[~,peakRate]=findpeaks(deriv_amp_env*(1/max(abs(deriv_amp_env))),'MinPeakDistance',400,'MinPeakProminence',0.02);

% remove if there is one very close to end, probably error
peakEnv( peakEnv > length(deriv_amp_env)-100) = [];
peakRate( peakRate > length(deriv_amp_env)-100) = [];

%% Plot results

if ifPlot
    figure('Units','normalized','Position', [0 0  .8 .5]);
    subplot(211)
    imagesc([1:length(spectrogram)]/Fs,CF,spectrogram)
    axis xy
    ylabel('Frequency (Hz)')
    set(gca, 'FontSize',13,'FontWeight','bold','XGrid','on','GridColor','w','GridLineStyle','--','GridAlpha',1,'LineWidth',.7);
    
    subplot(212)
    p1=plot([1:length(speech)]/Fs,speech*(1/max(abs(speech))),'Color',[0.651,0.651,0.651]);
    ylim=[-1 1];
    y=ylim;
    hold on
    p2=plot([1:length(amp_env)]/Fs,amp_env*(1/max(amp_env)),'r');
    if length(peakEnv) == 2
        plot([peakEnv(1)/Fs peakEnv(1)/Fs],[y(1) y(2)],'r','LineWidth',2);
        plot([peakEnv(2)/Fs peakEnv(2)/Fs],[y(1) y(2)],'r','LineWidth',2);
    elseif ~isempty(peakEnv)
        plot([peakEnv/Fs peakEnv/Fs],[y(1) y(2)],'r','LineWidth',2);
    end
    
    p4=plot([1:length(deriv_amp_env)]/Fs,deriv_amp_env*(1/max(abs(deriv_amp_env))),'k');
    if length(peakRate) == 2
        plot([peakRate(1)/Fs peakRate(1)/Fs],[y(1) y(2)],'k','LineWidth',2);
        plot([peakRate(2)/Fs peakRate(2)/Fs],[y(1) y(2)],'k','LineWidth',2);
    elseif ~isempty(peakRate)
        plot([peakRate/Fs peakRate/Fs],[y(1) y(2)],'k','LineWidth',2);
    end
    xlim([0 length(deriv_amp_env)/Fs])
    xlabel('Time (s)')
    ylabel('Normalized amplitude')
    legend([p1 p2 p4],{'Speech', 'Amplitude Envelope', 'Positive Envelope Rate of Change'},'Location','northoutside','Orientation','horizontal','FontSize',15);
    
    set(gca, 'FontSize',13,'FontWeight','bold','XGrid','on','GridColor','k','GridLineStyle','--','GridAlpha',1,'LineWidth',.7);
end

end

