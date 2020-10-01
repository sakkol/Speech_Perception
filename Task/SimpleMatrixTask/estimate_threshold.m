function [current_SNR] = estimate_threshold(slowVSfast,threshold_accr,threshLevel,save_filename)
% Threshold estimation to be used in master script or event_creator to
% decrease clutter over there.
%% Set default variables
vars=who;
if ~any(ismember(vars,'slowVSfast')) || isempty(slowVSfast)
    error('slowVSfast is a required input!')
end

if ~any(ismember(vars,'threshold_accr')) || isempty(threshold_accr)
    % Input accuracies for each SNR
    if slowVSfast == 1
        threshold_accr = inputdlg({'SNR = -8:','SNR = -4:','SNR = 0:'},'Please input how many words were correct for each SNR!',[1 70]);
        data_threshold =    [...
            -8,   str2double(threshold_accr{1}),   50.0000;...
            -4,   str2double(threshold_accr{2}),   50.0000;...
             0,   str2double(threshold_accr{3}),   50.0000];
    elseif slowVSfast == 2
        threshold_accr = inputdlg({'SNR = -6:','SNR = -2:','SNR = 2:'},'Please input how many words were correct for each SNR!',[1 70]);
        data_threshold =    [...
            -6,   str2double(threshold_accr{1}),   50.0000;...
            -2,   str2double(threshold_accr{2}),   50.0000;...
             2,   str2double(threshold_accr{3}),   50.0000];
    elseif slowVSfast == 3
        threshold_accr = inputdlg({'SNR = -2:','SNR = 2:','SNR = 6:'},'Please input how many words were correct for each SNR!',[1 70]);
        data_threshold =    [...
            -2,   str2double(threshold_accr{1}),   40.0000;...
             2,   str2double(threshold_accr{2}),   40.0000;...
             6,   str2double(threshold_accr{3}),   40.0000];
    end
elseif length(threshold_accr) ~= 3
    error('threshold_accr should be 3 variables in length!')
else % if only 3 values were given
    if slowVSfast == 1
        data_threshold =    [...
            -8,   threshold_accr(1),   50.0000;...
            -4,   threshold_accr(2),   50.0000;...
             0,   threshold_accr(3),   50.0000];
    elseif slowVSfast == 2
        data_threshold =    [...
            -6,   threshold_accr(1),   50.0000;...
            -2,   threshold_accr(2),   50.0000;...
             2,   threshold_accr(3),   50.0000];
    elseif slowVSfast == 3
        data_threshold =    [...
            -2,   threshold_accr(1),   40.0000;...
             2,   threshold_accr(2),   40.0000;...
             6,   threshold_accr(3),   40.0000];
    end
end

if ~any(ismember(vars,'threshLevel')) || isempty(threshLevel)
    threshLevel = .5;
end

if ~any(ismember(vars,'save_filename')) || isempty(save_filename)
    try
        save_filename = evalin('base','save_filename');
    catch
        save_filename = fullfile(pwd,'psignifit_res');
    end
end

%% Set options and run psignifit
options             = [];   % initialize as an empty struct
options.sigmoidName = 'norm';   % normal Gaussian curve as the sigmoid
% options.expType     = '4AFC'; % ours is not forced choice, so  no input should be given
options.threshPC       = threshLevel;

result = psignifit(data_threshold,options);
p1=plotPsych(result);
print([save_filename '_psignifit.jpg'],'-djpeg','-r300')
current_SNR = getThreshold(result,0.5)
% close p1

end