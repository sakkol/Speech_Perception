function [lineprops,out] = IL_powerpeak_plot(sampstr,whichfield,plot_cl)
% Helper function to plot the average power, hfa or erp data

% Check inputs
if ~isfield(sampstr.(whichfield), {'powspctrm','trial'})
    error('Input ftrip is not recognized. Make sure it has powspctrm or trial field.')
elseif isfield(sampstr.(whichfield), 'powspctrm')
    to_conv = sampstr.(whichfield).powspctrm;
    if contains(whichfield,'hfa')
        navg_dim = 1;
        xaxis = sampstr.(whichfield).time;
    elseif contains(whichfield,'wlt')
        navg_dim = 4;
        xaxis = sampstr.(whichfield).freq;
    end
elseif isfield(sampstr.(whichfield), 'trial')
    to_conv = sampstr.(whichfield).trial;
    navg_dim = 1;
    xaxis = sampstr.(whichfield).time;
end

% average in desired dimension
out=squeeze(nanmean(to_conv(:,1,:,:),navg_dim));

if isnan(out)
    out=nan([1,length(xaxis)]);
elseif isvector(out)
    out=out';
end

% plot
forplot = nanmean(out,1);
hold on;
lineprops=plot(xaxis,forplot,plot_cl);

end