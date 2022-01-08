function [significant,tresults] = IL_power_peak_stat(for_avg,freq,freq_to_check)

to_comp1 = mean(for_avg(:,freq<freq_to_check+0.25 & freq>freq_to_check-0.25),2);
[~,y]=(min(abs(freq - freq_to_check)));
to_comp2 = for_avg(:,y);

[significant,p,ci,stats] = ttest(to_comp1,to_comp2,'Alpha',0.05/250);
tresults.h=significant;
tresults.p=p;
tresults.ci=ci;
tresults.stats=stats;

end