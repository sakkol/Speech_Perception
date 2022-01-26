function [significant,ranksumresults] = IL_powerpeak_stat(for_avg,freq,lineinfo,freq_to_check)

% find closest 2 frequencies
[~,y] = mink(abs(freq-freq_to_check),1);
to_comp1 = mean(for_avg(:,y),2);

% average the 2 frequency bins above and beyond excluding the freq of interest +/-2bin
to_comp2 = mean(for_avg(:,[y-3:y-2,y+2:y+3]),2);%to_comp1 = mean(for_avg(:,inds),2);

if ~isnan(to_comp1)
    [p,significant,stats] = ranksum(to_comp1,to_comp2,'alpha',0.05);
    ranksumresults.h=significant;
    ranksumresults.p=p;
    ranksumresults.stats=stats;
else
    significant=0;p=1;stats=[];
    ranksumresults.h=significant;
    ranksumresults.p=p;
    ranksumresults.stats=stats;
end

if significant
    plot(freq_to_check,lineinfo.YData(nearest(freq,freq_to_check)),'*','MarkerSize',10,'MarkerEdgeColor',lineinfo.Color,'LineWidth',1)
end

end