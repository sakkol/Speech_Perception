function [classes] = IL_ranksum_SignClass(ranksumresults)

loop1={'iso','a'};
loop2={4,3,5};
loop3={'sentence','scrambled'};

if strcmp(comparison{2},'4')
    WPS_levels = {'24Hz','12Hz','06Hz'};
elseif strcmp(comparison{2},'3')
    WPS_levels = {'24Hz','08Hz'};
elseif strcmp(comparison{2},'5')
    WPS_levels = {'24Hz','048Hz'};
end

comps = fieldnames(ranksumresults);
classes={};
for c = 1:length(comps)
    if ranksumresults.(comps{c}).h
        classes{end+1,1} = comps{c};
    end
end

classes = replace(classes,'24Hz','W');
classes = replace(classes,'12Hz','P');
classes = replace(classes,{'06Hz','08Hz','048Hz'},'S');


end