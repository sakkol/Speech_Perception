function [classes] = IL_ranksum_class(ranksumresults,comparison)

to_comp = find(cellfun(@iscell,comparison));
if to_comp==1
    base = ['xxxx_' comparison{2} '_' comparison{3}];
elseif to_comp==2
    base = [comparison{1} '_xxxx_' comparison{3}];
elseif to_comp==3
    base = [comparison{1} '_' comparison{2} '_xxxx'];
else
    base = [comparison{1} '_' comparison{2} '_' comparison{3}];
end

if strcmp(comparison{2},'4')
    WPS_levels = {'24Hz','12Hz','06Hz'};
elseif strcmp(comparison{2},'3')
    WPS_levels = {'24Hz','08Hz'};
elseif strcmp(comparison{2},'5')
    WPS_levels = {'24Hz','048Hz'};
end

if isempty(to_comp)
    classes = cell(1,1);
else
    classes = cell(1,length(comparison{to_comp}));
end

for ll = 1:size(classes,2)
    if isempty(to_comp)
        main_look = '';
    else
        main_look = comparison{to_comp}{ll};
    end
    
    if length(WPS_levels)==3
        classes{1,ll} = [main_look '_',...
            repmat('W',1,ranksumresults.([replace(base,'xxxx',main_look) '_' WPS_levels{1}]).h),...
            repmat('P',1,ranksumresults.([replace(base,'xxxx',main_look) '_' WPS_levels{2}]).h),...
            repmat('S',1,ranksumresults.([replace(base,'xxxx',main_look) '_' WPS_levels{3}]).h)];
        
    else
        classes{1,ll} = [main_look '_',...
            repmat('W',1,ranksumresults.([replace(base,'xxxx',main_look) '_' WPS_levels{1}]).h),...
            repmat('S',1,ranksumresults.([replace(base,'xxxx',main_look) '_' WPS_levels{2}]).h)];
    end
    if startsWith(classes{1,ll},'_')
        classes{1,ll} = erase(classes{1,ll},'_');
    end
end
end