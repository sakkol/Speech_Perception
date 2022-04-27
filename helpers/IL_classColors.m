function [ClassColors,classes_and] = IL_classColors(classes,comparison)
% Assigning colors to different 

% create all the combinations
% first cell all the parts
non_to_comp = find(~cellfun(@iscell,comparison));cellcomparison=comparison;
for n=1:length(non_to_comp)
    cellcomparison{non_to_comp(n)} = comparison(non_to_comp(n));
end
% now columns to rows
to_comp = find(cellfun(@iscell,comparison));
for n=1:length(to_comp)
    cellcomparison{to_comp(n)} = cellcomparison{to_comp(n)}';
end
% create combinations
ncomb = prod(cellfun(@length,comparison(to_comp)));
[Ax,Bx,Cx,Dx] = ndgrid(1:numel(cellcomparison{1}),1:numel(cellcomparison{2}),1:numel(cellcomparison{3}),1:numel(cellcomparison{4}));
all_comb = strcat(cellcomparison{1}(Ax(:)),repmat({'_'},ncomb,1),cellcomparison{2}(Bx(:)),...
    repmat({'_'},ncomb,1),cellcomparison{3}(Cx(:)),repmat({'_'},ncomb,1),cellcomparison{4}(Dx(:)));

% find the significant classes of that are compared
for c = 1:size(classes,1)
    if strcmp(classes{c},'Non-resp/Other')
        continue
    elseif any(ismember(classes{c},all_comb))
        classes{c} = classes{c}(ismember(classes{c},all_comb));
    else
        classes{c} = {'Non-resp/Other'};
    end
end

% check which classes are found
unique_classes = unique(cellfun(@(x)join(x,{'&'}),classes));%unique(cellfun(@(x)join(x,{'&'}),classes(~cellfun(@(x)strcmp(x,'Non-resp/Other'),classes)),'UniformOutput',0));
% for u = 1:length(unique_classes)
%     sep_class = strsplit(unique_classes{u},'&');
%     if ~any(ismember(sep_class,all_comb))
%         unique_classes{u}='';
%     else
%         unique_classes{u} = strjoin({sep_class{ismember(sep_class,all_comb)}},'&');
%     end
% end
% unique_classes = unique(unique_classes(~cellfun(@isempty,unique_classes)));

% create the color palette
cmap=distinguishable_colors(length(unique_classes),{'w','k'});
cmap=sortrows(cmap,[1,2,3]); % easily distinguishable colors, sorted in similarity
% if ncomb <= 6
%     [cmap] = master_ColorMaps('plasma',6);
% elseif ncomb > 6 && ncomb <= 18
%     [cmap] = [master_ColorMaps('winter',6); master_ColorMaps('spring',6);master_ColorMaps('PuOr',6)];
%     cmap = cmap(1:ncomb,:);
% else
%     error('Need more colors')
% end

classes_and = cellfun(@(x)join(x,{'&'}),classes,'UniformOutput',0);

ClassColors = zeros(size(classes,1),3);
for c = 1:size(classes,1)
    if isempty(classes_and{c})
        continue
    else
        ClassColors(c,:) = cmap(strcmp(classes_and{c},unique_classes),:);
    end
end


end