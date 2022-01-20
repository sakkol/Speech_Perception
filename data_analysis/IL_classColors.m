function [ClassColors] = IL_classColors(classes,comparison)

to_comp = find(cellfun(@iscell,comparison));
ClassColors = zeros(size(classes,1),3);
class_types = {'W', 'P', 'S', 'WP', 'WS', 'WPS'};
if isempty(to_comp) % just 1 dimensional
    [cmap] = master_ColorMaps('plasma',6);
    for c = 1:length(classes)
        ClassColors(c,:) = cmap(strcmp(classes{c},class_types),:);
    end
    
elseif length(comparison{to_comp})==2
    
    
    [cmap1] = master_ColorMaps('winter',6);
    [cmap2] = master_ColorMaps('spring',6);
    [cmap3] = master_ColorMaps('PuOr',6);
    base1 = strcat(repmat(comparison{to_comp}(1),1,6),repmat({'_'},1,6),class_types);
    base2 = strcat(repmat(comparison{to_comp}(2),1,6),repmat({'_'},1,6),class_types);
    
    
    for c = 1:size(classes,1)
        
        
        if any(strcmp(classes{c,1},base1)) && ~any(strcmp(classes{c,2},base2)) 
            ClassColors(c,:) = cmap1(strcmp(classes{c,1},base1),:);
            
        elseif ~any(strcmp(classes{c,1},base1)) && any(strcmp(classes{c,2},base2)) 
            ClassColors(c,:) = cmap2(strcmp(classes{c,2},base1),:);
        elseif any(strcmp(classes{c,1},base1)) && any(strcmp(classes{c,2},base2)) 
            ClassColors(c,:) = cmap3(endsWith(classes{c,2},strcat(repmat({'_'},1,6),class_types)),:);
        end
        
        
        
        
        
        
        
        
    end
    
elseif length(comparison{to_comp})==3
    error('Not done yet')
    
end


end