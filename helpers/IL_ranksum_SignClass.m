function [classes] = IL_ranksum_SignClass(ranksumresults)

% loop1={'iso','a'};
% loop2={4,3,5};
% loop3={'sentence','scrambled'};
% 
% if strcmp(comparison{2},'4')
%     WPS_levels = {'24Hz','12Hz','06Hz'};
% elseif strcmp(comparison{2},'3')
%     WPS_levels = {'24Hz','08Hz'};
% elseif strcmp(comparison{2},'5')
%     WPS_levels = {'24Hz','048Hz'};
% end

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

for c = 1:size(classes,1)
    if any(contains(classes,classes{c}(1:end-1))) && any(find(contains(classes,classes{c}(1:end-1)))~=c) && ~isempty(classes{c})
        to_add=find(contains(classes,classes{c}(1:end-1)));
        to_add(to_add==c)=[];
        for a = 1:length(to_add)
            classes{c} = [classes{c},classes{to_add(a)}(end)];
            classes{to_add(a)}='';
        end
            
        
    end
end
classes(cellfun(@isempty,classes))=[];

end