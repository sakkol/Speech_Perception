

%% SECOND HALF OF FIRST BLOCK
%% SECOND HALF OF FIRST BLOCK::

%% Behavioral analysis: accuracy of each block (raw trial numbers)
load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']))
syllable_couples = {'pa','ba';'ta','da';'fa','va';'sa','za'};
unvoiced_syll = {'pa';'ta';'fa';'sa'};
voiced_syll = {'ba';'da';'va';'za'};

info.events.Accr = num2cell(strcmp(info.events.SyllablePresented,info.events.Response));
info.events.Accr(strcmp(info.events.Response,'')) = {'noresp'};

consonant = cell(size(info.events,1),1);
for t=1:size(info.events,1)
    if ismember(info.events.SyllablePresented{t},voiced_syll)
        consonant{t} = 'voiced';
    else
        consonant{t} = 'unvoiced';
    end
    if strcmp(info.events.SyllablePresented{t},info.events.Response{t})
        info.events.Accr{t} = num2str(strcmp(info.events.SyllablePresented{t},info.events.Response{t}));
    elseif strcmp(info.events.Response{t},'')
        info.events.Accr{t} = 'noresp';
    else
        info.events.Accr{t} = '0';
    end
end
info.events = [info.events,cell2table(consonant)];

if strcmp(Sbj_Metadata.sbj_ID,'NS163') && strcmp(curr_block,'SyllablePerception_EStim')
    info.events = info.events(~ismember(1:64,[1,11,14,16,21,22,25,31,32,33,37,41,43,51,54,58,61,62,63,64]),:);
end
if strcmp(Sbj_Metadata.sbj_ID,'NS163') && strcmp(curr_block,'SyllablePerception')
    info.events = info.events(~ismember(1:92,1:46),:);
end

voiced_corr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'1'));
voiced_incorr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'0'));
voiced_noresp = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'noresp'));

unvoiced_corr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'1'));
unvoiced_incorr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'0'));
unvoiced_noresp = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'noresp'));

figure('Position',[100 100 1000 1000])
to_bar = [voiced_corr voiced_incorr voiced_noresp;...
          unvoiced_corr unvoiced_incorr unvoiced_noresp];
hBar = bar(to_bar,'stacked');
set(gca,'xticklabels',{'Voiced','Unvoiced'},'FontSize',15)
grid on
legend('Correct response','Incorrect Response','No-response')                                                  % Return ‘bar’ Handle
for k1 = 1:size(to_bar,2)
    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, hBar(k1).XOffset');    % Note: ‘XOffset’ Is An Undocumented Feature, This Selects The ‘bar’ Centres
    ydt(k1,:) = hBar(k1).YData;                                     % Individual Bar Heights
    if k1==1
        text(ctr(k1,:), ydt(k1,:), sprintfc('%.0f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',18, 'Color','k')
    else
        text(ctr(k1,:), ydt(k1,:)+sum(ydt(1:k1-1,:),1), sprintfc('%.0f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',18, 'Color','k')
    end
end

% Observed data
n1 = voiced_corr; N1 = voiced_incorr+voiced_noresp+voiced_corr;
n2 = unvoiced_corr; N2 = unvoiced_incorr+unvoiced_noresp+unvoiced_corr;
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval] = crosstab(x1,x2)
title(['Comparing accuracy between voiced vs unvoiced, p-val: ' num2str(pval)])

print(fullfile(Sbj_Metadata.results,['2ndhalf_' curr_block '_accuracy_onlyclean.jpg']),'-djpeg','-r300')

%% Behavioral analysis: accuracy of each block (percentage)
load(fullfile(Sbj_Metadata.iEEG_data,curr_block,[curr_block '_info.mat']))
syllable_couples = {'pa','ba';'ta','da';'fa','va';'sa','za'};
unvoiced_syll = {'pa';'ta';'fa';'sa'};
voiced_syll = {'ba';'da';'va';'za'};

info.events.Accr = num2cell(strcmp(info.events.SyllablePresented,info.events.Response));
info.events.Accr(strcmp(info.events.Response,'')) = {'noresp'};

consonant = cell(size(info.events,1),1);
for t=1:size(info.events,1)
    if ismember(info.events.SyllablePresented{t},voiced_syll)
        consonant{t} = 'voiced';
    else
        consonant{t} = 'unvoiced';
    end
    if strcmp(info.events.SyllablePresented{t},info.events.Response{t})
        info.events.Accr{t} = num2str(strcmp(info.events.SyllablePresented{t},info.events.Response{t}));
    elseif strcmp(info.events.Response{t},'')
        info.events.Accr{t} = 'noresp';
    else
        info.events.Accr{t} = '0';
    end
end
info.events = [info.events,cell2table(consonant)];

if strcmp(Sbj_Metadata.sbj_ID,'NS163') && strcmp(curr_block,'SyllablePerception_EStim')
    info.events = info.events(~ismember(1:64,[1,11,14,16,21,22,25,31,32,33,37,41,43,51,54,58,61,62,63,64]),:);
end
if strcmp(Sbj_Metadata.sbj_ID,'NS163') && strcmp(curr_block,'SyllablePerception')
    info.events = info.events(~ismember(1:92,1:46),:);
end

voiced_corr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'1'));
voiced_incorr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'0'));
voiced_noresp = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'noresp'));
voiced_total = voiced_corr +voiced_incorr+voiced_noresp;
unvoiced_corr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'1'));
unvoiced_incorr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'0'));
unvoiced_noresp = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'noresp'));
unvoiced_total = unvoiced_corr + unvoiced_incorr+unvoiced_noresp;

figure('Position',[100 100 1000 1000])
to_bar = 100*[voiced_corr/voiced_total voiced_incorr/voiced_total voiced_noresp/voiced_total;...
          unvoiced_corr/unvoiced_total unvoiced_incorr/unvoiced_total unvoiced_noresp/unvoiced_total];
hBar = bar(to_bar,'stacked');
set(gca,'xticklabels',{'Voiced','Unvoiced'},'FontSize',15)
grid on
legend('Correct response','Incorrect Response','No-response','Location','southoutside','NumColumns',3)
ctr=[];ydt=[];
for k1 = 1:size(to_bar,2)
    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, hBar(k1).XOffset');    % Note: ‘XOffset’ Is An Undocumented Feature, This Selects The ‘bar’ Centres
    ydt(k1,:) = hBar(k1).YData;                                     % Individual Bar Heights
    if k1==1
        text(ctr(k1,:), ydt(k1,:), sprintfc('%.0f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',18, 'Color','k')
    else
        text(ctr(k1,:), ydt(k1,:)+sum(ydt(1:k1-1,:),1), sprintfc('%.0f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',18, 'Color','k')
    end
end

% Observed data
n1 = voiced_corr; N1 = voiced_incorr+voiced_noresp+voiced_corr;
n2 = unvoiced_corr; N2 = unvoiced_incorr+unvoiced_noresp+unvoiced_corr;
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval] = crosstab(x1,x2)
title(['Comparing accuracy between voiced vs unvoiced, p-val: ' num2str(pval)])
ylim([0 100])

print(fullfile(Sbj_Metadata.results,['2ndhalf_' curr_block '_accuracy_onlyclean_percentage.jpg']),'-djpeg','-r300')

%% Comparing between blocks raw trial numbers
syllable_couples = {'pa','ba';'ta','da';'fa','va';'sa','za'};
unvoiced_syll = {'pa';'ta';'fa';'sa'};
voiced_syll = {'ba';'da';'va';'za'};
block=cell(1,2);
for i=1:2
    load(fullfile(Sbj_Metadata.iEEG_data,Sbj_Metadata.BlockLists{i},[Sbj_Metadata.BlockLists{i} '_info.mat']));
    info.events.Accr = num2cell(strcmp(info.events.SyllablePresented,info.events.Response));
    info.events.Accr(strcmp(info.events.Response,'')) = {'noresp'};
    
    consonant = cell(size(info.events,1),1);
    for t=1:size(info.events,1)
        if ismember(info.events.SyllablePresented{t},voiced_syll)
            consonant{t} = 'voiced';
        else
            consonant{t} = 'unvoiced';
        end
        if strcmp(info.events.SyllablePresented{t},info.events.Response{t})
            info.events.Accr{t} = num2str(strcmp(info.events.SyllablePresented{t},info.events.Response{t}));
        elseif strcmp(info.events.Response{t},'')
            info.events.Accr{t} = 'noresp';
        else
            info.events.Accr{t} = '0';
        end
    end
    info.events = [info.events,cell2table(consonant)];
    if strcmp(Sbj_Metadata.sbj_ID,'NS163') && strcmp(Sbj_Metadata.BlockLists{i},'SyllablePerception_EStim')
        info.events = info.events(~ismember(1:64,[1,11,14,16,21,22,25,31,32,33,37,41,43,51,54,58,61,62,63,64]),:);
    end
    if strcmp(Sbj_Metadata.sbj_ID,'NS163') && strcmp(Sbj_Metadata.BlockLists{i},'SyllablePerception')
        info.events = info.events(~ismember(1:92,1:46),:);
    end
    
    voiced_corr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'1'));
    voiced_incorr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'0'));
    voiced_noresp = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'noresp'));
    
    unvoiced_corr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'1'));
    unvoiced_incorr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'0'));
    unvoiced_noresp = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'noresp'));
    
    block{i}.voiced_corr = voiced_corr;
    block{i}.voiced_incorr = voiced_incorr;
    block{i}.voiced_noresp = voiced_noresp;
    block{i}.unvoiced_corr = unvoiced_corr;
    block{i}.unvoiced_incorr = unvoiced_incorr;
    block{i}.unvoiced_noresp = unvoiced_noresp;
end


figure('Position',[100 100 1200 1000])
to_bar = [block{1}.voiced_corr block{1}.voiced_incorr block{1}.voiced_noresp;...
          block{2}.voiced_corr block{2}.voiced_incorr block{2}.voiced_noresp;...
          block{1}.unvoiced_corr block{1}.unvoiced_incorr block{1}.unvoiced_noresp;...
          block{2}.unvoiced_corr block{2}.unvoiced_incorr block{2}.unvoiced_noresp];
hBar = bar(to_bar,'stacked');
set(gca,'xticklabels',{'Voiced-NoEStim','Voiced-EStim','Unvoiced-NoEStim','Unvoiced-EStim'},'FontSize',15)
grid on
legend('Correct response','Incorrect Response','No-response')                                                  % Return ‘bar’ Handle
ctr=[];ydt=[];
for k1 = 1:size(to_bar,2)
    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, hBar(k1).XOffset');    % Note: ‘XOffset’ Is An Undocumented Feature, This Selects The ‘bar’ Centres
    ydt(k1,:) = hBar(k1).YData;                                     % Individual Bar Heights
    if k1==1
        text(ctr(k1,:), ydt(k1,:), sprintfc('%.0f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',18, 'Color','k')
    else
        text(ctr(k1,:), ydt(k1,:)+sum(ydt(1:k1-1,:),1), sprintfc('%.0f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',18, 'Color','k')
    end
end

% Observed data
n1 = block{1}.voiced_corr; N1 = block{1}.voiced_corr+block{1}.voiced_incorr+block{1}.voiced_noresp;
n2 = block{2}.voiced_corr; N2 = block{2}.voiced_corr+block{1}.voiced_incorr+block{1}.voiced_noresp;
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl1,chi2stat1,pval1] = crosstab(x1,x2)

n1 = block{1}.unvoiced_corr; N1 = block{1}.unvoiced_corr+block{1}.unvoiced_incorr+block{1}.unvoiced_noresp;
n2 = block{2}.unvoiced_corr; N2 = block{2}.unvoiced_corr+block{1}.unvoiced_incorr+block{1}.unvoiced_noresp;
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl2,chi2stat2,pval2] = crosstab(x1,x2)
title({['Comparing accuracy of voiced between no-EStim vs EStim, p-val: ' num2str(pval1)];...
       ['Comparing accuracy of unvoiced between no-EStim vs EStim, p-val: ' num2str(pval2)]})

print(fullfile(Sbj_Metadata.results,['2ndhalf_' curr_block '_block_comp_onlyclean.jpg']),'-djpeg','-r300')
% print(fullfile(Sbj_Metadata.results,[curr_block '_block_comp.jpg']),'-djpeg','-r300')

%% Comparing between blocks with percentages, NOT raw trial numbers
syllable_couples = {'pa','ba';'ta','da';'fa','va';'sa','za'};
unvoiced_syll = {'pa';'ta';'fa';'sa'};
voiced_syll = {'ba';'da';'va';'za'};
block=cell(1,2);
for i=1:2
    load(fullfile(Sbj_Metadata.iEEG_data,Sbj_Metadata.BlockLists{i},[Sbj_Metadata.BlockLists{i} '_info.mat']));
    info.events.Accr = num2cell(strcmp(info.events.SyllablePresented,info.events.Response));
    info.events.Accr(strcmp(info.events.Response,'')) = {'noresp'};
    
    consonant = cell(size(info.events,1),1);
    for t=1:size(info.events,1)
        if ismember(info.events.SyllablePresented{t},voiced_syll)
            consonant{t} = 'voiced';
        else
            consonant{t} = 'unvoiced';
        end
        if strcmp(info.events.SyllablePresented{t},info.events.Response{t})
            info.events.Accr{t} = num2str(strcmp(info.events.SyllablePresented{t},info.events.Response{t}));
        elseif strcmp(info.events.Response{t},'')
            info.events.Accr{t} = 'noresp';
        else
            info.events.Accr{t} = '0';
        end
    end
    info.events = [info.events,cell2table(consonant)];
    % remove bad trials in e-stim block (the ones that extend into speaking or the ones that are not stimulating during the stim presentation)
    if strcmp(Sbj_Metadata.sbj_ID,'NS163') && strcmp(Sbj_Metadata.BlockLists{i},'SyllablePerception_EStim')
        info.events = info.events(~ismember(1:64,[1,11,14,16,21,22,25,31,32,33,37,41,43,51,54,58,61,62,63,64]),:);
    end
    if strcmp(Sbj_Metadata.sbj_ID,'NS163') && strcmp(Sbj_Metadata.BlockLists{i},'SyllablePerception')
        info.events = info.events(~ismember(1:92,1:46),:);
    end
    
    voiced_corr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'1'));
    voiced_incorr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'0'));
    voiced_noresp = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'voiced')),'noresp'));
    
    unvoiced_corr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'1'));
    unvoiced_incorr = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'0'));
    unvoiced_noresp = sum(ismember(info.events.Accr(strcmp(info.events.consonant,'unvoiced')),'noresp'));
    
    voiced_total = voiced_corr+voiced_incorr+voiced_noresp;
    unvoiced_total = unvoiced_corr+unvoiced_incorr+unvoiced_noresp;
    block{i}.voiced_corr = voiced_corr;
    block{i}.voiced_incorr = voiced_incorr;
    block{i}.voiced_noresp = voiced_noresp;
    block{i}.voiced_total = voiced_total;
    block{i}.unvoiced_corr = unvoiced_corr;
    block{i}.unvoiced_incorr = unvoiced_incorr;
    block{i}.unvoiced_noresp = unvoiced_noresp;
    block{i}.unvoiced_total = unvoiced_total;
end


figure('Position',[100 100 1200 1000])
to_bar = [block{1}.voiced_corr/block{1}.voiced_total block{1}.voiced_incorr/block{1}.voiced_total block{1}.voiced_noresp/block{1}.voiced_total;...
          block{2}.voiced_corr/block{2}.voiced_total block{2}.voiced_incorr/block{2}.voiced_total block{2}.voiced_noresp/block{2}.voiced_total;...
          block{1}.unvoiced_corr/block{1}.unvoiced_total block{1}.unvoiced_incorr/block{1}.unvoiced_total block{1}.unvoiced_noresp/block{1}.unvoiced_total;...
          block{2}.unvoiced_corr/block{2}.unvoiced_total block{2}.unvoiced_incorr/block{2}.unvoiced_total block{2}.unvoiced_noresp/block{2}.unvoiced_total]*100;
hBar = bar(to_bar,'stacked');
set(gca,'xticklabels',{'Voiced-NoEStim','Voiced-EStim','Unvoiced-NoEStim','Unvoiced-EStim'},'FontSize',15)
grid on
legend('Correct response','Incorrect Response','No-response','Location','southoutside','NumColumns',3)
ctr=[];ydt=[];
for k1 = 1:size(to_bar,2)
    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, hBar(k1).XOffset');    % Note: ‘XOffset’ Is An Undocumented Feature, This Selects The ‘bar’ Centres
    ydt(k1,:) = hBar(k1).YData;                                     % Individual Bar Heights
    if k1==1
        text(ctr(k1,:), ydt(k1,:), sprintfc('%.0f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',18, 'Color','k')
    else
        text(ctr(k1,:), ydt(k1,:)+sum(ydt(1:k1-1,:),1), sprintfc('%.0f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',18, 'Color','k')
    end
end

% Observed data
n1 = block{1}.voiced_corr; N1 = block{1}.voiced_corr+block{1}.voiced_incorr+block{1}.voiced_noresp;
n2 = block{2}.voiced_corr; N2 = block{2}.voiced_corr+block{1}.voiced_incorr+block{1}.voiced_noresp;
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl1,chi2stat1,pval1] = crosstab(x1,x2)

n1 = block{1}.unvoiced_corr; N1 = block{1}.unvoiced_corr+block{1}.unvoiced_incorr+block{1}.unvoiced_noresp;
n2 = block{2}.unvoiced_corr; N2 = block{2}.unvoiced_corr+block{1}.unvoiced_incorr+block{1}.unvoiced_noresp;
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl2,chi2stat2,pval2] = crosstab(x1,x2)
title({['Comparing accuracy of voiced between no-EStim vs EStim, p-val: ' num2str(pval1)];...
       ['Comparing accuracy of unvoiced between no-EStim vs EStim, p-val: ' num2str(pval2)]})

print(fullfile(Sbj_Metadata.results,['2ndhalf_' curr_block '_block_comp_onlyclean_percent.jpg']),'-djpeg','-r300')
% print(fullfile(Sbj_Metadata.results,[curr_block '_block_comp.jpg']),'-djpeg','-r300')
