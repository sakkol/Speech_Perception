function IL_color_elecs_wData(Sbj_Metadata,data,AllSubElecNames,AllSubElecCoords,contVSdiscr,MinMax,ColorBarTitle,alleleccols,edgeColoredElecs)
% This function plots data from electrodes on individual or fsaverage
% brain with a colorbar next to it. 
% Sbj_Metadata can be a struct, 'fsaverage' or a directory. For better
% understanding, first part of the code is pretty easy tu understand:
%       If a struct, 'fsaverage_Dir' field is directory to fsaverage
%       Freesurfer folder (e.g.'/media/sakkol/HDD1/HBML/DERIVATIVES/freesurfer/fsaverage')
%       and 'fsname' field is subject's Freesurfer name (e.g. 'NS120_02').
%       If 'fsaverage', default Freesurfer directory will be used (taken from shell).
%       If a directory, it will be accepted as Freesurfer main directory.
% data is either an array or a matrix of data where each row is the data
%       for each electrode. If matrix, data will be averaged in columns (so
%       that there will be 1 value for each electrode. If data is a cell,
%       it will be assumed that these are labels of each electrode in
%       discrete fashion.
% AllSubElecNames and AllSubElecCoords are needed if average brain or if
%       not all electrodes are not used for the subject and can be obtained
%       with this function:
%       [AllSubElecNames,AllSubElecCoords] = gather_allelecsinfo(sbjs_elecs);
% contVSdiscr defines if the colors of data should be given as continuous
%       colorbar (if data is [0.1 0.3 1.2 -2]) or discrete legend points
%       (if data is [1 2 5 4]).
% MinMax defines the higher and lower boundaries of the plot so that the
%       colors won't exceed these limits. Lower or higher values will be
%       assigned with min or max values. (e.g. [2 15])
% ColorBarTitle defines what the colorbar title will be (possible for
%       continuous data). (e.g. 't-values')
% ColorMaps is to select which colormap is going to be used. Default value
%       is inferno from matplotlib.
% edgeColoredElecs is a cell matrix (nx2) that contains electrode names in
%       first column and the color of the edge of the electrodes in the
%       second column. If a list of electrodes are given, those electrodes'
%       edges will be green. Edges whose electrodes that are not present in
%       the list will be same as the electrodes'. This option can be used
%       to show significant electrodes.
%
% Additional notes:
% - There will be no title, but it can be added using this line of code:
% text(gca,.5,1.07,'Your beautiful title','Units','normalized','FontSize',18,'FontWeight','bold','HorizontalAlignment','center')
% 
% Example use:
% sbjs_elecs{1,1} = Sbj_Metadata.fsname;
% sbjs_elecs{1,2} = ecog_avg.ftrip.label;
% [AllSubElecNames,AllSubElecCoords] = gather_allelecsinfo(sbjs_elecs,'LEPTO');
% edgeColoredElecs = ecog_avg.ftrip.label(pval < (.001/length(ecog_avg.ftrip.label)));
% color_elecs_wData(Sbj_Metadata,tval(ismember(ecog_avg.ftrip.label,AllSubElecNames))',AllSubElecNames,AllSubElecCoords,'continuous',[3 18],'t-values','inferno',edgeColoredElecs)
% text(gca,.5,1.07,'t-values of HFA comparing silent and sound periods during Auditory Localizer','Units','normalized','FontSize',18,'FontWeight','bold','HorizontalAlignment','center')
% print(fullfile(Sbj_Metadata.results,[Sbj_Metadata.sbj_ID '_tvals_brain_inferno.jpg']),'-djpeg','-r300')
%
% 
% Serdar Akkol, May 2020, HBML
%
% Updates:
% - MinMax and ColorBarTitle options are added for improved functionality.
% (July 2020). 
% - ColorMaps option added with new colormaps. crameri.m and
% select_colormaps.mat files are needed for these. (December 2020)
% - edgeColoredElecs option added. (December 2020)

%% Set default inputs
if isstruct(Sbj_Metadata)
    fsDir = erase(Sbj_Metadata.fsaverage_Dir,'fsaverage');
    sub = Sbj_Metadata.fsname;
elseif strcmp(Sbj_Metadata,'fsaverage')
    fsDir =  getFsurfSubDir();
    sub = 'fsaverage';
else  % for example if Sbj_Metadata is a directory
    fsDir = Sbj_Metadata;
    sub = 'fsaverage';
end
if ~exist('contVSdiscr','var') || isempty(contVSdiscr)
    contVSdiscr = 'continuous';
end
if ~exist('ColorBarTitle','var') || isempty(ColorBarTitle)
    ColorBarTitle = 'Data range';
end

%% Arrange individual brain electrode coordinates
if ~strcmp(sub,'fsaverage') && (~exist('AllSubElecNames','var') || isempty(AllSubElecNames))
    % Find the locations of electrodes
    elecInfoFname=fullfile(fsDir,sub,'elec_recon',[sub '.electrodeNames']);
    elecInfo=csv2Cell(elecInfoFname,' ',2);
    elec_toplot_names = elecInfo(:,1);
    if isempty(AllSubElecNames) % if the electrode names are not given as input.
        AllSubElecNames = elec_toplot_names;
        if length(AllSubElecNames) ~= size(data,1)
            error('There should be equal number of electrodes and data to plot. Please provide an input for electrode names.')
        end
    end
    
    % Load electrode coordinates from LEPTO
    elecCoordFname=fullfile(fsDir,sub,'elec_recon',[sub '.LEPTO']);
    if ~exist(elecCoordFname,'file'),error('No LEPTO for %s\n',sub),end
    elecCoordCsv=csv2Cell(elecCoordFname,' ',2);
    fullElecCoord=zeros(size(elecCoordCsv));
    
    % Convert coordinates from string to #
    for a=1:size(elecCoordCsv,1)
        for b=1:3
            fullElecCoord(a,b)=str2double(elecCoordCsv{a,b});
        end
    end
    
    clear allelecCoords allelecnames
    AllSubElecCoords(size(AllSubElecNames,1),1:3) = 0;
    allelecnames{size(AllSubElecNames,1),1} = '';
    for el = 1:size(AllSubElecNames,1)
        AllSubElecCoords(el,1:3) = fullElecCoord(strcmpi(elec_toplot_names,AllSubElecNames{el}),:);
        AllSubElecCoords(el,4) = strcmpi(elecInfo{strcmpi(elec_toplot_names,AllSubElecNames{el}),3},'L');
        allelecnames{el,1} = AllSubElecNames{el};
    end
end

%% Set color for each data point using the selected colormap
if size(data,1) == 1 && (size(data,2) == 1 || size(data,2) == 3)
    error('There is only 1 data point for only 1 electrode, quiting!')
elseif iscell(data)
    contVSdiscr = 'discrete';
    for d = 1:length(data)
        data{d} = char(data{d});
    end
    uni_data = unique(data);
    x = ones(length(data),1);
    for d = 1:length(uni_data)
        x(strcmp(data,uni_data(d)))=d;
    end
elseif size(data,1) == 1 && size(data,2) > 1
    x = data';
elseif size(data,1) > 1 && size(data,2) > 1
    x = zeros(size(data,1),1);
    for i = 1:size(data,1) % if data is 2D, this will average it in columns
        x(i,1) = mean(data(i,:),2); %data to be plotted
    end
    fprintf('You may get an error if matrix is in data x elecs size.\n')
else
    x=data;
end
% check data length
if size(x,1) ~= length(AllSubElecNames)
    error('Numbers of data and electrode names do not match. Please check the input.')
end

% create edge colors
if ~exist('edgeColoredElecs','var') || isempty(edgeColoredElecs)
    edgeColoredElecs = [AllSubElecNames,num2cell(alleleccols,2)];
    markerVSsphere='sphere';
    elecSize=2;
elseif exist('edgeColoredElecs','var') && size(edgeColoredElecs,2) == 1
    edgeColoredElecs(:,2) = {'g'};
    markerVSsphere='marker';
    elecSize=13;
else
    markerVSsphere='marker';
    elecSize=13;
end
alledgecolors = cell(length(AllSubElecNames),2);
for el = 1:length(AllSubElecNames)
    alledgecolors{el,1} = AllSubElecNames{el};
    if ~isempty(edgeColoredElecs) && ismember(AllSubElecNames{el},edgeColoredElecs(:,1))
        alledgecolors{el,2} = edgeColoredElecs{ismember(edgeColoredElecs(:,1),AllSubElecNames{el}),2};
    else
        alledgecolors{el,2} = alleleccols(el,:);
    end
end

%% Start plotting
% arrange cfg and plotPialSurf
cfg=[];
if any(AllSubElecCoords(:,4)) && any(~AllSubElecCoords(:,4)) % 1 means Left, 0 means Right
    cfg.view='omni';
    fprintf('\n\tPlotting data of %d electrodes in both hemispheres.\n',length(AllSubElecNames))
elseif any(AllSubElecCoords(:,4)) && ~any(~AllSubElecCoords(:,4))
    cfg.view='lomni';
    fprintf('\n\tPlotting data of %d electrodes only in left hemispheres.\n',length(AllSubElecNames))
elseif ~any(AllSubElecCoords(:,4)) && any(~AllSubElecCoords(:,4))
    cfg.view='romni';
    fprintf('\n\tPlotting data of %d electrodes only in right hemispheres.\n',length(AllSubElecNames))
end
cfg.elecShape=markerVSsphere;
cfg.elecColors = alleleccols;
cfg.elecCbar='n';
cfg.elecCoord = AllSubElecCoords;
cfg.elecNames = AllSubElecNames;
cfg.elecSize = elecSize;
if ~isempty(alledgecolors),cfg.edgeColors = alledgecolors;end
cfg.title='';
cfg.ignoreDepthElec = 'n';
cfg.pullOut = 0;
cfg.showLabels = 'n';
cfg.opaqueness = .25;
cfg.fsurfSubDir = fsDir;
% cfg.surfType = 'inflated';
plotPialSurf(sub,cfg);

% general settings
set(gcf,'Units','normalized');
set(gcf,'Position', [0 0 1 1]);
set(gcf,'color','w')
set(gcf, 'InvertHardcopy', 'off') % to set saving white color as white

if strcmp(contVSdiscr,'discrete')
    % Create and delete new axes to plot discrete data as a legend
    ax=axes;
    uni_x=unique(x);
    for i=1:length(uni_x)
        tmpcol = alleleccols(x==uni_x(i),:);
        plot(uni_x(i),'o','MarkerFaceColor',tmpcol(1,:),'MarkerEdgeColor',tmpcol(1,:),'MarkerSize',20)
        hold on
    end
    xlim([-11550 -11549])
    ylim([-11550 -11549])
    
    if ~exist('uni_data','var')
        l=legend(cellfun(@num2str,num2cell(uni_x),'UniformOutput',false),'Interpreter','none');
    else
        l=legend(cellfun(@num2str,uni_data,'UniformOutput',false),'Interpreter','none');
    end
    l.FontSize=20;
    l.LineWidth = 1;
    l.Location = 'southoutside';
    l.Color = [0.9412    0.9412    0.9412]; % change background color so that whites would be seen
    
    a=get(l); %gets properties of colorbar
    a = a.Position; %gets the positon and size of the legend
    set(l,'Position',[a(1) a(2)-0.07 a(3) a(4)])% To change position of legend: put it right out size
    ax.Visible = 'off';
    
else
    error('There is something wrong with ''contVSdiscr'' input')
end


end
