



















fsSub=subj_list{s};

% Grab electrode names and hemisphere from subject dir
elecInfoFname=fullfile(fsDir,fsSub,'elec_recon',[fsSub '.electrodeNames']);
elecInfo=csv2Cell(elecInfoFname,' ',2);
elec_toplot_names = elecInfo(:,1);

% Load electrode coordinates from LEPTO
elecCoordFname=fullfile(fsDir,fsSub,'elec_recon',[fsSub '.LEPTO']);
elecCoordCsv=csv2Cell(elecCoordFname,' ',2);
fullElecCoord=zeros(size(elecCoordCsv));
% Convert coordinates from string to #
for a=1:size(elecCoordCsv,1),
    for b=1:3,
        fullElecCoord(a,b)=str2double(elecCoordCsv{a,b});
    end
end