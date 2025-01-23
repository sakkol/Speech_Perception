% Calculating distance of cathode/anode to STG

%% NS144_02 (S1)
anode1 = [-47.982000 -3.559000 65.252200]; % Ref1
anode2 = [-55.617000 -2.411000 58.084100]; % Ref2

cathode1 = [-69.136000 -1.395000 23.782000]; % Ref6
cathode2 = [-71.075000 -0.851000 13.739000]; % Ref7

% calculate mid points of stim electrodes
mid_anode = mean([anode1;anode2],1);
mid_cathode = mean([cathode1;cathode2],1);

% distances to point of interest
POI = [-67.762657 -26.354834 13.395566]; % LGrid22

dist_anode = norm(mid_anode - POI)
dist_cathode = norm(mid_cathode - POI)

distances = [cell2table({'NS144_02'},"VariableNames",{'Session'}),array2table([dist_anode,dist_cathode],'VariableNames',{'dist_anode','dist_cathode'})];

%% NS148 (S2)
anode = [62.000000 5.000000 -5.000000]; % RTs_bolt

cathode = [37.850000 26.850000 56.300000]; % RFa_bolt

% distances to point of interest
POI = [55.333300 2.000000 -4.667000]; % RTs16

dist_anode = norm(anode - POI)
dist_cathode = norm(cathode - POI)

distances = [distances; [cell2table({'NS148'},"VariableNames",{'Session'}),array2table([dist_anode,dist_cathode],'VariableNames',{'dist_anode','dist_cathode'})]];

%% NS148_02 (S3)
anode1 = [57.742400 13.318000 38.151500]; % Ref8Blue6
anode2 = [58.868000 15.139000 28.125000]; % Ref8Blue7

cathode1 = [2.524000 -36.905000 69.000000]; % Ref4Blue2
cathode2 = [8.028000 -44.972000 65.777800]; % Ref4Blue3

% calculate mid points of stim electrodes
mid_anode = mean([anode1;anode2],1);
mid_cathode = mean([cathode1;cathode2],1);

% distances to point of interest
POI = [62.416939 -15.737971 4.548935]; % RTi5

dist_anode = norm(mid_anode - POI)
dist_cathode = norm(mid_cathode - POI)

distances = [distances; [cell2table({'NS148_02'},"VariableNames",{'Session'}),array2table([dist_anode,dist_cathode],'VariableNames',{'dist_anode','dist_cathode'})]];

%% NS150 (S4)
anode = [38.9582  -48.0964   53.7647]; % RIp_bolt1

cathode = [81.9089   -5.2299   -4.8790]; % RDh_bolt1

% distances to point of interest
POI = [61.973000 -0.490000 8.128000]; % RTs14

dist_anode = norm(anode - POI)
dist_cathode = norm(cathode - POI)

distances = [distances; [cell2table({'NS150'},"VariableNames",{'Session'}),array2table([dist_anode,dist_cathode],'VariableNames',{'dist_anode','dist_cathode'})]];

%% NS170 (S5)
anode1 = [-70.146000 2.922000 -24.204000]; % Ref2_1
anode2 = [-71.968000 -6.075000 -20.075000]; % Ref2_2

cathode1 = [-72.094000 -23.812000 -11.553000]; % Ref2_4
cathode2 = [-70.737000 -32.579000 -6.958000]; % Ref2_5

% calculate mid points of stim electrodes
mid_anode = mean([anode1;anode2],1);
mid_cathode = mean([cathode1;cathode2],1);

% distances to point of interest
POI = [-63.324928 -19.807814 -6.132242]; % LGrid38

dist_anode = norm(mid_anode - POI)
dist_cathode = norm(mid_cathode - POI)

distances = [distances; [cell2table({'NS170'},"VariableNames",{'Session'}),array2table([dist_anode,dist_cathode],'VariableNames',{'dist_anode','dist_cathode'})]];

%% NS174_02 (S6)
anode1 = [-70.043000 -1.006000 -9.784000]; % Ref1
anode2 = [-70.019000 -7.115000 -1.683000]; % Ref2

cathode1 = [-72.593000 -17.074000 25.241000]; % Ref5
cathode2 = [-67.726000 -19.630000 33.205500]; % Ref6

% calculate mid points of stim electrodes
mid_anode = mean([anode1;anode2],1);
mid_cathode = mean([cathode1;cathode2],1);

% distances to point of interest
POI = [-62.571000 -4.857000 5.286000]; % LTs16

dist_anode = norm(mid_anode - POI)
dist_cathode = norm(mid_cathode - POI)

distances = [distances; [cell2table({'NS174_02'},"VariableNames",{'Session'}),array2table([dist_anode,dist_cathode],'VariableNames',{'dist_anode','dist_cathode'})]];


%%
%%
%% Interpolating
RTs15 = [53.666600 1.333500 -4.417000];
RTs16 = [55.333300 2.000000 -4.667000];
RTsbolt = [62.000000 5.000000 -5.000000];
k = (RTsbolt-RTs16) ./ (RTs16-RTs15);

last_one = [62.142800 -3.785000 -6.925000];
last_minone = [57.201200 -3.464000 -8.461000];
last_one + (k .* (last_one-last_minone))

