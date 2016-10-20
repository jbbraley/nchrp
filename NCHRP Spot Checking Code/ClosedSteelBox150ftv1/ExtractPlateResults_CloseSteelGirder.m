%% Extract Plate Results
% For NCHRP spot checking task

clear
clc

%% INPUTS

% Model Path
modelPath = 'C:\Users\DomWirk\Documents\git_jbbraley\nchrp\Models';
modelFile = 'ClosedSteelBox150ftv1.st7';
fullPath{1} = fullfile(modelPath,modelFile);

% Result Path
resultPath = 'C:\Users\DomWirk\Documents\git_jbbraley\nchrp\Models';
resultFile = 'ClosedSteelBox150ftv1.LSA';
fullPath{2} = fullfile(resultPath,resultFile);

% Result Case Vector
% For DL, resultCaseNum = 1
% For LL, resultCaseNum = vector of result case numbers that correspond to
% each responsevariable/location. First Row is moment result cases, second
% row is shear result cases
resultCaseNums = 29;

% Support Node numbers for reactions at supports
supportNodes(1) = 67;
supportNodes(2) = 70;
supportNodes(3) = 15511;
supportNodes(4) = 16108;
supportNodes(5) = 15611;
supportNodes(6) = 16208;

% Specify which locations are over the pier
pierLocations = [3 4];

% plate numbers for each location
% plateNums{Location #} =
% {[Deck plate #s];[Haunch plate #s];[Top Flange plate #s];[Web 1 plate #s];[Web plate 2 #s];[Bottom Flange plate #s]}
% ***Web plates must be orderd top to bottom
plateNums{1} = {[17694 1774 14311 14510 14709 1973 17893 18092 18291];...% [Deck plate #s]      10530   int M+1 case 12
                [3963 4162 4361 580];...% [Haunch plate #s]
                [7147 7346 7545 5555];...% [Top Flange plate #s]
                [12321 24261 24062 23863 23664 23465];...% [Web 1 plate #s]
                [12520 25256 25057 24858 24659 24460];...% [Web 2 plate #s]
                [10331 10530 10729 8739]}; % [Bottom Flange plate #s]
plateNums{2} = {[18490 2172 14908 15107 15306 2371 18689 18888 19087 19286 2570];...% [Deck plate #s]      11127   ext M+1 case 14
                [4560 4759 4958 779];...% [Haunch plate #s]
                [7744 7943 8142 5754];...% [Top Flange plate #s]
                [12719 26251 26052 25853 25654 25455];...% [Web 1 plate #s]
                [12918 27246 27047 26848 26649 26450];...% [Web 2 plate #s]
                [10928 11127 11326 8938]}; % [Bottom Flange plate #s]
plateNums{3} = {[17748 1828 14365 14564 14763 2027 17947 18146 18345];...% [Deck plate #s]     10584 int M- case 4
                [4017 4216 4415 634];...% [Haunch plate #s]
                [7201 7400 7599 5609];...% [Top Flange plate #s]
                [12375 24315 24116 23917 23718 23519];...% [Web 1 plate #s]
                [12574 25310 25111 24912 24713 24514];...% [Web 2 plate #s]
                [10385 10584 10783 8793]}; % [Bottom Flange plate #s]
plateNums{4} = {[18544 2226 14962 15161 15360 2425 18743 18942 19141 19340 2624];...% [Deck plate #s]         11181 ext M- case 5
                [4614 4813 5012 833];...% [Haunch plate #s]
                [7798 7997 8196 5808];...% [Top Flange plate #s]
                [12773 26305 26106 25907 25708 25509];...% [Web 1 plate #s]
                [12972 27300 27101 26902 26703 26504];...% [Web 2 plate #s]
                [10982 11181 11380 8992]}; % [Bottom Flange plate #s]
plateNums{5} = {[17803 1883 14420 14619 14818 2082 18002 18201 18400];...% [Deck plate #s]         10639 int M+2 case 13
                [4072 4217 4470 689];...% [Haunch plate #s]
                [7256 7455 7654 5664];...% [Top Flange plate #s]
                [12430 24370 24171 23972 23773 23574];...% [Web 1 plate #s]
                [12629 25365 25166 24967 24768 24569];...% [Web 2 plate #s]
                [10440 10639 10838 8848]}; % [Bottom Flange plate #s]
plateNums{6} = {[18599 2281 15017 15216 15415 2480 18798 18997 19196 19395 2679];...% [Deck plate #s]          11236 ext M+2 case 15
                [4669 4868 5067 888];...% [Haunch plate #s]
                [7853 8052 8251 5863];...% [Top Flange plate #s]
                [12828 26360 26161 25962 25763 25564];...% [Web 1 plate #s]
                [13027 27355 27156 26957 26758 26559];...% [Web 2 plate #s]
                [11037 11236 11435 9047]}; % [Bottom Flange plate #s]


% Specify which locations are exterior web
extLocations = [2 4 6];

% Average Length and Thickness of Elements
DeckPlateLength = 13;
DeckPlateThickness = 9.5;
HaunchPlateLength = 13;
HaunchPlateThickness = 3.5;
TopFlangePlateLength = 13;
TopFlangePlateThickness = 2.375;
WebPlateLength = 12;
WebPlateThickness = 0.625;
BottomFlangePlateLength = 13;
BottomFlangePlateThickness = 2.375;

% Web Angle
intWebAngle = 0;
extWebAngle = 0;


%% Setup/Housekeeping
global rtPlateMoment stPlateGlobal AtCentroid psPlateMidPlane rtPlateForce rtPlateStress rtNodeReact
resultType = {rtPlateMoment;rtPlateForce;rtPlateStress;rtNodeReact};
resultSubType = stPlateGlobal;
sampleLocation = AtCentroid;
surface = psPlateMidPlane;
layer = 1;
NumPoints = 1;
plateResults = cell(1,length(plateNums));


for ii = 1:length(plateNums)

    % Get Number of elements
    numDeckPlates = length(plateNums{ii}{1});
    numHaunchPlates = length(plateNums{ii}{2});
    numTopFlangePlates = length(plateNums{ii}{3});
    numLeftWebPlates = length(plateNums{ii}{4});
    numRightWebPlates = length(plateNums{ii}{5});
    numBottomFlangePlates = length(plateNums{ii}{6});

    % Average (approximate) length for each plate
    l{ii} = {DeckPlateLength*ones(1,numDeckPlates);...
        HaunchPlateLength*ones(1,numHaunchPlates);...
        TopFlangePlateLength*ones(1,numTopFlangePlates);...
        WebPlateLength*ones(1,numLeftWebPlates);...
        WebPlateLength*ones(1,numRightWebPlates);...
        TopFlangePlateLength*ones(1,numBottomFlangePlates)};

    % Get Web Angle
    if  any(ii == extLocations)
        webAngle = extWebAngle;
    else
        webAngle = intWebAngle;
    end

    % Find centroid of each cross-section
    DeckArea = numDeckPlates*DeckPlateLength*DeckPlateThickness;
    DeckY = 0.5*DeckPlateThickness + HaunchPlateThickness +  TopFlangePlateThickness...
        + numRightWebPlates*WebPlateLength*cosd(webAngle) + BottomFlangePlateThickness;

    HaunchArea = numHaunchPlates*HaunchPlateLength*HaunchPlateThickness;
    HaunchY = 0.5*HaunchPlateThickness +  TopFlangePlateThickness...
        + numRightWebPlates*WebPlateLength*cosd(webAngle) + BottomFlangePlateThickness;

    TopFlangeArea = numTopFlangePlates*TopFlangePlateLength*TopFlangePlateThickness;
    TopFlangeY = 0.5*TopFlangePlateThickness + numRightWebPlates*WebPlateLength*cosd(webAngle) + BottomFlangePlateThickness;

    RightWebArea = numRightWebPlates*WebPlateLength*WebPlateThickness;
    RightWebY = numRightWebPlates*WebPlateLength*cosd(webAngle) + BottomFlangePlateThickness;

    LeftWebArea = numLeftWebPlates*WebPlateLength*WebPlateThickness;
    LeftWebY = numLeftWebPlates*WebPlateLength*cosd(webAngle) + BottomFlangePlateThickness;

    BottomFlangeArea = numBottomFlangePlates*BottomFlangePlateLength*BottomFlangePlateThickness;
    BottomFlangeY = 0.5*BottomFlangePlateThickness;

    yBar = (DeckArea*DeckY + HaunchArea*HaunchY + TopFlangeArea*TopFlangeY...
        + RightWebArea*RightWebY + LeftWebArea*LeftWebY + BottomFlangeArea*BottomFlangeY)...
        /(DeckArea+HaunchArea+TopFlangeArea+RightWebArea+LeftWebArea+BottomFlangeArea);

    % moment arms for each plate
    for jj = 1:length(plateNums{ii}) % For each cross-section component

        if jj == 1 % Deck
            y{ii}{jj} = (yBar - DeckY) * ones(1,numDeckPlates);
        elseif jj == 2 % Haunch
            y{ii}{jj} = (yBar - HaunchY) * ones(1,numHaunchPlates);
        elseif jj == 3 % Top Flange
            y{ii}{jj} = (yBar - TopFlangeY) * ones(1,numTopFlangePlates);
        elseif jj == 4 || jj == 5 % Web
            for kk = 1:size(plateNums{ii}{jj},2) % For each individual plate
                vect(kk) = yBar - (kk*WebPlateLength*cosd(webAngle)) + 0.5*BottomFlangePlateThickness;
            end
            y{ii}{jj} = vect(size(plateNums{ii}{jj},2):-1:1); % Rearrange so order is top to bottom
        elseif jj == 6 % Bottom Flange
            y{ii}{jj} = (yBar - BottomFlangeY) * ones(1,numBottomFlangePlates);
        end

    end

end

%% Initalize API
St7Start = 1;
ScratchPath = 'C:\Temp';
uID = 1;
InitializeSt7();

%% Open Model & Result Files
St7OpenModelFile(uID, fullPath{1}, ScratchPath);
St7OpenResultFile(uID, fullPath{2});

%% Get plate element results

for ii = 1:length(plateNums) % For each location of interest

    for jj = 1:length(plateNums{ii}) % For each cross-section component

        % Pre-allocate
        % Double 6xNumber of plates for cross-section component
        PlateRes = zeros(3,size(plateNums{ii}{jj},2));

        for kk = 1:size(plateNums{ii}{jj},2) % For each individual plate

            % Pre-allcoate for individual plate results
            pRes = zeros(6,1);

            % Get single shell number
            plateNum = plateNums{ii}{jj}(kk);

            % Get result case number
            if size(resultCaseNums,2) > 1 % Multiple result cases
                resultCaseNum = resultCaseNums(2,ii);
            else % Single result case
                resultCaseNum = resultCaseNums;
            end

            % Get moment response for individual plate
            [iErr,NumPoints,~,pRes] = calllib('St7API',...
                'St7GetPlateResultArray', uID, resultType{1}, resultSubType,...
                plateNum, resultCaseNum, sampleLocation, surface,...
                layer, NumPoints, 6, pRes);
            HandleError(iErr);

            % assign
            PlateRes(1,kk) = pRes(1);

            % Get force response for individual plate
            [iErr,NumPoints,~,pRes] = calllib('St7API',...
                'St7GetPlateResultArray', uID, resultType{2}, resultSubType,...
                plateNum, resultCaseNum, sampleLocation, surface,...
                layer, NumPoints, 6, pRes);
            HandleError(iErr);

            % assign
            PlateRes(2,kk) = pRes(1);

            % Get stress for individual plate
            [iErr,NumPoints,~,pRes] = calllib('St7API',...
                'St7GetPlateResultArray', uID, resultType{3}, resultSubType,...
                plateNum, resultCaseNum, sampleLocation, surface,...
                layer, NumPoints, 6, pRes);
            HandleError(iErr);

            % assign
            PlateRes(3,kk) = pRes(1);

        end

        % assign
        plateResults{ii}{jj,1} = PlateRes;


    end
end
%% Get Support node reactions
for ii = 1:length(supportNodes) % For each location of interest

    % Pre-allcoate for individual node results
    nodeRes = zeros(6,1);

    % Get single node number
    supportNode = supportNodes(ii);

    % Get result case number
    if size(resultCaseNums,2) > 1 % Multiple result cases
        resultCaseNum = resultCaseNums(2,ii);
    else % Single result case
        resultCaseNum = resultCaseNums;
    end

    [iErr, nodeRes] = calllib('St7API',...
            'St7GetNodeResult', uID, resultType{4}, supportNode,...
            resultCaseNum, nodeRes);

    nodeResults(:,ii) = nodeRes;

end

%% Close Model & Result File, Unload API
CloseAndUnload(uID);

%% Calculate total moment and fitler results

for ii = 1:length(plateResults) % For each location of interest
    for jj = 1:length(plateResults{ii}) % For each cross-section component

        % Moment from axial forces
        MomentF(jj,ii) = sum(abs(plateResults{ii}{jj}(2,:).*l{ii}{jj}.*y{ii}{jj}));

        % Moment in flanges
        if jj ~= 2 %(exclude web)
            MomentM(jj,ii) = sum(abs(plateResults{ii}{jj}(1,:).*l{ii}{jj}));
        end

        % Stress [Top Flange; Web; Bottom Flange]
        [~,ind] = max(abs(plateResults{ii}{jj}(3,:)));
        Stress(jj,ii) = max(abs(plateResults{ii}{jj}(3,:))).*sign(plateResults{ii}{jj}(3,ind));

    end
end

% Shear
Shear = nodeResults(3,:);
Shear(pierLocations) = Shear(pierLocations)/2;

% Total Moment
TotalMoment = sum(MomentF+MomentM,1); % Total Absolute Moment

clearvars -except Stress TotalMoment MomentF MomentM plateNums l y...
    plateResults modelFile resultFile resultCaseNums Shear nodeResults

uisave(who,resultFile(1:end-4))
