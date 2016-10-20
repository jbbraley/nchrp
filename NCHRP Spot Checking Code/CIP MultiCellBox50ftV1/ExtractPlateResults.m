%% Extract Plate Results
% For NCHRP spot checking task

clear
clc

%% INPUTS

% Model Path
modelPath = 'C:\Users\DomWirk\Documents\git_jbbraley\nchrp\Models';
modelFile = 'CIPMultiCellBox50ftv1.st7';
fullPath{1} = fullfile(modelPath,modelFile);

% Result Path
resultPath = 'C:\Users\DomWirk\Documents\git_jbbraley\nchrp\Models';
resultFile = 'CIPMultiCellBox50ftv1.LSA';
fullPath{2} = fullfile(resultPath,resultFile);

% Result Case Vector
% For DL, resultCaseNum = 1
% For LL, resultCaseNum = vector of result case numbers that correspond to
% each responsevariable/location.
resultCaseNums = [1];

% plate numbers for each location
% plateNums{Location #} =
% {[Top Flange plate #s];[Web plate #s];[Bottom Flange plate #s]}
% ***Web plates must be orderd top to bottom
plateNums{1} = {[4888 4789 4690 4591 334 2710 2611 2512];[9145 9244 9343 136];[6670 6769 6868 1423 6967 7066 7165 7264]}; % mid M+1
plateNums{2} = {[4294 4195 4096 3997 829 3898 3799 3700];[9442 9541 9640 1819];[7858 7957 8056 1621 8452 8551]};          % 2nd M+1
%plateNums{3} = {[3601 3502 3403 235 8848 8749];[433 10927 10828 10729 10630];[8650 1720]};                                % ext M+1
plateNums{3} = {[4916 4817 4718 4619 362 2738 2639 2540];[9173 9272 9371 164];[6698 6797 6896 1451 6995 7094 7193 7292]}; % mid M-
plateNums{4} = {[4322 4223 4124 4025 857 3926 3827 3728];[9470 9569 9668 1847];[7886 7985 8084 1649 8480 8579]};          % 2nd M-
%plateNums{6} = {[3629 3530 3431 263 8876 8777];[461 10955 10856 10757 10658];[8678 1748]};                                % ext M-
plateNums{5} = {[4943 4844 4745 4646 389 2765 2666 2567];[9200 9299 9398 191];[6725 6824 6923 1478 7022 7121 7220 7319]}; % mid M+2
plateNums{6} = {[4349 4250 4151 4052 884 3953 3854 3755];[9497 9596 9695 1874];[7912 8012 8111 1676 8507 8606]};          % 2nd M+2
%plateNums{9} = {[3656 3557 3458 290 8903 8804];[488 10982 10883 10784 10685];[8705 1775]};                                % ext M+2

% Average Length of Elements
TopFlangePlateLength = 12;
WebPlateLength = 13;
BottomFlangePlateLength = 12;

% Moment arms for each element
TopFlange_y = -26;
Web_y = [-13*1.5 -13*0.5 13*0.5 13*1.5];
BottomFlange_y = 26;

%% Setup/Housekeeping
global rtPlateMoment stPlateGlobal AtCentroid psPlateMidPlane rtPlateForce rtPlateStress
resultType = {rtPlateMoment;rtPlateForce;rtPlateStress};
resultSubType = stPlateGlobal;
sampleLocation = AtCentroid;
surface = psPlateMidPlane;
layer = 1;
NumPoints = 1;
plateResults = cell(1,length(plateNums));

for ii = 1:length(plateNums)

    numTopFlangePlates = length(plateNums{ii}{1});
    numWebPlates = length(plateNums{ii}{2});
    numBottomFlangePlates = length(plateNums{ii}{3});

    % Average (approximate) length for each plate
    l{ii} = {TopFlangePlateLength*ones(1,numTopFlangePlates);...
        WebPlateLength*ones(1,numWebPlates);...
        TopFlangePlateLength*ones(1,numBottomFlangePlates)};

    % moment arms for each plate
    y{ii} = {TopFlange_y*ones(1,numTopFlangePlates);...
        Web_y;...
        BottomFlange_y*ones(1,numBottomFlangePlates)};
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
            if length(resultCaseNums) == 1 % Single result case
                resultCaseNum = resultCaseNums;
            else
                resultCaseNum = resultCaseNums(ii);
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

% Close Model & Result File, Unload API
CloseAndUnload(uID);

%% Filter Results, Calculate total moment

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

TotalMoment = sum(MomentF+MomentM,1); % Total Absolute Moment

clearvars -except Stress TotalMoment MomentF MomentM plateNums l y...
    plateResults modelFile resultFile resultCaseNums

uisave(who,resultFile(1:end-4))
