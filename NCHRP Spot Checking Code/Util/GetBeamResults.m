function beamResults = GetBeamResults(modelPath,resultPath,beamNums,resultCaseNum)

% Load St7 API
uID = 1;
ScratchPath = 'C:\Temp';
InitializeSt7(); 

% Open model file
St7OpenModelFile2(uID, modelPath, ScratchPath)

% Open result file
St7OpenResultFile(uID, resultPath)

% Get Beam Element Results
beamResults = St7GetBeamEndResults(uID, beamNums, resultCaseNum);

% Close and UNload
CloseAndUnload(uID);
    
% % Filter Beam Results
% M1 = max(abs(beamResults(4,:)),abs(beamResults(10,:)));
% P = max(abs(beamResults(1,:)),abs(beamResults(7,:)));
% V = max(abs(beamResults(5,:)),abs(beamResults(9,:)));


end