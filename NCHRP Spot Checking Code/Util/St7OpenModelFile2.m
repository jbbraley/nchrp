function St7OpenModelFile2(uID, modelPath, ScratchPath)

iErr = calllib('St7API', 'St7OpenFile', uID, modelPath, ScratchPath);
if iErr == 1
    try
        % close file
        iErr = calllib('St7API','St7CloseFile', uID);
        HandleError(iErr);
    catch
    end
elseif iErr ~= 0 
    HandleError(iErr);
end

end