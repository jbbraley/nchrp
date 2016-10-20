function St7OpenModelFile(uID, fullPath, ScratchPath)

iErr = calllib('St7API', 'St7OpenFile', uID, fullPath, ScratchPath);
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