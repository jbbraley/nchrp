function St7CloseResultFile(uID)

try
iErr= calllib('St7API', 'St7CloseResultFile', uID);
HandleError(iErr);
catch
end
    
end %St7CloseResultFile