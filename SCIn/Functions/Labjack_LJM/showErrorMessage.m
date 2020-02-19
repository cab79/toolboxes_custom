function showErrorMessage(e)
% SHOWERRORMESSAGE  Displays the LJM or .NET error from a MATLAB exception.

if(isa(e, 'NET.NetException'))
    eNet = e.ExceptionObject;
    if(isa(eNet, 'LabJack.LJM.LJMException'))
        disp(['LJM Error: ' char(eNet.ToString())])
    else
        disp([char(class(eNet)) ': ' char(eNet.ToString())])
    end
end
disp(getReport(e))
end  % showErrorMessage end
