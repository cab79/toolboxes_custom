%
% Demonstrates how to use the eWriteNames (LJM_eWriteNames) function using .NET.
%
% support@labjack.com
%

clc  % Clear the MATLAB command window
clear  % Clear the MATLAB variables

% Make the LJM .NET assembly visible in MATLAB
ljmAsm = NET.addAssembly('LabJack.LJM');

% Creating an object to nested class LabJack.LJM.CONSTANTS
t = ljmAsm.AssemblyHandle.GetType('LabJack.LJM+CONSTANTS');
LJM_CONSTANTS = System.Activator.CreateInstance(t);

handle = 0;

try
    % Open first found LabJack

    % Any device, Any connection, Any identifier
    [ljmError, handle] = LabJack.LJM.OpenS('ANY', 'ANY', 'ANY', handle);

    % T7 device, Any connection, Any identifier
    % [ljmError, handle] = LabJack.LJM.OpenS('T7', 'ANY', 'ANY', handle);

    % T4 device, Any connection, Any identifier
    % [ljmError, handle] = LabJack.LJM.OpenS('T4', 'ANY', 'ANY', handle);

    % Any device, Any connection, Any identifier
    % [ljmError, handle] = LabJack.LJM.Open(LJM_CONSTANTS.dtANY, ...
    %     LJM_CONSTANTS.ctANY, 'ANY', handle);

    showDeviceInfo(handle);

    % Setup and call eWriteNames to write values.
    numFrames = 2;
    aNames = NET.createArray('System.String', numFrames);
    aNames(1) = 'DAC0';
    aNames(2) = 'TEST_UINT16';
    aValues = NET.createArray('System.Double', numFrames);
    aValues(1) = 2.5;  % 2.5 V
    aValues(2) = 12345;  % 12345
    LabJack.LJM.eWriteNames(handle, numFrames, aNames, aValues, 0);

    disp('eWriteNames:');
    for i=1:numFrames,
        disp(['  Name: ' char(aNames(i)) ', Value: ' num2str(aValues(i))])
    end
catch e
    showErrorMessage(e)
    LabJack.LJM.CloseAll();
    return
end

try
    % Close handle
    LabJack.LJM.Close(handle);
catch e
    showErrorMessage(e)
end
