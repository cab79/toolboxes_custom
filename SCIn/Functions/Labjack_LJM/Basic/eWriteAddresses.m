%
% Demonstrates how to use the eWriteAddresses (LJM_eWriteAddresses)
% function using .NET.
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

    % Setup and call eWriteAddresses to write values.
    numFrames = 2;
    aAddresses = NET.createArray('System.Int32', numFrames);
    aAddresses(1) = 1000;  % DAC0
    aAddresses(2) = 55110;  % TEST_UINT16
    aTypes = NET.createArray('System.Int32', numFrames);
    aTypes(1) = LJM_CONSTANTS.FLOAT32;
    aTypes(2) = LJM_CONSTANTS.UINT16;
    aValues = NET.createArray('System.Double', numFrames);
    aValues(1) = 2.5;  % 2.5 V
    aValues(2) = 12345;
    LabJack.LJM.eWriteAddresses(handle, numFrames, aAddresses, aTypes, ...
        aValues, 0);

    disp('eWriteAddresses:')
    for i=1:numFrames,
        disp(['  Address: ' num2str(aAddresses(i)) ', Data Type: ' ...
              num2str(aTypes(i)) ', Value: ' num2str(aValues(i))])
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
