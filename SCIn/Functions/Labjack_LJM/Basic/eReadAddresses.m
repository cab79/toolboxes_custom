%
% Demonstrates how to use the eReadAddresses (LJM_eReadAddresses) function using
% .NET.
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

    % Setup and call eReadAddresses to read values.
    numFrames = 3;
    aAddresses = NET.createArray('System.Int32', numFrames);
    aAddresses(1) = 60028;  % Serial number
    aAddresses(2) = 60000;  % Product ID
    aAddresses(3) = 60004;  % Firmware version
    aTypes = NET.createArray('System.Int32', numFrames);
    aTypes(1) = LJM_CONSTANTS.UINT32;
    aTypes(2) = LJM_CONSTANTS.FLOAT32;
    aTypes(3) = LJM_CONSTANTS.FLOAT32;
    aValues = NET.createArray('System.Double', numFrames);
    LabJack.LJM.eReadAddresses(handle, numFrames, aAddresses, aTypes, ...
        aValues, 0);

    disp('eReadAddresses results:')
    for i = 1:numFrames
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
