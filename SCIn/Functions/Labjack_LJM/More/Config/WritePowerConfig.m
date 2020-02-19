%
% Demonstrates how to configure default power settings on a LabJack .NET.
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

    % Setup and call eWriteNames to write configuration values.
    numFrames = 4;
    aNames = NET.createArray('System.String', numFrames);
    aNames(1) = 'POWER_ETHERNET_DEFAULT';
    aNames(2) = 'POWER_WIFI_DEFAULT';
    aNames(3) = 'POWER_AIN_DEFAULT';
    aNames(4) = 'POWER_LED_DEFAULT';
    aValues = NET.createArray('System.Double', numFrames);
    aValues(1) = 1;  % Ethernet on
    aValues(2) = 0;  % WiFi off
    aValues(3) = 1;  % AIN on
    aValues(4) = 1;  % LED on
    LabJack.LJM.eWriteNames(handle, numFrames, aNames, aValues, 0);

    disp('Set configuration settings:');
    for i = 1:numFrames
        disp(['  ' char(aNames(i)) ' : ' num2str(aValues(i))])
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
