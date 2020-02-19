%
% Demonstrates how to read configuration settings on a LabJack using .NET.
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

    % Setup and call eReadNames to read configuration values.
    if getDeviceType(handle) == LJM_CONSTANTS.dtT4
        % LabJack T4 configuration to read
        numFrames = 8;
        aNames = NET.createArray('System.String', numFrames);
        aNames(1) = 'PRODUCT_ID';
        aNames(2) = 'HARDWARE_VERSION';
        aNames(3) = 'FIRMWARE_VERSION';
        aNames(4) = 'BOOTLOADER_VERSION';
        aNames(5) = 'SERIAL_NUMBER';
        aNames(6) = 'POWER_ETHERNET_DEFAULT';
        aNames(7) = 'POWER_AIN_DEFAULT';
        aNames(8) = 'POWER_LED_DEFAULT';
    else
        % LabJack T7 and other devices configuration to read
        numFrames = 10;
        aNames = NET.createArray('System.String', numFrames);
        aNames(1) = 'PRODUCT_ID';
        aNames(2) = 'HARDWARE_VERSION';
        aNames(3) = 'FIRMWARE_VERSION';
        aNames(4) = 'BOOTLOADER_VERSION';
        aNames(5) = 'WIFI_VERSION';
        aNames(6) = 'SERIAL_NUMBER';
        aNames(7) = 'POWER_ETHERNET_DEFAULT';
        aNames(8) = 'POWER_WIFI_DEFAULT';
        aNames(9) = 'POWER_AIN_DEFAULT';
        aNames(10) = 'POWER_LED_DEFAULT';
    end
    aValues = NET.createArray('System.Double', numFrames);
    LabJack.LJM.eReadNames(handle, numFrames, aNames, aValues, 0);

    disp('Configuration settings:')
    for i = 1:numFrames
        disp([' ' char(aNames(i)) ': ' num2str(aValues(i))])
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
