%
% Demonstrates how to read the WiFi configuration from a LabJack using
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

    % Any device, Any connection, Any identifier
    % [ljmError, handle] = LabJack.LJM.Open(LJM_CONSTANTS.dtANY, ...
    %     LJM_CONSTANTS.ctANY, 'ANY', handle);

    showDeviceInfo(handle);

    if getDeviceType(handle) == LJM_CONSTANTS.dtT4
        disp('The LabJack T4 does not support WiFi.')
        LabJack.LJM.Close(handle);
        return
    end

    % Setup and call eReadNames to read WiFi configuration.
    numFrames = 9;
    aNames = NET.createArray('System.String', numFrames);
    aNames(1) = 'WIFI_IP';
    aNames(2) = 'WIFI_SUBNET';
    aNames(3) = 'WIFI_GATEWAY';
    aNames(4) = 'WIFI_DHCP_ENABLE';
    aNames(5) = 'WIFI_IP_DEFAULT';
    aNames(6) = 'WIFI_SUBNET_DEFAULT';
    aNames(7) = 'WIFI_GATEWAY_DEFAULT';
    aNames(8) = 'WIFI_DHCP_ENABLE_DEFAULT';
    aNames(9) = 'WIFI_STATUS';
    aValues = NET.createArray('System.Double', numFrames);
    LabJack.LJM.eReadNames(handle, numFrames, aNames, aValues, 0);

    disp('WiFi configuration:')
    str = '';
    for i = 1:numFrames
        k1 = strfind(char(aNames(i)), 'WIFI_STATUS');
        k2 = strfind(char(aNames(i)), 'WIFI_DHCP_ENABLE');
        if ~isempty(k1) || ~isempty(k2)
            disp(['    ' char(aNames(i)) ' : ' num2str(aValues(i))])
        else
            [ljmError, str] = LabJack.LJM.NumberToIP( ...
                typecast(uint32(aValues(i)), 'int32'), str);
            disp(['    ' char(aNames(i)) ' : ' num2str(aValues(i)) ...
                  ' - ' char(str)])
        end
    end

    % Setup and call eReadNameString to read the WiFi SSID string.
    name = 'WIFI_SSID';
    str = '';
    [ljmError, str] = LabJack.LJM.eReadNameString(handle, name, str);

    disp(['    ' name ' : ' char(str)])

    % Close handle
    LabJack.LJM.Close(handle);
catch e
    showErrorMessage(e)
    LabJack.LJM.CloseAll();
end
