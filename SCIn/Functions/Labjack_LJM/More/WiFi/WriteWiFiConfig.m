%
% Demonstrates how to configure the WiFi settings on a LabJack using .NET.
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

    % Setup and call eWriteNames to configure WiFi default settings.
    numFrames = 3;
    aNames = NET.createArray('System.String', numFrames);
    aNames(1) = 'WIFI_IP_DEFAULT';
    aNames(2) = 'WIFI_SUBNET_DEFAULT';
    aNames(3) = 'WIFI_GATEWAY_DEFAULT';
    aValues = NET.createArray('System.Double', numFrames);
    [ljmError, aValues(1)] = LabJack.LJM.IPToNumber('192.168.1.207', 0);
    aValues(1) = typecast(int32(aValues(1)), 'uint32');
    [ljmError, aValues(2)] = LabJack.LJM.IPToNumber('255.255.255.0', 0);
    aValues(2) = typecast(int32(aValues(2)), 'uint32');
    [ljmError, aValues(3)] = LabJack.LJM.IPToNumber('192.168.1.1', 0);
    aValues(3) = typecast(int32(aValues(1)), 'uint32');
    LabJack.LJM.eWriteNames(handle, numFrames, aNames, aValues, 0);

    disp('Set WiFi configuration:')
    str = '';
    for i = 1:numFrames
        [ljmError, str] = LabJack.LJM.NumberToIP( ...
            typecast(uint32(aValues(i)), 'int32'), str);
        disp(['    ' char(aNames(i)) ' : ' num2str(aValues(i)) ...
              ' - ' char(str)])
    end

    % Setup and call eWriteString to configure the default WiFi SSID.
    name = 'WIFI_SSID_DEFAULT';
    str = 'LJOpen';
    LabJack.LJM.eWriteNameString(handle, name, str);
    disp(['    ' name ' : ' char(str)])

    % Setup and call eWriteString to configure the default WiFi password.
    name = 'WIFI_PASSWORD_DEFAULT';
    str = 'none';
    LabJack.LJM.eWriteNameString(handle, name, str);
    disp(['    ' name ' : ' char(str)])

    % Setup and call eWriteName to apply the new WiFi configuration.
    name = 'WIFI_APPLY_SETTINGS';
    value = 1;  % 1 = apply
    LabJack.LJM.eWriteName(handle, name, value);
    disp(['    ' name ' : ' num2str(value)])

    % Close handle
    LabJack.LJM.Close(handle);
catch e
    showErrorMessage(e)
    LabJack.LJM.CloseAll();
end
