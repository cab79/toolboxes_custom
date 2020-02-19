%
% Demonstrates how to read the ethernet configuration settings from a
% LabJack using .NET.
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

    %Setup and call eReadNames to read ethernet configuration.
    numFrames = 8;
    aNames = NET.createArray('System.String', numFrames);
    aNames(1) = 'ETHERNET_IP';
    aNames(2) = 'ETHERNET_SUBNET';
    aNames(3) = 'ETHERNET_GATEWAY';
    aNames(4) = 'ETHERNET_IP_DEFAULT';
    aNames(5) = 'ETHERNET_SUBNET_DEFAULT';
    aNames(6) = 'ETHERNET_GATEWAY_DEFAULT';
    aNames(7) = 'ETHERNET_DHCP_ENABLE';
    aNames(8) = 'ETHERNET_DHCP_ENABLE_DEFAULT';
    aValues = NET.createArray('System.Double', numFrames);
    LabJack.LJM.eReadNames(handle, numFrames, aNames, aValues, -1);

    disp('Ethernet configuration:')
    str = '';
    for i = 1:numFrames
        k = strfind(char(aNames(i)), 'ETHERNET_DHCP_ENABLE');
        if isempty(k)
            number = typecast(uint32(aValues(i)), 'int32');
            [ljmError, str] = LabJack.LJM.NumberToIP(number, str);
            disp(['    ' char(aNames(i)) ' : ' num2str(aValues(i)) ...
                  ' - ' char(str)])
        else
            disp(['    ' char(aNames(i)) ' : ' num2str(aValues(i))])
        end
    end

    % Close handle
    LabJack.LJM.Close(handle);
catch e
    showErrorMessage(e)
    LabJack.LJM.CloseAll();
end
