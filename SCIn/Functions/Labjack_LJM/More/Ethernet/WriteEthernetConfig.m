%
% Demonstrates how to set ethernet configuration settings on a LabJack
% using .NET.
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

    % Setup and call eWriteNames to set the ethernet configuration.
    numFrames = 4;
    aNames = NET.createArray('System.String', numFrames);
    aNames(1) = 'ETHERNET_IP_DEFAULT';
    aNames(2) = 'ETHERNET_SUBNET_DEFAULT';
    aNames(3) = 'ETHERNET_GATEWAY_DEFAULT';
    aNames(4) = 'ETHERNET_DHCP_ENABLE_DEFAULT';
    aValues = NET.createArray('System.Double', numFrames);
    [ljmError, aValues(1)] = LabJack.LJM.IPToNumber('192.168.1.207', 0);
    aValues(1) = typecast(int32(aValues(1)), 'uint32');
    [ljmError, aValues(2)] = LabJack.LJM.IPToNumber('255.255.255.0', 0);
    aValues(2) = typecast(int32(aValues(2)), 'uint32');
    [ljmError, aValues(3)] = LabJack.LJM.IPToNumber('192.168.1.1', 0);
    aValues(3) = typecast(int32(aValues(3)), 'uint32');
    aValues(4) = 1;
    LabJack.LJM.eWriteNames(handle, numFrames, aNames, aValues, 0);

    disp('Set ethernet configuration:')
    str = '';
    for i = 1:numFrames
        k = strfind(char(aNames(i)), 'ETHERNET_DHCP_ENABLE_DEFAULT');
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
