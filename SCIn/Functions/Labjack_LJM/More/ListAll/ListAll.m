%
% Demonstrates usage of the ListAll functions (LJM_ListAll) which scans for
% LabJack devices and returns information about the found devices. This
% will only find LabJack devices supported by the LJM library.
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

try
    % Device constant to string map
    devMap = containers.Map( ...
        {LJM_CONSTANTS.dtT7, LJM_CONSTANTS.dtT4, ...
         LJM_CONSTANTS.dtDIGIT}, ...
        {'T7', 'T4', 'Digit'});
    % Connection constant to string map
    connMap = containers.Map( ...
        {LJM_CONSTANTS.ctUSB, LJM_CONSTANTS.ctTCP, ...
         LJM_CONSTANTS.ctETHERNET, LJM_CONSTANTS.ctWIFI}, ...
        {'USB', 'TCP', 'Ethernet', 'WiFi'});

    MAX_LIST_SIZE = LJM_CONSTANTS.LIST_ALL_SIZE;
    numFound = 0;
    aDeviceTypes = NET.createArray('System.Int32', MAX_LIST_SIZE);
    aDeviceTypes = NET.createArray('System.Int32', MAX_LIST_SIZE);
    aConnectionTypes = NET.createArray('System.Int32', MAX_LIST_SIZE);
    aSerialNumbers = NET.createArray('System.Int32', MAX_LIST_SIZE);
    aIPAddresses = NET.createArray('System.Int32', MAX_LIST_SIZE);

    % Find and display LabJack devices with ListAllS.
    [~, numFound] = LabJack.LJM.ListAllS('ANY', 'ANY', 0, aDeviceTypes, ...
        aConnectionTypes, aSerialNumbers, aIPAddresses);
    disp(['ListAllS found ' num2str(numFound) ' LabJacks:'])

    % Find and display LabJack devices with ListAll.
    % [~, numFound] = LabJack.LJM.ListAll(LJM_CONSTANTS.dtANY, ...
    %     LJM_CONSTANTS.ctANY, 0, aDeviceTypes, aConnectionTypes, ...
    %     aSerialNumbers, aIPAddresses);
    % disp(['ListAll found ' num2str(numFound) ' LabJacks:'])

    fprintf('%-18s%-18s%-18s%-18s\n', 'Connection Type', 'Device Type', ...
            'Serial Number', 'IP Address')
    for i = 1:numFound
        if isKey(devMap, aDeviceTypes(i))
            devStr = devMap(aDeviceTypes(i));
        else
            devStr = num2str(aDeviceTypes(i));
        end
        if isKey(connMap, aConnectionTypes(i))
            connStr = connMap(aConnectionTypes(i));
        else
            connStr = num2str(aConnectionTypes(i));
        end
        [~, ipAddrStr] = LabJack.LJM.NumberToIP(aIPAddresses(i), '');
        fprintf('%-18s%-18s%-18d%-18s\n', devStr, connStr, ...
                aSerialNumbers, char(ipAddrStr));
    end
catch e
    showErrorMessage(e)
end
