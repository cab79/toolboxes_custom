%
% Demonstrates how to read the WiFi MAC from a LabJack using .NET.
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

    % Call eReadAddressByteArray to read the WiFi MAC (address 60024).
    % We are reading a byte array which is the big endian binary
    % representation of the 64-bit MAC.
    aBytes = NET.createArray('System.Byte', 8);
    LabJack.LJM.eReadAddressByteArray(handle, 60024, 8, aBytes, -1);

    % Convert returned .NET bytes to MATLAB bytes
    macBytes = uint8(aBytes);

    % Convert big endian byte array to a 64-bit signed integer value
    [computerType, maxSize, endian] = computer;
    if endian == 'L'
        macBytes = macBytes(end:-1:1);
    end
    macNumber = typecast(macBytes, 'int64');

    % Convert the MAC value/number to its string representation
    macString = '';
    [ljmError, macString] = LabJack.LJM.NumberToMAC(macNumber, macString);

    disp(['WiFi MAC : ' num2str(macNumber) ' - ' char(macString)])

    % Close handle
    LabJack.LJM.Close(handle);
catch e
    showErrorMessage(e)
    LabJack.LJM.CloseAll();
end
