%
% Demonstrates how to read a single digital input on a LabJack using .NET.
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

    % Setup and call eReadName to read the DIO state.
    if getDeviceType(handle) == LJM_CONSTANTS.dtT4
        % Reading from FIO5 on the LabJack T4. FIO0-FIO3 are reserved for
        % AIN0-AIN3.
        % Note: Reading a single digital I/O will change the line from
        % analog to digital input.
        name = 'FIO5';
    else
        % Reading from FIO1 on the LabJack T7 and other devices.
        name = 'FIO1';
    end
    state = 0;
    [ljmerror, state] = LabJack.LJM.eReadName(handle, name, state);

    disp([name ' state: ' num2str(state)])
catch e
    showErrorMessage(e)
    LabJack.LJM.CloseAll();
    return
end

try
    % Close handle
    LabJack.LJM.Close(handle);
catch 
    showErrorMessage(e)
end
