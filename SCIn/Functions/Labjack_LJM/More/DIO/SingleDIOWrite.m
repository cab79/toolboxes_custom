%
% Demonstrates how to set a single digital output on a LabJack using .NET.
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

    % Setup and call eWriteName to set the DIO state.
    if getDeviceType(handle) == LJM_CONSTANTS.dtT4
        % Setting FIO5 on the LabJack T4. FIO0-FIO3 are reserved for
        % AIN0-AIN3.
        name = 'FIO4';

        % If the FIO/EIO line is an analog input, it needs to first be
        % changed to a digital I/O by reading from the line or setting it
        % to digital I/O with the DIO_ANALOG_ENABLE register.

        % Reading from the digital line in case it was previously an analog
        % input.
        LabJack.LJM.eReadName(handle, name, 0);
    else
        % Setting FIO0 on the LabJack T7 and other devices.
        name = 'FIO0';
    end
    state = 0;  % 0 = output-low, 1 = output-high
    LabJack.LJM.eWriteName(handle, name, state);

    disp(['Set ' name ' state: ' num2str(state)])
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
