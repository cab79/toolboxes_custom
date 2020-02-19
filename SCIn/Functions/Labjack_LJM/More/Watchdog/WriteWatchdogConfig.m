%
% Demonstrates how to configure the Watchdog on a LabJack using .NET.
%
% support@labjack.com
%

clc  % Clear the MATLAB command window
clear  % Clear the MATLAB variables

%Make the LJM .NET assembly visible in MATLAB
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

    % Setup and call eWriteNames to configure the Watchdog. Disable the
    % Watchdog first before any other configuration.
    numFrames = 16;
    aNames = NET.createArray('System.String', numFrames);
    aNames(1) = 'WATCHDOG_ENABLE_DEFAULT';
    aNames(2) = 'WATCHDOG_ADVANCED_DEFAULT';
    aNames(3) = 'WATCHDOG_TIMEOUT_S_DEFAULT';
    aNames(4) = 'WATCHDOG_STARTUP_DELAY_S_DEFAULT';
    aNames(5) = 'WATCHDOG_STRICT_ENABLE_DEFAULT';
    aNames(6) = 'WATCHDOG_STRICT_KEY_DEFAULT';
    aNames(7) = 'WATCHDOG_RESET_ENABLE_DEFAULT';
    aNames(8) = 'WATCHDOG_DIO_ENABLE_DEFAULT';
    aNames(9) = 'WATCHDOG_DIO_STATE_DEFAULT';
    aNames(10) = 'WATCHDOG_DIO_DIRECTION_DEFAULT';
    aNames(11) = 'WATCHDOG_DIO_INHIBIT_DEFAULT';
    aNames(12) = 'WATCHDOG_DAC0_ENABLE_DEFAULT';
    aNames(13) = 'WATCHDOG_DAC0_DEFAULT';
    aNames(14) = 'WATCHDOG_DAC1_ENABLE_DEFAULT';
    aNames(15) = 'WATCHDOG_DAC1_DEFAULT';
    aNames(16) = 'WATCHDOG_ENABLE_DEFAULT';
    aValues = NET.createArray('System.Double', numFrames);
    % Disable (0) WATCHDOG_ENABLE_DEFAULT before configuring.
    aValues(1) = 0;
    aValues(2) = 0;
    % Set WATCHDOG_TIMEOUT_S_DEFAULT to 20 seconds.
    aValues(3) = 20;
    aValues(4) = 0;
    aValues(5) = 0;
    aValues(6) = 0;
    % Enable (1) WATCHDOG_RESET_ENABLE_DEFAULT.
    aValues(7) = 1;
    aValues(8) = 0; 
    aValues(9) = 0;
    aValues(10) = 0;
    aValues(11) = 0;
    aValues(12) = 0;
    aValues(13) = 0;
    aValues(14) = 0;
    aValues(15) = 0;
    % Change to 1 to enable WATCHDOG_ENABLE_DEFAULT.
    aValues(16) = 0;
    LabJack.LJM.eWriteNames(handle, numFrames, aNames, aValues, 0);

    disp('Set Watchdog configuration:')
    for i = 1:numFrames
        disp(['    ' char(aNames(i)) ' : ' num2str(aValues(i))])
    end

    % Close handle
    LabJack.LJM.Close(handle);
catch e
    showErrorMessage(e)
    LabJack.LJM.CloseAll();
end
