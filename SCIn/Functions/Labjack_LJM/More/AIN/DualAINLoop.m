%
% Demonstrates reading 2 analog inputs (AINs) in a loop using .NET.
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

    % Setup and call eWriteNames to configure AINs.
    if getDeviceType(handle) == LJM_CONSTANTS.dtT4
        % LabJack T4 configuration

        % AIN0 and AIN1:
        %     Range = +/-10 V. Only AIN0-AIN3 can support +/-10 V range.
        %     Resolution index = 0 (default)
        %     Settling = 0 (auto)
        numFrames = 6;
        aNames = NET.createArray('System.String', numFrames);
        aNames(1) = 'AIN0_RANGE';
        aNames(2) = 'AIN0_RESOLUTION_INDEX';
        aNames(3) = 'AIN0_SETTLING_US';
        aNames(4) = 'AIN1_RANGE';
        aNames(5) = 'AIN1_RESOLUTION_INDEX';
        aNames(6) = 'AIN1_SETTLING_US';
        aValues = NET.createArray('System.Double', numFrames);
        aValues(1) = 10;
        aValues(2) = 0;
        aValues(3) = 0;
        aValues(4) = 10;
        aValues(5) = 0;
        aValues(6) = 0;
    else
        % LabJack T7 and other devices configuration

        % AIN0 and AIN1:
        %     Negative Channel = 199 (Single-ended)
        %     Range = +/-10 V
        %     Resolution index = 0 (default)
        %     Settling = 0 (auto)
        numFrames = 8;
        aNames = NET.createArray('System.String', numFrames);
        aNames(1) = 'AIN0_NEGATIVE_CH';
        aNames(2) = 'AIN0_RANGE';
        aNames(3) = 'AIN0_RESOLUTION_INDEX';
        aNames(4) = 'AIN0_SETTLING_US';
        aNames(5) = 'AIN1_NEGATIVE_CH';
        aNames(6) = 'AIN1_RANGE';
        aNames(7) = 'AIN1_RESOLUTION_INDEX';
        aNames(8) = 'AIN1_SETTLING_US';
        aValues = NET.createArray('System.Double', numFrames);
        aValues(1) = 199;
        aValues(2) = 10;
        aValues(3) = 0;
        aValues(4) = 0;
        aValues(5) = 199;
        aValues(6) = 10;
        aValues(7) = 0;
        aValues(8) = 0;
    end
    LabJack.LJM.eWriteNames(handle, numFrames, aNames, aValues, 0);

    disp('Set configuration:');
    for i = 1:numFrames
        disp(['  ' char(aNames(i)) ': ' num2str(aValues(i))])
    end

    % Setup and call eReadNames to read AINs.
    numFrames = 2;
    aNames = NET.createArray('System.String', numFrames);
    aNames(1) = 'AIN0';
    aNames(2) = 'AIN1';
    aValues = NET.createArray('System.Double', numFrames);

    numReadings = 10;
    delay = 1;  % Delay (in sec.) between readings
    disp(['Performing ' num2str(numReadings) ' AIN0 and AIN1 readings ' ...
          'with ' num2str(delay) ' second delay between readings:']);

    for i = 1:numReadings
        LabJack.LJM.eReadNames(handle, numFrames, aNames, aValues, 0);
        disp(['  ' char(aNames(1)) ': ' num2str(aValues(1)) ' V, ' ...
              char(aNames(2)) ': ' num2str(aValues(2)) ' V'])
        pause(delay);
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
