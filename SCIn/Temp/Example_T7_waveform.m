% T7 CONTROL EXAMPLE
% https://labjack.com/support/datasheets/t-series/communication/stream-mode/stream-out
% https://labjack.com/support/datasheets/t-series/communication/stream-mode/stream-out/stream-out-description

clc  % Clear the MATLAB command window
clear  % Clear the MATLAB variables

addpath(genpath('Q:\MATLAB\toolboxes_external\Labjack_LJM'))

%Make the LJM .NET assembly visible in MATLAB
ljmAsm = NET.addAssembly('LabJack.LJM');

% Creating an object to nested class LabJack.LJM.CONSTANTS
t = ljmAsm.AssemblyHandle.GetType('LabJack.LJM+CONSTANTS');
LJM_CONSTANTS = System.Activator.CreateInstance(t);


%Precision needed is 0.005s (5ms) so aim for 0.001 as standard. 
% For a 1s train this is 2048, since each digital output is 2 bytes. 
% Max train length is therefore 8s at this precision.
%The standard DS8R can be triggered at a maximum frequency of 1000Hz


try
    %Open first found LabJack
    [ljmError, handle] = LabJack.LJM.OpenS('ANY', 'ANY', 'ANY', 0);
    
    showDeviceInfo(handle);
    
    %Call eReadNames to read the serial number from the LabJack.
    name = 'SERIAL_NUMBER';
    [ljmError, value] = LabJack.LJM.eReadName(handle, name, 0);
    
    disp('eReadName result:')
    disp(['  ' name ' = ' num2str(value)])
catch e
    showErrorMessage(e)
end

%Call eWriteNames
try
    % Setup stream-out
    numAddressesOut = 1;
    aNamesOut = NET.createArray('System.String', numAddressesOut);
    aNamesOut(1) = 'FIO_STATE';
    aAddressesOut = NET.createArray('System.Int32', numAddressesOut);
    aTypesOut = NET.createArray('System.Int32', numAddressesOut);  % Dummy
    LabJack.LJM.NamesToAddresses(numAddressesOut, aNamesOut, ...
        aAddressesOut, aTypesOut);
    
    % Allocate memory for the stream-out buffer
    buffer_size = 2048;
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_ENABLE', 0);
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_TARGET', aAddressesOut(1));
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_SIZE', buffer_size);
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_ENABLE', 1);
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_LOOP_NUM_VALUES', 0);
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_SET_LOOP', 1);
    
catch e
    showErrorMessage(e)
    LabJack.LJM.CloseAll();
    return
end

% Load the waveform data points
digital_waveform = repmat([0 5],1,10);
LJMError = LabJack.LJM.eWriteNameArray(handle, 'STREAM_OUT0_BUFFER_U16', length(digital_waveform), digital_waveform, 0);


% Setup stream-out
numAddressesOut = 1;
aNamesOut = NET.createArray('System.String', numAddressesOut);
aNamesOut(1) = 'STREAM_OUT0';
aAddressesOut = NET.createArray('System.Int32', numAddressesOut);
aTypesOut = NET.createArray('System.Int32', numAddressesOut);  % Dummy
    LabJack.LJM.NamesToAddresses(numAddressesOut, aNamesOut, ...
        aAddressesOut, aTypesOut);

% Configure and start stream
aScanList = aAddressesOut(1); % FIO0
scanRate = 20; %1000; % Hz
scansPerRead = scanRate / 2; % eStreamRead frequency = ScanRate / ScansPerRead
numAddresses = length(aScanList);

[~, scanRate] = LabJack.LJM.eStreamStart(handle, scansPerRead, ...
    numAddresses, aScanList, scanRate);
disp(['Stream started with a scan rate of ' num2str(scanRate) ...
              ' Hz.'])

pause(length(digital_waveform)/scanRate)

% stop stream
disp('Stop Stream')
LabJack.LJM.eStreamStop(handle);

% Close handle
LabJack.LJM.Close(handle);
