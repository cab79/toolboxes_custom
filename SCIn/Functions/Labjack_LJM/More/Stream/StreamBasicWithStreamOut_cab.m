%
% Demonstrates setting up stream-in and stream-out together, then reading
% stream-in values using .NET.
%
% Connect a wire from AIN0 to DAC0 to see the effect of stream-out on
% stream-in channel 0.
%
% support@labjack.com
%

clc  % Clear the MATLAB command window
clear  % Clear the MATLAB variables

addpath(genpath('Q:\MATLAB\toolboxes_external\Labjack_LJM'))

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

    % Setup stream-out
    numAddressesOut = 1;
    aNamesOut = NET.createArray('System.String', numAddressesOut);
    aNamesOut(1) = 'DAC0';
    aAddressesOut = NET.createArray('System.Int32', numAddressesOut);
    aTypesOut = NET.createArray('System.Int32', numAddressesOut);  % Dummy
    LabJack.LJM.NamesToAddresses(numAddressesOut, aNamesOut, ...
        aAddressesOut, aTypesOut);
    
    buffer_size = 32;%max(32,pow2(floor(log2(16)))*2);

    % Allocate memory for the stream-out buffer
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_TARGET', aAddressesOut(1));
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_SIZE', buffer_size);
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_ENABLE', 1);

    % Write values to the stream-out buffer
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_LOOP_SIZE', 0);
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 5.0);  % 0.0 V
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 5.0);  % 0.0 V
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 0);  % 1.0 V
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 0);  % 1.0 V
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 5.0);  % 2.0 V
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 5.0);  % 2.0 V
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 0);  % 3.0 V
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 0);  % 3.0 V
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 5.0);  % 4.0 V
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 5.0);  % 4.0 V
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 0);  % 5.0 V
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 0);  % 5.0 V
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 5.0);  % 4.0 V
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 5.0);  % 4.0 V
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 0);  % 5.0 V
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 0);  % 5.0 V

    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_SET_LOOP', 1);

%     [~, value] = LabJack.LJM.eReadName(handle, ...
%         'STREAM_OUT0_BUFFER_STATUS', 0);
%     disp(['STREAM_OUT0_BUFFER_STATUS = ' num2str(value)])

    % Stream  configuration

    % Scan list names to stream-in
%     numAddressesIn = 2;
%     aScanListNames = NET.createArray('System.String', numAddressesIn);
%     aScanListNames(1) = 'AIN0';
%     aScanListNames(2) = 'AIN1';

    % Scan list addresses to stream
%     aScanList = NET.createArray('System.Int32', ...
%         (numAddressesIn + numAddressesOut));

    % Get stream-in addresses
%     aTypes = NET.createArray('System.Int32', numAddressesIn); % Dummy
%     LabJack.LJM.NamesToAddresses(numAddressesIn, aScanListNames, ...
%         aScanList, aTypes);
    
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

    % Add the scan list outputs to the end of the scan list.
    % STREAM_OUT0 = 4800, STREAM_OUT1 = 4801, ...
%     aScanList(numAddressesIn+1) = 4800;  % STREAM_OUT0
    % If we had more STREAM_OUTs
    % aScanList[numAddressesIn+2] = 4801;  % STREAM_OUT1
    % etc.

    scanRate = double(2);  % Scans per second
    scansPerRead = scanRate/2;  % Scans returned by eStreamRead call

    % Stream reads will be stored in aData. Needs to be at least
    % numAddresses*scansPerRead in size.
%     aData = NET.createArray('System.Double', numAddressesIn*scansPerRead);

    try
        % When streaming, negative channels and ranges can be configured for
        % individual analog inputs, but the stream has only one settling time
        % and resolution.

        if getDeviceType(handle) == LJM_CONSTANTS.dtT4
            % LabJack T4 configuration

            % AIN0 and AIN1 ranges are +/-10 V, stream settling is 0 (default) and
            % stream resolution index is 0 (default).
            numFrames = 4;
            aNames = NET.createArray('System.String', numFrames);
            aNames(1) = 'AIN0_RANGE';
            aNames(2) = 'AIN1_RANGE';
            aNames(3) = 'STREAM_SETTLING_US';
            aNames(4) = 'STREAM_RESOLUTION_INDEX';
            aValues = NET.createArray('System.Double', numFrames);
            aValues(1) = 10.0;
            aValues(2) = 10.0;
            aValues(3) = 0;
            aValues(4) = 0;
        else
            % LabJack T7 and other devices configuration

%             % Ensure triggered stream is disabled.
%             LabJack.LJM.eWriteName(handle, 'STREAM_TRIGGER_INDEX', 0);
% 
%             % Enabling internally-clocked stream.
%             LabJack.LJM.eWriteName(handle, 'STREAM_CLOCK_SOURCE', 0);

            % All negative channels are single-ended, AIN0 and AIN1 ranges are
            % +/-10 V, stream settling is 0 (default) and stream resolution index
            % is 0 (default).
%             numFrames = 5;
%             aNames = NET.createArray('System.String', numFrames);
%             aNames(1) = 'AIN_ALL_NEGATIVE_CH';
%             aNames(2) = 'AIN0_RANGE';
%             aNames(3) = 'AIN1_RANGE';
%             aNames(4) = 'STREAM_SETTLING_US';
%             aNames(5) = 'STREAM_RESOLUTION_INDEX';
%             aValues = NET.createArray('System.Double', numFrames);
%             aValues(1) = LJM_CONSTANTS.GND;
%             aValues(2) = 10.0;
%             aValues(3) = 10.0;
%             aValues(4) = 0;
%             aValues(5) = 0;
        end
        % Write the analog inputs' negative channels (when applicable), ranges
        % stream settling time and stream resolution configuration.
%         LabJack.LJM.eWriteNames(handle, numFrames, aNames, aValues, -1);

        % Configure and start stream
        numAddresses = 1;
        [~, scanRate] = LabJack.LJM.eStreamStart(handle, scansPerRead, ...
            numAddresses, aScanList, scanRate);

        disp(['Stream started with a scan rate of ' num2str(scanRate) ...
              ' Hz.'])
          pause(10)
% 
%         % The number of eStreamRead calls to perform in the stream read loop
%         maxRequests = 20;
% 
%         % Make a cell array out of scan list names. Helps performance when
%         % converting and displaying in the loop.
%         aScanListNamesML = cell(aScanListNames);
% 
%         tic
% 
%         disp(['Performing ' num2str(maxRequests) ' stream reads.'])
% 
%         totalScans = 0;
%         curSkippedSamples = 0;
%         totalSkippedSamples = 0;
% 
%         for i = 1:maxRequests
%             [~, devScanBL, ljmScanBL] = LabJack.LJM.eStreamRead( ...
%                 handle, aData, 0, 0);
% 
%             totalScans = totalScans + scansPerRead;
% 
%             % Count the skipped samples which are indicated by -9999
%             % values. Skipped samples occur after a device's stream buffer
%             % overflows and are reported after auto-recover mode ends.
%             % When streaming at faster scan rates in MATLAB, try counting
%             % the skipped packets outside your eStreamRead loop if you are
%             % getting skipped samples/scan.
%             curSkippedSamples = sum(double(aData) == -9999.0);
%             totalSkippedSamples = totalSkippedSamples + curSkippedSamples;
% 
%             disp(['eStreamRead ' num2str(i)])
%             slIndex = 1;
%             for j = 1:scansPerRead*numAddressesIn
%                 fprintf('  %s = %.4f,', ...
%                         char(aScanListNamesML(slIndex)), aData(j))
%                 slIndex = slIndex + 1;
%                 if slIndex > numAddressesIn
%                     slIndex = 1;
%                     fprintf('\n')
%                 end
%             end
%             disp(['  Scans Skipped = ' ...
%                   num2str(curSkippedSamples/numAddressesIn) ...
%                   ', Scan Backlogs: Device = ' num2str(devScanBL) ...
%                   ', LJM = ' num2str(ljmScanBL)])
%         end
% 
%         timeElapsed = toc;
% 
%         disp(['Total scans = ' num2str(totalScans)])
%         disp(['Skipped Scans = ' ...
%               num2str(totalSkippedSamples/numAddressesIn)])
%         disp(['Time Taken = ' num2str(timeElapsed) ' seconds'])
%         disp(['LJM Scan Rate = ' num2str(scanRate) ' scans/second'])
%         disp(['Timed Scan Rate = ' num2str(totalScans/timeElapsed) ...
%               ' scans/second'])
%         disp(['Sample Rate = ' ...
%               num2str(numAddressesIn*totalScans/timeElapsed) ...
%               ' samples/second'])
    catch e
        showErrorMessage(e)
    end

    disp('Stop Stream')
    LabJack.LJM.eStreamStop(handle);

    % Close handle
    LabJack.LJM.Close(handle);
catch e
    showErrorMessage(e)
    LabJack.LJM.CloseAll();
end
