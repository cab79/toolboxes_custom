%
% Demonstrates how to stream a range of sequential analog inputs using the
% eStream functions. Useful when streaming many analog inputs. AIN channel
% scan list is firstAINChan to firstAINChan + numAINs - 1.
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

    % The number of eStreamRead calls to perform in the stream read loop
    maxRequests = 50;

    numAINs = 8;  % Number of AINs to stream.
    firstAINChan = 0;  % Starting AIN channel. 0 = AIN0.

    % When streaming, negative channels and ranges can be configured
    % for individual analog inputs, but the stream has only one settling
    % time and resolution.

    if getDeviceType(handle) == LJM_CONSTANTS.dtT4
        % T4 configuration

        % Configure the channels to analog input or digital I/O.

        % Update all I/O channels for analog or digital.
        % b1 = Ignored. b0 = Affected.
        dioInhibit = 0;  % b00000000000000000000
        % Set AIN0-AIN3 and AIN firstAINChan to firstAINChan+numAINs-1 as
        % analog inputs (b1), the rest as digital I/O (b0).
        dioAnalogEnable = bitshift((2^numAINs - 1), firstAINChan);
        dioAnalogEnable = bitor(dioAnalogEnable, 15);
        numFrames = 2;
        aNames = NET.createArray('System.String', numFrames);
        aNames(1) = 'DIO_INHIBIT';
        aNames(2) = 'DIO_ANALOG_ENABLE';
        aValues = NET.createArray('System.Double', numFrames);
        aValues(1) = dioInhibit;
        aValues(2) = dioAnalogEnable;
        LabJack.LJM.eWriteNames(handle, 2, aNames, aValues, -1);

        % Configure the analog input ranges
        numFrames = numAINs;
        rangeAINHV = 10.0;  % HV channels range (AIN0-AIN3)
        rangeAINLV = 2.5;  % LV channels range (AIN4+)
        aNames = NET.createArray('System.String', numFrames);
        aValues = NET.createArray('System.Double', numFrames);
        for i = 1:numFrames
            chan = firstAINChan + i - 1;
            aNames(i) = ['AIN' num2str(chan) '_RANGE'];
            if chan < 4
                aValues(i) = rangeAINHV;
            else
                aValues(i) = rangeAINLV;
            end
        end

        % Configure the stream settling times and stream resolution index.
        numFrames = 2;
        aNames = NET.createArray('System.String', numFrames);
        aNames(1) = 'STREAM_SETTLING_US';
        aNames(2) = 'STREAM_RESOLUTION_INDEX';
        aValues = NET.createArray('System.Double', numFrames);
        aValues(1) = 0;  % 0 (default)
        aValues(2) = 0;  % 0 (default)
        LabJack.LJM.eWriteNames(handle, numFrames, aNames, aValues, -1);
    else
        % T7 and other devices configuration

        % Configure the analog input negative channels, ranges, stream
        % settling times and stream resolution index.
        numFrames = 4;
        aNames = NET.createArray('System.String', numFrames);
        aNames(1) = 'AIN_ALL_NEGATIVE_CH';
        aNames(2) = 'AIN_ALL_RANGE';
        aNames(3) = 'STREAM_SETTLING_US';
        aNames(4) = 'STREAM_RESOLUTION_INDEX';
        aValues = NET.createArray('System.Double', numFrames);
        aValues(1) = LJM_CONSTANTS.GND;  % Single-ended
        aValues(2) = 10.0;  % +/-10 V
        aValues(3) = 0;  % 0 (default)
        aValues(4) = 0;  % 0 (default)
        LabJack.LJM.eWriteNames(handle, numFrames, aNames, aValues, -1);
    end

    % Stream Configuration

    % Scan list names to stream
    numAddresses = numAINs;
    aScanListNames = NET.createArray('System.String', numAddresses);
    for i = 1:numAddresses
        chan = firstAINChan + i - 1;
        aScanListNames(i) = ['AIN' num2str(chan)];
    end

    % Scan list addresses to stream
    aScanList = NET.createArray('System.Int32', numAddresses);
    aTypes = NET.createArray('System.Int32', numAddresses);  % Dummy array
    LabJack.LJM.NamesToAddresses(numAddresses, aScanListNames, ...
        aScanList, aTypes);

    scanRate = double(1000);  % Scans per second
    scansPerRead = int32(scanRate/2);
    % Stream reads will be stored in aData. Needs to be at least
    % numAINs*scansPerRead in size.
    aData = NET.createArray('System.Double', numAINs*scansPerRead);

    % Configure and start stream
    numAddresses = aScanList.Length;
    [~, scanRate] = LabJack.LJM.eStreamStart(handle, scansPerRead, ...
        numAddresses, aScanList, scanRate);

    disp(['Stream started with a scan rate of ' ...
          num2str(scanRate) ' Hz.'])

    tic

    try
        disp(['Performing ' num2str(maxRequests) ' stream reads.'])

        totalScans = 0;
        curSkippedSamples = 0;
        totalSkippedSamples = 0;

        for i = 1:maxRequests
            [~, devScanBL, ljmScanBL] = LabJack.LJM.eStreamRead( ...
                handle, aData, 0, 0);

            totalScans = totalScans + scansPerRead;

            % Count the skipped samples which are indicated by -9999
            % values. Skipped samples occur after a device's stream buffer
            % overflows and are reported after auto-recover mode ends.
            % When streaming at faster scan rates in MATLAB, try counting
            % the skipped packets outside your eStreamRead loop if you are
            % getting skipped samples/scan.
            curSkippedSamples = sum(double(aData) == -9999.0);
            totalSkippedSamples = totalSkippedSamples + curSkippedSamples;

            fprintf('eStreamRead %d\n  1st scan out of %d:\n  ', i, ...
                    scansPerRead)
            for j = 1:numAddresses
                fprintf('%s = %.4f, ', char(aScanListNames(j)), aData(j))
                if mod(j, 4) == 0 || j >= numAddresses
                    % Newline after 4 readings or last reading
                    fprintf('\n  ')
                end
            end
            fprintf(['Scans Skipped = %d, Scan Backlogs: Device = %d,' ...
                    ' LJM = %d\n'], (curSkippedSamples/numAddresses), ...
                    devScanBL, ljmScanBL)
        end

        timeElapsed = toc;

        disp(['Total scans = ' num2str(totalScans)])
        disp(['Skipped Scans = ' ...
              num2str(totalSkippedSamples/numAddresses)])
        disp(['Time Taken = ' num2str(timeElapsed) ' seconds'])
        disp(['LJM Scan Rate = ' num2str(scanRate) ' scans/second'])
        disp(['Timed Scan Rate = ' num2str(totalScans/timeElapsed) ...
              ' scans/second'])
        disp(['Sample Rate = ' ...
              num2str(numAddresses*totalScans/timeElapsed) ...
              ' samples/second'])
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
