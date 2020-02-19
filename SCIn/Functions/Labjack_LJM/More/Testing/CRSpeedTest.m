%
% Performs LabJack operations in a loop and reports the timing statistics
% for the operations.
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

    devType = getDeviceType(handle);

    % Number of iterations to perform in the loop
    numIterations = 1000;

    % Analog input settings
    numAIN = 1;  % Number of analog inputs to read
    rangeAIN = 10.0;  % T7 AIN range
    rangeAINHV = 10.0;  % T4 HV range
    rangeAINLV = 2.5;  % T4 LV range
    resolutionAIN = 1.0;

    % Digital settings
    readDigital = false;
    writeDigital = false;

    % Analog output settings
    writeDACs = false;

    % Use eAddresses (true) or eNames (false) in the operations loop.
    useAddresses = true;

    % Variables for LJM library calls
    numFrames = 0;

    if devType == LJM_CONSTANTS.dtT4
        % For the T4, configure the channels to analog input or digital
        % I/O.

        % Update all digital I/O channels.
        % b1 = Ignored. b0 = Affected.
        dioInhibit = 0;  % b00000000000000000000
        % Set AIN 0 to numAIN-1 as analog inputs (b1), the rest as digital
        % I/O (b0).
        dioAnalogEnable = 2^numAIN - 1;
        numFrames = 2;
        aNames = NET.createArray('System.String', numFrames);
        aNames(1) = 'DIO_INHIBIT';
        aNames(2) = 'DIO_ANALOG_ENABLE';
        aValues = NET.createArray('System.Double', numFrames);
        aValues(1) = dioInhibit;
        aValues(2) = dioAnalogEnable;
        LabJack.LJM.eWriteNames(handle, 2, aNames, aValues, -1);
        if writeDigital
            % Update only digital I/O channels in future digital write
            % calls. b1 = Ignored. b0 = Affected.
            dioInhibit = dioAnalogEnable;
            LabJack.LJM.eWriteName(handle, 'DIO_INHIBIT', dioInhibit);
        end
    end

    if numAIN > 0
        % Configure analog input settings
        numFrames = numAIN*2;
        aNames = NET.createArray('System.String', numFrames);
        aValues = NET.createArray('System.Double', numFrames);
        for i = 0:numAIN-1
            aNames(i*2 + 1) = ['AIN' num2str(i) '_RANGE'];
            if devType == LJM_CONSTANTS.dtT4
                % T4 range
                if i < 4
                    aValues(i*2 + 1) = rangeAINHV;  % HV line
                else
                    aValues(i*2 + 1) = rangeAINLV;  % LV line
                end
            else
                % T7 range
                aValues(i*2 + 1) = rangeAIN;
            end
            aNames(i*2 + 2) = ['AIN' num2str(i) '_RESOLUTION_INDEX'];
            aValues(i*2 + 2) = resolutionAIN;
        end
        LabJack.LJM.eWriteNames(handle, numFrames, aNames, aValues, -1);
    else
        numAIN = 0;
    end

    % Initialize and configure eNames parameters for loop's eNames or
    % eAddresses call
    numFrames = numAIN + readDigital + writeDigital + writeDACs*2;
    aNames = NET.createArray('System.String', numFrames);
    aWrites = NET.createArray('System.Int32', numFrames);
    aNumValues = NET.createArray('System.Int32', numFrames);
    % In this case numFrames is the size of aValue
    aValues = NET.createArray('System.Double', numFrames);

    i = 1;

    % Add analog input reads (AIN 0 to numAIN-1)
    while i <= numAIN
        aNames(i) = ['AIN' num2str(i-1)];
        aWrites(i) = LJM_CONSTANTS.READ;
        aNumValues(i) = 1;
        aValues(i) = 0;
        i = i + 1;
    end

    if readDigital
        % Add digital read
        aNames(i) = 'DIO_STATE';
        aWrites(i) = LJM_CONSTANTS.READ;
        aNumValues(i) = 1;
        aValues(i) = 0;
        i = i + 1;
    end

    if writeDigital
        % Add digital write
        aNames(i) = 'DIO_STATE';
        aWrites(i) = LJM_CONSTANTS.WRITE;
        aNumValues(i) = 1;
        aValues(i) = 0;  % output-low
        i = i + 1;
    end

    if writeDACs
        % Add analog output writes (DAC0-1)
        for j = 0:1
            aNames(i) = ['DAC' num2str(j)];
            aWrites(i) = LJM_CONSTANTS.WRITE;
            aNumValues(i) = 1;
            aValues(i) = 0.0;  % 0.0 V
            i = i + 1;
        end
    end

    % Make arrays of addresses and data types for eAddresses.
    aAddresses = NET.createArray('System.Int32', numFrames);
    aTypes = NET.createArray('System.Int32', numFrames);
    LabJack.LJM.NamesToAddresses(numFrames, aNames, aAddresses, aTypes);

    disp('Test frames:')

    wrStr = '';
    for i = 1:numFrames
        if aWrites(i) == LJM_CONSTANTS.READ
            wrStr = 'READ';
        else
            wrStr = 'WRITE';
        end
        disp(['    ' wrStr ' ' char(aNames(i))])
    end
    disp(['Performing ' num2str(numIterations) ' iterations...']);

    % Initialize time variables
    maxMS = 0;
    minMS = 0;
    totalMS = 0;
    curMS = 0;

    tic;

    % eAddresses/eNames operations loop
    for i = 1:numIterations
        st2 = tic;
        if useAddresses
            LabJack.LJM.eAddresses(handle, numFrames, aAddresses, ...
                aTypes, aWrites, aNumValues, aValues, -1);
        else
            LabJack.LJM.eNames(handle, numFrames, aNames, aWrites, ...
                aNumValues, aValues, -1);
        end
        curMS = toc(st2)*1000;
        if minMS == 0
            minMS = curMS;
        end
        minMS = min(curMS, minMS);
        maxMS = max(curMS, maxMS);
    end

    totalMS = toc*1000;

    disp(['    ' num2str(numIterations) ' iterations performed:'])
    disp(['    Time taken: ' num2str(totalMS) ' ms']);
    disp(['    Average time per iteration: ' ...
          num2str(totalMS/numIterations) ' ms'])
    disp(['    Min / Max time for one iteration: ' num2str(minMS) ...
          ' ms / ' num2str(maxMS) ' ms'])

    if useAddresses
        disp('Last eAddresses results:')
    else
        disp('Last eNames results:')
    end
    for i = 1:numFrames
        if aWrites(i) == LJM_CONSTANTS.READ
            wrStr = 'READ';
        else
            wrStr = 'WRITE';
        end
        disp(['    ' char(aNames(i)) ' (' num2str(aAddresses(i)) ') ' ...
              wrStr ' value : ' num2str(aValues(i))])
    end

    % Close handle
    LabJack.LJM.Close(handle);
catch e
    showErrorMessage(e)
    LabJack.LJM.CloseAll();
end
