%
% Demonstrates SPI communication.
%
% You can short MOSI to MISO for testing.
%
% T7:
%     MOSI    FIO2
%     MISO    FIO3
%     CLK     FIO0
%     CS      FIO1
%
% T4:
%     MOSI    FIO6
%     MISO    FIO7
%     CLK     FIO4
%     CS      FIO5
%
% If you short MISO to MOSI, then you will read back the same bytes that
% you write. If you short MISO to GND, then you will read back zeros. If
% you short MISO to VS or leave it unconnected, you will read back 255s.
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

    if getDeviceType(handle) == LJM_CONSTANTS.dtT4
        % Configure FIO4 to FIO7 as digital I/O.
        LabJack.LJM.eWriteName(handle, 'DIO_INHIBIT', hex2dec('FFF0F'));
        LabJack.LJM.eWriteName(handle, 'DIO_ANALOG_ENABLE', 0);

        % Setting CS, CLK, MISO, and MOSI lines for the T4. FIO0 to FIO3
        % are reserved for analog inputs, and SPI requires digital lines.

        % CS is FIO5
        LabJack.LJM.eWriteName(handle, 'SPI_CS_DIONUM', 5);
        % CLK is FIO4
        LabJack.LJM.eWriteName(handle, 'SPI_CLK_DIONUM', 4);
        % MISO is FIO7
        LabJack.LJM.eWriteName(handle, 'SPI_MISO_DIONUM', 7);
        % MOSI is FIO6
        LabJack.LJM.eWriteName(handle, 'SPI_MOSI_DIONUM', 6);
    else
        % Setting CS, CLK, MISO, and MOSI lines for the T7 and other
        % devices.

        % CS is FIO1
        LabJack.LJM.eWriteName(handle, 'SPI_CS_DIONUM', 1);
        % CLK is FIO0
        LabJack.LJM.eWriteName(handle, 'SPI_CLK_DIONUM', 0);
        % MISO is FIO3
        LabJack.LJM.eWriteName(handle, 'SPI_MISO_DIONUM', 3);
        % MOSI is FIO2
        LabJack.LJM.eWriteName(handle, 'SPI_MOSI_DIONUM', 2);
    end

    % Selecting Mode CPHA=1 (bit 0), CPOL=1 (bit 1)
    LabJack.LJM.eWriteName(handle, 'SPI_MODE', 3);

    % Speed Throttle:
    % Valid speed throttle values are 1 to 65536 where 0 = 65536.
    % Configuring Max. Speed (~800 kHz) = 0
    LabJack.LJM.eWriteName(handle, 'SPI_SPEED_THROTTLE', 0);

    % Options
    % bit 0:
    %     0 = Active low clock select enabled
    %     1 = Active low clock select disabled.
    % bit 1:
    %     0 = DIO directions are automatically changed
    %     1 = DIO directions are not automatically changed.
    % bits 2-3: Reserved
    % bits 4-7: Number of bits in the last byte. 0 = 8.
    % bits 8-15: Reserved

    % Enabling active low clock select pin
    LabJack.LJM.eWriteName(handle, 'SPI_OPTIONS', 0);

    % Read back and display the SPI settings
    numFrames = 7;
    aNames = NET.createArray('System.String', numFrames);
    aNames(1) = 'SPI_CS_DIONUM';
    aNames(2) = 'SPI_CLK_DIONUM';
    aNames(3) = 'SPI_MISO_DIONUM';
    aNames(4) = 'SPI_MOSI_DIONUM';
    aNames(5) = 'SPI_MODE';
    aNames(6) = 'SPI_SPEED_THROTTLE';
    aNames(7) = 'SPI_OPTIONS';
    aValues = NET.createArray('System.Double', numFrames);
    LabJack.LJM.eReadNames(handle, numFrames, aNames, aValues, 0);

    disp('SPI Configuration:')
    for i = 1:numFrames
        disp(['  ' char(aNames(i)) ' = ' num2str(aValues(i))]);
    end

    % Write(TX)/Read(TX) 4 bytes
    numBytes = 4;
    LabJack.LJM.eWriteName(handle, 'SPI_NUM_BYTES', numBytes);
    aBytes = NET.createArray('System.Byte', numBytes);

    % Write the bytes
    % Setting array to random values
    for i = 1:numBytes
        aBytes(i) = randi(255);  % 1 to 255
    end
    LabJack.LJM.eWriteNameByteArray(handle, 'SPI_DATA_TX', numBytes, ...
        aBytes, -1);
    % Do the SPI communications
    LabJack.LJM.eWriteName(handle, 'SPI_GO', 1);

    % Display the bytes written
    disp('Write (TX): ')
    disp(uint8(aBytes))

    % Read the bytes
    % Setting array to 0 values
    for i = 1:numBytes
        aBytes(i) = 0;
    end
    LabJack.LJM.eReadNameByteArray(handle, 'SPI_DATA_RX', numBytes, ...
        aBytes, -1);
    % Do the SPI communications
    LabJack.LJM.eWriteName(handle, 'SPI_GO', 1);

    % Display the bytes read
    disp('Read (RX): ')
    disp(uint8(aBytes));

    % Close handle
    LabJack.LJM.Close(handle);
catch e
    showErrorMessage(e)
    LabJack.LJM.CloseAll();
end
