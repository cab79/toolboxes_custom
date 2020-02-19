%
% Demonstrates I2C communications using a LabJack and .NET. The
% demonstration uses a LJTick-DAC connected to FIO0/FIO1 for the T7 or
% FIO4/FIO5 for the T4, and configures the I2C settings. Then a read,
% write and again a read are performed on the LJTick-DAC EEPROM.
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
        % Configure FIO4 and FIO5 as digital I/O.
        LabJack.LJM.eWriteName(handle, 'DIO_INHIBIT', hex2dec('FFFCF'));
        LabJack.LJM.eWriteName(handle, 'DIO_ANALOG_ENABLE', 0);

        % For the T4, using FIO4 and FIO5 for SCL and SDA pins. FIO0 to
        % FIO3 are reserved for analog inputs, and digital lines are
        % required.

        % SDA pin number = 5 (FIO5)
        LabJack.LJM.eWriteName(handle, 'I2C_SDA_DIONUM', 5);
        % SCL pin number = 4 (FIO4)
        LabJack.LJM.eWriteName(handle, 'I2C_SCL_DIONUM', 4);
    else
        % For the T7 and other devices, using FIO0 and FIO1 for the SCL and
        % SDA pins.

        % SDA pin number = 1 (FIO1)
        LabJack.LJM.eWriteName(handle, 'I2C_SDA_DIONUM', 1);
        % SCL pin number = 0 (FIO0)
        LabJack.LJM.eWriteName(handle, 'I2C_SCL_DIONUM', 0);
    end

    % Speed throttle is inversely proportional to clock frequency. 0 = max.
    % Speed throttle = 65516 (~100 kHz)
    LabJack.LJM.eWriteName(handle, 'I2C_SPEED_THROTTLE', 65516);

    % Options bits:
    %     bit0: Reset the I2C bus.
    %     bit1: Restart w/o stop
    %     bit2: Disable clock stretching.
    % Options = 0
    LabJack.LJM.eWriteName(handle, 'I2C_OPTIONS', 0);

    % Slave Address of the I2C chip = 80 (0x50)
    LabJack.LJM.eWriteName(handle, 'I2C_SLAVE_ADDRESS', 80);

    % Initial read of EEPROM bytes 0-3 in the user memory area. We need a
    % single I2C transmission that writes the chip's memory pointer and
    % then reads the data.

    % Set the number of bytes to transmit
    LabJack.LJM.eWriteName(handle, 'I2C_NUM_BYTES_TX', 1);
    % Set the number of bytes to receive
    LabJack.LJM.eWriteName(handle, 'I2C_NUM_BYTES_RX', 4);

    % TX/RX bytes will go here
    aBytes = NET.createArray('System.Byte', 5); 

    % Set the TX bytes. We are sending 1 byte for the address.
    numBytes = 1;
    aBytes(1) = 0;  % Byte 0: Memory pointer = 0
    LabJack.LJM.eWriteNameByteArray(handle, 'I2C_DATA_TX', numBytes, ...
        aBytes, -1);

    % Do the I2C communications
    LabJack.LJM.eWriteName(handle, 'I2C_GO', 1);

    % Read the RX bytes
    numBytes = 4;
    % aBytes(1) to aBytes(4) will contain the data
    for i = 1:numBytes
        aBytes(i) = 0;
    end
    LabJack.LJM.eReadNameByteArray(handle, 'I2C_DATA_RX', numBytes, ...
        aBytes, -1);

    disp('Read User Memory = ')
    aBytesM = uint8(aBytes);
    disp(aBytesM(1:4))

    % Write EEPROM bytes 0-3 in the user memory area, using the page write
    % technique. Note that page writes are limited to 16 bytes max, and
    % must be aligned with the 16-byte page intervals. For instance, if
    % you start writing at address 14, you can only write two bytes
    % because byte 16 is the start of a new page.

    % Set the number of bytes to transmit
    LabJack.LJM.eWriteName(handle, 'I2C_NUM_BYTES_TX', 5);
    % Set the number of bytes to receive
    LabJack.LJM.eWriteName(handle, 'I2C_NUM_BYTES_RX', 0);

    % Set the TX bytes
    numBytes = 5;
    aBytes(1) = 0;  % Byte 0: Memory pointer = 0
    % Create 4 new random numbers to write (aBytes(2:5)).
    for i = 2:5
        aBytes(i) = randi(255);  % 1 to 255
    end
    LabJack.LJM.eWriteNameByteArray(handle, 'I2C_DATA_TX', numBytes, ...
        aBytes, -1);

    % Do the I2C communications.
    LabJack.LJM.eWriteName(handle, 'I2C_GO', 1);

    disp('Write User Memory = ');
    aBytesM = uint8(aBytes);
    disp(aBytesM(2:5))

    % Final read of EEPROM bytes 0-3 in the user memory area. We need a
    % single I2C transmission that writes the address and then reads the
    % data.

    %Set the number of bytes to transmit
    LabJack.LJM.eWriteName(handle, 'I2C_NUM_BYTES_TX', 1);
    % Set the number of bytes to receive
    LabJack.LJM.eWriteName(handle, 'I2C_NUM_BYTES_RX', 4);

    % Set the TX bytes. We are sending 1 byte for the address.
    numBytes = 1;
    aBytes(1) = 0;  % Byte 0: Memory pointer = 0
    LabJack.LJM.eWriteNameByteArray(handle, 'I2C_DATA_TX', numBytes, ...
        aBytes, -1);

    % Do the I2C communications.
    LabJack.LJM.eWriteName(handle, 'I2C_GO', 1);

    % Read the RX bytes.
    numBytes = 4;
    % aBytes(1) to aBytes(4) will contain the data
    for i = 1:4
        aBytes(i) = 0;
    end
    LabJack.LJM.eReadNameByteArray(handle, 'I2C_DATA_RX', numBytes, ...
        aBytes, -1);

    disp('Read User Memory = ');
    aBytesM = uint8(aBytes);
    disp(aBytesM(1:4))

    % Close handle
    LabJack.LJM.Close(handle);
catch e
    showErrorMessage(e)
    LabJack.LJM.CloseAll();
end
