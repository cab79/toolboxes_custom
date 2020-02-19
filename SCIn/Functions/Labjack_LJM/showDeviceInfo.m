function showDeviceInfo(handle)
% SHOWDEVICEINFO  Displays the device's information based on the passed
%   device handle.

[~, devType, connType, serNum, ipAddr, port, maxBytesMB] = ...
    LabJack.LJM.GetHandleInfo(handle, 0, 0, 0, 0, 0, 0);
ipAddrStr = '';
[~, ipAddrStr] = LabJack.LJM.NumberToIP(ipAddr, ipAddrStr);
disp(['Opened a LabJack with Device type: ' num2str(devType) ', ' ...
      'Connection type: ' num2str(connType) ','])
disp(['Serial number: ' num2str(serNum) ', IP address: ' ...
      char(ipAddrStr) ', Port: ' num2str(port) ','])
disp(['Max bytes per MB: ' num2str(maxBytesMB)])
end  % showDeviceInfo end
