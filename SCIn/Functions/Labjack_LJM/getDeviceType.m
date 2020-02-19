function deviceType = getDeviceType(handle)
% GETDEVICETYPE  Returns the device type value based on the passed device
%   handle.

[~, deviceType, ~, ~, ~, ~, ~] = ...
    LabJack.LJM.GetHandleInfo(handle, 0, 0, 0, 0, 0, 0);
end  % getDeviceType end
