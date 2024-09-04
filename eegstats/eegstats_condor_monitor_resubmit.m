function eegstats_condor_monitor_resubmit(varargin) % local version

pth=pwd;
pthData=fullfile(pwd,'Data');

% monitor for Condor job submission to complete
eegstats_condor_monitor_outputs;
 
quit

function [in,Z] = eegstats_condor_monitor_outputs

nowtime = now;
in = dir('input*.mat');

fin=0;
while fin==0
    Z = dir('output*.mat');
    complete = [Z(:).bytes]>0;
    if length(Z)>=length(in) && all(complete)
        fin=1;
        break
    end
    disp([num2str(sum(complete)) ' complete out of ' num2str(length(complete)) ' generated from ' num2str(length(in)) ' inputs' ])
    pause(300)
end