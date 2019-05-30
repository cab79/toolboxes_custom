function timing = batBeatTimes(excerpt, ITI_dir)

% Retrieves beat times for BAT stimuli
%
% REQUIRED INPUTS:
%   * excerpt - name of song without WAV extension, or iso350, iso600, iso850
%   * ITI_dir - file directory where beat timing data for each stim live
%         - Note: if creating iso stims, pass any string as second param
%         (e.g., ' ') as this input will be ignored
% OUTPUT:
%   * timing - structure array witht the following fields:
%       - excerpt = excerpt name       
%       - IBI     = mean IBI, ms
%       - t       = beat times (as tapped by JRI, mean of multiple runs), ms
%       - tISO    = isochronized beat times (used in BAT stimuli), ms
%
% -------------------------   
% written by JRI 5/12/08
% 
% Brian K. Hurley 5 Oct 2016 - modified so that ITI directory is passed 
% in as second parameter rather than hard-coded.

if strncmp(excerpt,'iso',3),
  ioi = str2num(excerpt(4:end));
  require(any(ioi==[350 600 850]),'bad ITI value')
  
  stimLength = 30*1000;
  times = 0:ioi:stimLength;
  
  timing.excerpt = excerpt;
  timing.IBI = ioi;
  timing.t = times;
  timing.tISO = times; %t is isochronous
  
else %songs
  
  
  try
    mfilename = fullfile(ITI_dir, ['ji_' excerpt '2.txt']);
    times = load(mfilename);
  catch
    timing.excerpt = excerpt;
    timing.IBI = [];
    timing.t = [];
    timing.tISO = [];
    return
  end
  
  timing.excerpt = excerpt;
  timing.IBI = mean(diff(times));
  timing.t = times';
  
  %anchor at first tap after 5 second
  tBeepStart = 5000;
  
  times(times<tBeepStart) = [];
  
  timing.tISO = times(1):timing.IBI:times(end);
  
end