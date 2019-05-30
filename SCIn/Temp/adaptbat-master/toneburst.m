function [burst, ttl, am] = toneburst(pitch, dur, tag, fs, risefall,type, envelope,pad)
% toneburst  create a toneburst audio signal, with optional trigger and AM
%
% burst = toneburst(pitch, dur, tag, fs, risefall,type)
%
% [burst, ttl, am] = toneburst(pitch, dur, tag, fs, risefall,type,env)
%
% pitch = tone frequency in Hz
% dur   = burst duration in seconds
% tag   = [] or 0: no AM modulation; otherwise AM frequency
% fs    = sample rate (Hz)
% risefall = duration of rise & fall
% type  = carrier type: 'sine', 'square' (default), 'pulse' (single pulse with duration risefall)
% env   = 'exp' applies an exponential damping
%
% JRI 4/2/09

% there's a choice here: burst dur from 50% of rise/fall or from start
% start is simpler as far as duration & tagging, so use that

rfTrig = 0;
%rfTrig = 0.5;

if nargin==0,
  eval(['help ' mfilename])
  return
end

if nargin < 6,
  type = 'square';
end

if tag==0,
    tag = [];
end

%calc durations in samples
dur_p = round(dur * fs);
risefall_p = round(risefall*fs);
top_p = dur_p - 2*(1-rfTrig)*risefall_p;

%burst envelope
bramp_down = cos(0.5 * pi * [1:risefall_p]/risefall_p).^2;
bramp_up = bramp_down(end:-1:1);
top = ones(1, top_p); 
burstEnv = [bramp_up top bramp_down];

% carrier
pp = [0:length(burstEnv)-1];
tt = pp / fs;
switch type,
  case 'square',
    % a low-passed square wave
    N = 7; %uses first 7 non-zero harmonics
    square = square_series(2*pi*pitch*tt, N);
    carrier = square * (1/max(abs(square))); %normalize peak to 1
    headroom = 1;
  case 'sine'
    carrier = sin(2*pi*pitch*tt);
    headroom = 1;
  case 'pulse',
    % currently this only makes one pulse...
    carrier = ones(size(burstEnv));
    burstEnv = zeros(size(carrier));
    burstEnv(1:risefall_p) = 1;
    headroom = 1;
end

% burst
burst = carrier .* burstEnv * headroom;

% ttl
ttl = ones(size(burst));
startFloorIdx = 1:floor(rfTrig*risefall_p);
%minEndFloor = round(0.025 * fs); %need at least this length offset to enable onset
minEndFloor = round(0.001 * fs); %need at least this length offset to enable onset
endFloor = max(minEndFloor, floor(rfTrig*risefall_p));
endFloorIdx = (length(burst)-endFloor+1):length(burst);
ttl(startFloorIdx) = 0;
ttl(endFloorIdx) = 0;

% tag
if ~isempty(tag),
  am = -cos(2*pi*tag*tt);
  %amdepth = 1;
  amdepth = 0.75;
  am = am*(amdepth/2) + (1 - amdepth/2);
  burst = burst .* am;
else
  am = [];
end

%envelope
if nargin > 6,
  range = (0:length(tt)-1)/length(tt)*5; %exp(-5)=0.007
  env = exp(-range);
  burst = burst.*env;
end

function x = square_series(t, N)
% square_series  create a square wave using harmonic series
%
%   x = square_series(t, N)
%
%   use N harmonics, and Feijer window for truncated series
%       note: N is coeficient number, not number of odd
%       terms. 
%
%   e.g. x = square_series(2*pi*f*t, 5), includes 1, 3, 5th
%   harmonics
%
% JRI 4/14/03 (iversen@nsi.edu)

% standard square: sum 1/k*sin(2pi*f*k*t) for k odd
% feijer windows: multipy coefficient by (N-k)/N, for N
% harmonics

x = zeros(size(t));
for k = 1:2:N,
    x = x + (N-k)/N * (1/k) * sin(k*t);
end


