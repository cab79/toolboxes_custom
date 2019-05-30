function mkbatstim(source, marker, times, fname)
% mkbatstim(source, marker, timing, fname)
%
%  generate stim by placing marker sound at times, over a source signal, or silence
%
% source is either a filename, or scalar indicating total duration of signals
%   e.g. when making a metronome  without underlying sound
% marker is either a filename for percussive sound, a vector (arbitrary sound), or a scalar
%   indicating to use a beep with this frequency
% times are times to place marker (in seconds)
% fname is file path for output file
%
% fs = 44100 assumed

fs = 44100;

%handle source
if ischar(source),
  [snd, wavFs, n] = wavread(source);
  if wavFs ~= fs,
    error('source sample rate must be 44100')
  end
  %rms normalize = 1
  snd = snd(:,1);
  snd = snd - mean(snd); %zero mean
  snd = snd / std(snd,1); %rms
else
  pnts = fs*source; %duration of silence in points
  snd = zeros(pnts,1);
  n = 16;
end
stimlen = length(snd) / fs;

%handle marker sound
if ischar(marker),
  [marksnd, wavFs, n] = wavread(marker);
  marksnd = marksnd(:,1);
  if wavFs ~= fs,
    error('marker sample rate must be 44100')
  end
  marksnd = marksnd * 10; %boost relative amplitude
elseif length(marker)==1,
  toneDur = .100;
  riseFall = .005;
  marksnd = toneburst(marker, toneDur, 0, fs, riseFall, 'sine');
  marksnd = marksnd * 4;
else
  %use snd as is
end
%rms normalize = 1
marksnd = marksnd(:);
% marksnd = marksnd - mean(marksnd); %zero mean
% marksnd = marksnd / std(marksnd,1); %rms
markSamp = length(marksnd);

%truncate
times(times>stimlen-0.5) = []; %remove any beeps within 1/2s of end of stim

% add timed beeps to excerpt
for iT = 1:length(times),
  tBeat = times(iT);
  iBeat = floor(tBeat * fs) + 1;
  snd(iBeat:iBeat+markSamp-1) = snd(iBeat:iBeat+markSamp-1) + marksnd;
end

% normalize
% snd = snd - mean(snd); %zero mean
% snd = snd / std(snd,1); %rms

% scale sound to rms of 0.1
% rms = sqrt(mean(snd(:).*snd(:)))
% targetrms = 0.1;
% scale = targetrms / rms;
% snd = snd * scale;

%rescale within +/- 1
snd = snd / max(abs(snd));

% save
wavwrite(snd, fs, n, fname);
