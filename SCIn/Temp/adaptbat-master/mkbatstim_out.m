function bat_stim = mkbatstim_out(source, marker, times, fname)

% Generate stim by placing marker sound at times, over a source signal, or silence
%
% REQUIRED INPUTS:
%   * source - either a filename, or scalar indicating total duration of signals
%       e.g. when making a metronome  without underlying sound
%   * marker - either a filename for percussive sound, a vector (arbitrary sound), or a scalar
%   indicating to use a beep with this frequency
%   * times - times to place marker (in seconds)
%   * fname is file path for output file
%
% fs = 44100 assumed
%
% ----------------------
% Original version (mkbatstim.m )written by John Iversen.
% Current version modified by Brian K. Hurley to optionally return stimulus
% waveform. Saving stim audio files to disk also made to be optional.

RETURN_WAVFORM = 1; % 1 if returning stimulus waveform, otherwise 0
WRITE_WAVFILE = 0; % 1 if saving stim to disk, otherwise 0

fs = 44100;

%handle source
if ischar(source),
    [snd, wavFs] = audioread(source);
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
    [marksnd, wavFs] = audioread(marker);
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

marksnd = marksnd(:);
markSamp = length(marksnd);

%truncate
times(times>stimlen-0.5) = []; %remove any beeps within 1/2s of end of stim

% add timed beeps to excerpt
for iT = 1:length(times),
    tBeat = times(iT);
    iBeat = floor(tBeat * fs) + 1;
    snd(iBeat:iBeat+markSamp-1) = snd(iBeat:iBeat+markSamp-1) + marksnd;
end

%rescale within +/- 1
bat_stim = snd / max(abs(snd));

if WRITE_WAVFILE
    % save
    wavwrite(bat_stim, fs, n, fname);
end

if ~RETURN_WAVFORM
    bat_stim = [];
end
