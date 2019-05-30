function bat_stim = mkbat3(deviation, type, excerpts, params)

% Generates BAT stimuli. Either returns stimulus waveform or writes stim
% audio file to disk depending on input.
%
% Stimulus metronomes are constructed using actual ITIs from JI's tapping
% to define 'on the beat.'
%
% Tempo cases are scaled. For faster tempo additional beeps are added at
% end of the beep train using mean ITI as the IBI.
%
% Phase uses relative phase of surrounding ITIs.
%
% REQUIRED INPUTS:
%   * deviation - percent of perturbation (integer); e.g., -10 = 10%
%       shorter IBI (if tempo test) or 10% negative shift in phase (if phase test)
%   * type - type of manipulation to be performed. 'B' for tempo or 'P' for phase
%   * excerpts - names of musical stimuli (without .WAV extention) to be operated on
%   * params - structure array containing the following fields:
%       - ITI_fpath = file directory where beat timing data for each stim live
%       - stim_fpath = file directory where stimulus audio files live
%
% OUTPUTS:
%   * bat_stim - stimulus audio waveform with metronome beeps superimposed
% ----------------
% Majority of code written by John Iversen for mkbat2.m
% Version 3 updates written 30 Sept 2016 by Brian K. Hurley

toneFreq = 1000;

% transform TYPE and EXCERPTS parameters to cells if not already
if ~iscell(type)
    type = {type};
end

if ~iscell(excerpts)
    excerpts = {excerpts};
end
n_excerpts = length(excerpts);

% iterate through excerpts
for iE = 1:n_excerpts
    
    % get stimulus beat times
    timing = batBeatTimes(excerpts{iE},params.ITI_fpath);
    taps = timing.t;
    ITIm = timing.IBI;
    
    %start taps after 5000ms
    taps(taps<5000) = [];
    
    %create BAT stimuli
    stimlen = 17; %set this to the longest excerpt length
    
    % iterate through deviations
    for iS = 1:length(deviation)
        
        fprintf(2,'%s: %c%s\n',excerpts{iE}, type{iS}, num2str(deviation(iS)));
        
        %determine timing of beeps
        switch type{iS},
            case 'B', %tempo manipulation
                if deviation(iS) == 0 %on beat
                    beats = taps;
                else
                    %scale tap times by deviation
                    tmp = taps - taps(1); %anchor on initial beat
                    tmp = tmp * (1+deviation(iS)/100);
                    beats = tmp + taps(1);
                    if deviation(iS) > 0, %slower, so truncate
                        beats(beats > stimlen*1000) = [];
                    else %faster, add extra beeps at mean ITI
                        extra = (beats(end)+ITIm):ITIm:stimlen*1000;
                        beats = [beats extra];
                    end
                end
                
            case 'P', %phase manipulation
                ITI = diff(taps);
                if deviation(iS) > 0, %shift after the beat
                    shift = [[ITI(1:end) ITIm]*deviation(iS)/100];
                else %shift before the beat
                    shift = [ [ITIm ITI(1:end)]*deviation(iS)/100];
                end
                beats = taps + shift;
            otherwise
        end
        
        beats = beats/1000; %ms to s
        excerptFname = fullfile(params.stim_fpath,[excerpts{iE} '.WAV']);
        
        % generate stimulus audio file name (only used if we are saving audio files)
        outfname = fullfile(params.stim_fpath,[excerpts{iE} '_' type{iS} ...
            num2str(deviation(iS)) '_v3.wav']);
        
        %superimpose beeps at times in beats vector
        bat_stim = mkbatstim_out(excerptFname, toneFreq, beats, outfname);
        
    end %loop over stimulus variants
end %loop over excerpt