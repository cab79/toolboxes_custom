function params = bat_params

% Parameters for adaptive BAT task.
% Populates various parameters for adaptive BAT PsychToolbox experiment. 
% Filehandling, stimulus, test version, and threshold preferences are 
% specified in this file.
% Users implementing adaptive BAT for the first time on a new computer
% must change file paths in the File Handling section.
%
% written 4 Oct 2016 by Brian K. Hurley

%% File Handling
% 
% NOTE: the below file paths need to be changed when running BAT on a new computer
% location of stimulus files
params.stim_fpath = 'C:\Users\janatalab\Documents\BAT\stimuli'; 
% location where data should be written
params.data_fpath = 'C:\Users\janatalab\Documents\BAT\subject_data\attmap_streamattend';
% location of file containing metronome event (ITI) times
params.ITI_fpath = 'C:\Users\janatalab\Documents\svn\private\matlab\projects\attmap\bat\ITI';

% stims
params.stim_names = {
    'ACL.wav', 'SAS.wav', 'OWA.wav', 'KPS.WAV', 'SMA.WAV', 'PAN.wav', ...
    'HSG.WAV', 'NYN.WAV', 'TJU.wav', 'RRW.WAV', 'HTH.wav', 'OCJ.WAV'
    };
params.freq = 48000;
params.nrchannels = 2;
params.repetitions = 1; % should be >1 if stims should loop during playback

% specify which tests to administer
params.perform = {'ibi', 'phase'};

% populate IBI params
params.ibi = params;
params.ibi.outdata_fname = 'bat_ibi_zest.csv';

%% ZEST Params

% params for initial p.d.f.
params.ibi.zest.zestA = 1;
params.ibi.zest.zestB = 2.50;
params.ibi.zest.zestC = 2.50;
params.ibi.zest.zestmaxrange = log10(25); % log of highest threshold value possible; 
params.ibi.zest.zestminrange = log10(0.01); % log of lowest threshold value possible; 
params.ibi.zest.zest_init_dl = 10; % initial difference to be presented (% change in IBI)

% Parameters for Response function (Weibull function) (P(yes) given stimulus)
params.ibi.zest.zestfa = 0.10; % false alarm rate
params.ibi.zest.zestmiss = 0.02; % miss rate
params.ibi.zest.zestbeta = 6; % slope of response function
params.ibi.zest.zesteta = 0; % "sweat factor" or response criterion parameter

% stopping criteria
params.ibi.zest.sd_stop = 1.15; % threshold algorithm stops when SD of PDF < or = this value
params.ibi.zest.max_trials = 20; % if SD value criterion not met, stops after 20 trials

% variables that should be reverted back to non-log units
params.ibi.zest.zestconvert = {'delta_L', 'sd_pdf'};

%% Populate phase params by copying from IBI params & changing as needed

params.phase = params.ibi;
params.phase.outdata_fname = 'bat_phase_zest.csv';

% larger starting deviation and larger deviation range, as this one is harder
params.phase.zest.zestmaxrange = log10(50);
params.phase.zest.zest_init_dl = 30; % initial difference (% change in phase)

end