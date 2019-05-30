function params = bat_analysis_params
% parameters for BAT task psychtoolbox.
% populates parameters for experiment filhandling and stimuli as well as
% for the ZEST algorithm that estimates subjects' thresholds.
%
% 4 Oct 2016 Brian Hurley

% file handling
% NOTE: file paths need to be changed when running BAT on a new computer
params.stim_fpath = '/data0/stimuli/audio/bat/stimuli/'; 
params.data_fpath = '/data1/attmap/attmap_altstim/subject_data/bat_data/'; 
params.fig_path = '/data1/attmap/attmap_altstim/figures/bat/';
params.IBI_fpath = '~/svn/private/matlab/projects/attmap/bat/ITI/'; 
params.matpath = '/data1/attmap/attmap_altstim/matfiles/';
params.tables = '/data1/attmap/attmap_altstim/bat/tables/';

% stims
params.stim_names = {
    'ACL.wav', 'SAS.wav', 'OWA.wav', 'KPS.WAV', 'SMA.WAV', 'PAN.wav', ...
    'HSG.WAV', 'NYN.WAV', 'TJU.wav', 'RRW.WAV', 'HTH.wav', 'OCJ.WAV'
    };
% params.freq = 48000;
% params.nrchannels = 2;
% params.repetitions = 1; % should be >1 if repeating stim during playback
% 
% % specify which tests to administer
% params.perform = {'ibi', 'phase'};

% populate IBI params
params.ibi = params;
params.ibi.outdata_fname = 'bat_ibi_zest.csv';

% %% zest params
% % params for initial p.d.f.
% params.ibi.zest.zestA = 1;
% params.ibi.zest.zestB = 2.50;
% params.ibi.zest.zestC = 2.50;
% params.ibi.zest.zestmaxrange = log10(25); % log of highest threshold value possible; 
% params.ibi.zest.zestminrange = log10(0.01); % log of lowest threshold value possible; 
% params.ibi.zest.zest_init_dl = 10; % initial difference limen (% change in IBI)
% 
% % Parameters for Response function (Weibull function) (P(yes) given stimulus)
% params.ibi.zest.zestfa = 0.10; %gamma in the text, false alarm rate (guess rate for 2AFC)
% params.ibi.zest.zestmiss = 0.02; %delta in the text, miss rate (1/2 inattention rate for 2AFC)
% params.ibi.zest.zestbeta = 6; %beta in the text, slope of response function
% params.ibi.zest.zesteta = 0; %eta in the text, "sweat factor" or response criterion parameter
% 
% % stopping criteria
% params.ibi.zest.sd_stop = 1.15; % threshold algorithm stops when SD of PDF < or = this value
% params.ibi.zest.max_trials = 20; % if SD value criterion not met, stops after 20 trials
% 
% % variables that should be reverted back to non-log units
% params.ibi.zest.zestconvert = {'delta_L', 'sd_pdf'};

%% populate phase params by copying from IBI params & changing as needed
params.phase = params.ibi;
params.phase.outdata_fname = 'bat_phase_zest.csv';

% % assume different starting diff limen w/ phase test, as this one is harder
% % set thresholds to max at 40
% params.phase.zest.zestmaxrange = log10(50);
% params.phase.zest.zest_init_dl = 30; % initial difference limen (% change in phase)