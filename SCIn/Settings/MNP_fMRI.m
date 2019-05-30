function h = MNP_fMRI(h,opt)

switch opt
    
    case 'setoptions'
        
    % settings options
    h.SettingsOptions = {'practiceTACTILE','TACTILE_EXPT'};
    
    
    case 'TACTILE_EXPT'

    % set general options
    h = setgeneral(h);
    
    % FILENAME OF SEQUENCE CREATION FUNCTION (without .m)
    h.SeqFun = 'CreateDesign';
    
    %% TRIALS or CONTINUOUS?
    h.Settings.design = 'trials';
    % if continuous, how many trials ahead should be in the player schedule?
    % (applied to stimulation via soundcard only)
    h.Settings.ntrialsahead = 0;  %0 = all trials
    
    %% EXPERIMENTAL CONDIITIONS
    % name the settings that define orthogonal conditions at a different row
    h.Settings.conds = {'oddprob'};
    
    %% Output options
    % save sinwave from all trials as part of stim sequence file
    %h.Settings.savesinwave = 0;
    
    %% BLOCKING/RUN OPTIONS
    % 'divide' = equally divide trials by nblocks; 
    % 'cond' = separate block for each condition
    h.Settings.blockopt = 'divide'; % cond or divide
    % further options for 'divide':
        % number of blocks (containing multiple conditions)
        h.Settings.nblocks = 2; % must integer-divide each value in h.Settings.cond_rep_init
        %distribute conditions equally among blocks
        h.Settings.distblocks = 1;
    % options to start sequence at beginning of every run
    % 'msgbox', 'labjack', 'buttonpress', 'audio' - can have more than one in
    % cell array
    h.Settings.blockstart = {'buttonpress','scannertrig'};%'buttonpress','audio','scannertrig'}; % audio,labjack,audio
    h.Settings.pauseeachblock = 0; % pause after every block?
    % names of any audiofiles
    h.Settings.audiofile = {'instruct.wav'};
    h.Settings.stimcontrol='PsychPortAudio';
    h.Settings.atten_instruct = -10; % add on further attentuation to any other attenuation for audio instructions
    
    % scanner triggers
    h.Settings.num_scanner_trig = 4;
    % wait time between scanner triggers (s)
    h.Settings.waittime_scanner_trig = 1;
    % key corresponding to scanner trigger
    h.Settings.triggeropt = '7&';
    
    %% Condition-independent stimulus parameters - can be superceded by condition-dependent parameters
    % duration of stimulus sequence in seconds
    h.Settings.totdur = 0; 
    % duration of trial in seconds
    h.Settings.trialdur = 1; % if 0, consecutive stimuli will occur with no gap
    % Tactile: number of pulses per trial
    h.Settings.nstim_trial = 1; 
    % Tactile: within-trial wait
    h.Settings.wait=[0]; % one value per nstim 
    
    %% stimulus: tactile
    h.Settings.stim(1).dur = [0.5 0.5]; % duration of stimulus in seconds; modified by oddball settings
    h.Settings.stim(1).patternmethod = 'intensity';% Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.stim(1).patternvalue = [1 0]; % one per stimdur in each cell; one cell per oddball value
    h.Settings.stim(1).durtype = 'reg'; % not needed unless 'rand'
    h.Settings.stim(1).inten = 1; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).f0 = 23; % pitch
    h.Settings.stim(1).inten_type = 'abs'; % either 'dB' or 'abs'
    h.Settings.stim(1).df = 0;
    h.Settings.stim(1).atten = {}; % attenuation level in decibels
    h.Settings.stim(1).attenchan = []; % apply attenuation (e.g. during thresholding) to these chans
    h.Settings.stim(1).control='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    h.Settings.stim(1).chan = [3 4]; 
    h.Settings.stim(1).nrchannels = 8; % total number of channels, e.g. on sound card
    h.Settings.stim(1).Tukey = 0; % Apply Tukey window?
    h.Settings.stim(1).Tukeytype = 1; % 1 = apply to each tone within pattern; 2 = apply to whole pattern
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: do the stimtypes indexed by
    % h.Settings.stimtypeouts_init (oddball method) differ by any other
    % characteristic? E.g. different modalities may requrie different pitches
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% for each number in h.Settings.stimtypeouts_init,
    % Oddball method: intensity, pitch, channel
    h.Settings.oddballmethod = 'channel'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.oddballvalue = {3,4}; % [7 8 3 4]
    h.Settings.monostereo = 'mono'; 
    
    % Oddball method: intensity, index, pitch, channel
    h.Settings.PL.oddballmethod = 'index'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.PL.oddballvalue = {[1 2], [1 2]}; % values to go into h.Seq.signal. One per oddprob row, or leave blank if determined from GUI
    h.Settings.PL.oddballtype = 'roving'; % options: 'roving', 'classical'

     %% SEQUENCES: PL
    % Change probablity (CP): each condition is in rows
    h.Settings.PL.oddprob = [
        % standard (left) vs oddball (right)
        0.6 0.4
        0.8 0.2
        ];
    % index of cols of oddprob that are standards and oddballs. Can be
    % overridden by h.Settings.oddballvalue if using index
    h.Settings.PL.standardind = 1; 
    h.Settings.PL.oddind = 2; 
    % keep oddball trials apart by at least sep_odd standards
    h.Settings.PL.sep_odd = [0 1];%[0 2 0 2 0 2 0 2 0 2 0 2]; % for each CP condition
    % for sep_odd, which indices of h.Settings.oddballvalue to consider
    % each time? (each list will be considered separately)
    h.Settings.PL.sep_odd_ind = {[1 2],[1 2]};
    h.Settings.PL.sep_odd_tol = [1 1]; % set these to be as high as possible (max 1)
    % for each set, ensure a number of leading standards 
    h.Settings.PL.std_lead = [0 0]; % for each CP condition
    % number of sets to randomise together
    h.Settings.PL.n_set = []; % Leave blank to calculate automatically; or one nunmber per CP condition
    % min number of oddballs within each CP condition
    h.Settings.PL.n_odd = [40 40]; % overrides h.Settings.totdur
    % min number of oddballs per randomised set, per CP
    h.Settings.PL.n_odd_set = [20 20]; % overrides h.Settings.totdur
    % randomise sets?
    h.Settings.PL.rand_set = [1 1]; 
    % condition numbers
    h.Settings.PL.condnum = [
        1 2 
        3 4 
        ]; 
    
    %% RESPONSE PARAMETERS
    % record responses during experiment? 0 or 1
    h.Settings.record_response = 1;
    % how to record responses?
    h.Settings.record_response_type = {'all'}; %options: 'all','thistrial','previoustrial'
    %intensity difference
    % buttonpress options: key: keyboard inputs. Blank for no button press
    h.Settings.buttontype='key';
    % range of keyboard presses indicating a recordable response
    h.Settings.buttonopt = {'1!','2@','3#','4$','7&'}; 
    % how early after start of trial can button press trigger the next trial? Empty if programmed
    % ISI
    h.Settings.response_nexttrialmin = 0;
    % when does next trial starts after button press? Empty if programmed
    % ISI
    h.Settings.response_nexttrialwait = 0;
    
    %display RT, 1=yes, 0=no
    h.Settings.displayRT=0;
    % response probe
    h.Settings.RPconds=[2 4 6 8]; % condition numbers to apply this
    h.Settings.RPprob=0.5;
    h.Settings.RPmethod='intensity';
    h.Settings.RPdur = [0.15 0.2 0.15 0.5];
    h.Settings.RPvalue=[1 0 1 0];

    %% THRESHOLDING
    % starting level and step size
    %h.Settings.threshold.startinglevel = 2; % for intensity)
    %h.Settings.threshold.step = 2;
    
    case 'thresholdAUDIO'

    % set general options
    h = setgeneral(h);
    
    % FILENAME OF SEQUENCE CREATION FUNCTION (without .m)
    h.SeqFun = 'SimpleSequence';
    
    %% TRIALS or CONTINUOUS?
    h.Settings.design = 'trials';
    % if continuous, how many trials ahead should be in the player schedule?
    % (applied to stimulation via soundcard only)
    h.Settings.ntrialsahead = 0;  %0 = all trials
    
    %% EXPERIMENTAL CONDIITIONS
    % name the settings that define orthogonal condtions at a different row
    h.Settings.conds = {};
    
    %% Output options
    % save sinwave from all trials as part of stim sequence file
    %h.Settings.savesinwave = 0;
    
    %% BLOCKING/RUN OPTIONS
    % 'divide' = equally divide trials by nblocks; 
    % 'cond' = separate block for each condition
    h.Settings.blockopt = 'cond';
    % further options for 'divide':
        % number of blocks (containing multiple conditions)
    %    h.Settings.nblocks = 2; % must integer-divide each value in h.Settings.cond_rep_init
        %distribute conditions equally among blocks
    %    h.Settings.distblocks = 1;
    % options to start sequence at beginning of every run
    % 'msgbox', 'labjack', 'buttonpress', 'audio' - can have more than one in
    % cell array
    h.Settings.blockstart = {'buttonpress'}; % audio,labjack,audio
    h.Settings.pauseeachblock = 0; % pause after every block?
    % names of any audiofiles
    h.Settings.audiofile = {};
    
    %% Condition-independent stimulus parameters - can be superceded by condition-dependent parameters
    % duration of stimulus sequence in seconds
    h.Settings.totdur = 0; 
    % duration of trial in seconds
    h.Settings.trialdur = 2; % if 0, consecutive stimuli will occur with no gap
    h.Settings.nstim_trial = 2; % set to zero to be determined by stimdur
    h.Settings.wait=[0.3 0]; % within-trial frequency (Hz); one value per nstim 
    
    %% first stimulus: colour dot
    % duration of stimulus in seconds
    h.Settings.stim(1).dur = 0; % modified by oddball settings
    h.Settings.stim(1).patternmethod = '';% Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.stim(1).patternvalue = {}; % one per stimdur in each cell; one cell per oddball value
    h.Settings.stim(1).durtype = 'reg'; % not needed unless 'rand'
    h.Settings.stim(1).inten = 1; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).maxinten = 1; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).f0 = 0; % freq
    h.Settings.stim(1).inten_type = 'abs'; % either 'dB' or 'abs'
    %h.Settings.stim(1).df = 0;
    %h.Settings.stim(1).atten = 0; % attenuation level in decibels OR cell with variable and multiplier
    %h.Settings.stim(1).attenchan = [1 2]; % apply attenuation (e.g. during thresholding) to these chans
    h.Settings.stim(1).control='ptb_visual'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt, ptb_visual
    %h.Settings.stim(1).chan = [1 2]; 
    %h.Settings.stim(1).nrchannels = 2; % total number of channels, e.g. on sound card
    %h.Settings.stim(1).Tukey = 0.25; % Apply Tukey window?
    %h.Settings.stim(1).Tukeytype = 2; % 1 = apply to each tone within pattern; 2 = apply to whole pattern
    h.Settings.stim(1).rectColor = [1 0 0; 1 1 1];
    h.Settings.stim(1).rectSize = [0 0 200 200; 0 0 20 20];
    
    %% first stimulus: fixation dot
    % duration of stimulus in seconds
    h.Settings.stim(2).dur = 0; % modified by oddball settings
    h.Settings.stim(2).patternmethod = '';% Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.stim(2).patternvalue = {}; % one per stimdur in each cell; one cell per oddball value
    h.Settings.stim(2).durtype = 'reg'; % not needed unless 'rand'
    h.Settings.stim(2).inten = 1; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(2).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(2).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(2).maxinten = 1; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(2).f0 = 0; % freq
    h.Settings.stim(2).inten_type = 'abs'; % either 'dB' or 'abs'
    %h.Settings.stim(2).df = 0;
    %h.Settings.stim(2).atten = 0; % attenuation level in decibels OR cell with variable and multiplier
    %h.Settings.stim(2).attenchan = [1 2]; % apply attenuation (e.g. during thresholding) to these chans
    h.Settings.stim(2).control='ptb_visual'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt, ptb_visual
    %h.Settings.stim(2).chan = [1 2]; 
    %h.Settings.stim(2).nrchannels = 2; % total number of channels, e.g. on sound card
    %h.Settings.stim(2).Tukey = 0.25; % Apply Tukey window?
    %h.Settings.stim(2).Tukeytype = 2; % 1 = apply to each tone within pattern; 2 = apply to whole pattern
    h.Settings.stim(2).rectColor = [1 1 1];
    h.Settings.stim(2).rectSize = [0 0 20 20];
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, pitch, channel
    h.Settings.oddballmethod = ''; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.oddballvalue = [];
    h.Settings.oddballtype = ''; % options: 'roving', 'classical'

    
    %% SEQUENCE
    h.Settings.ntrials = 500;
    
    %% RESPONSE PARAMETERS
    % record responses during experiment? 0 or 1
    h.Settings.record_response = 1;
    % buttonpress options: key: keyboard inputs. Blank for no button press
    h.Settings.buttontype='key';
    % range of keyboard presses indicating a recordable response
    h.Settings.buttonopt = {'UpArrow','DownArrow'}; 
    
    %% THRESHOLDING
    % starting level and step size
    h.Settings.threshold.type = 'intensity'; % for intensity
    h.Settings.threshold.stim = 1; % which stim to run adaptive on?
    h.Settings.threshold.stimpart = 1; % which part of that stim to run adaptive on?
    h.Settings.threshold.startinglevel = 0; % for intensity, in dB (e.g. 10); for pitch, in Hz (e.g. 100) 
    h.Settings.threshold.step = 0.01;
    h.Settings.threshold.signalval = [1 2]; % 1 = carrying on increasing; 2 = decrease
    h.Settings.threshold.maxinten = 1; % 0 is the max
    
    
end

function h = setgeneral(h)

%% EQUIPMENT CONTROL
% record EEG, NS: netstation, BV: brainvision, 'serial': serial port
% serial port
%h.Settings.serial = 'COM1';
%h.Settings.record_EEG='labjack_DB15';
%h.Settings.EEGport = 'COM3'; % only needed for 'serial' EEG triggers
%h.Settings.EEGMarkPattern = 0; % mark EEG for every change in stimulus pattern (0 = start of trial only)
%h.Settings.labjack=1; % Use labjack for controlling any equipment?
%h.Settings.stimcontrol='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
% if using PsychPortAudio with more than 2 channels:
    %7.1 soundcard channel numbers:
    %front = 1,2 - can't use these for 7.1 without new connectors
    %C-sub = 3,4 - connect to 1/2
    %rear = 5,6 - connect to 3/4
    %side = 7,8
%h.Settings.labjack_DACport = 0;
%h.Settings.DAC_multiply = 0; % multiply DAC output by this (e.g. to get mA on DS8R)
%h.Settings.stimchan = [5 6]; 
%h.Settings.stimchanforLJ = 0;
%h.Settings.nrchannels = 8; % total number of channels, e.g. on sound card
% channels on stimulator to use; use differenr rows for different pairs
% (e.g. for different conditions). If labjack, this can refer to output
% port(s).

% Tactile: number of pulses in a train
%h.Settings.npulses_train = 1; % set to zero to be determined by stimdur
% Tactile: within-train frequency (Hz)
%h.Settings.p_freq=[];
% Use timer to control frequency of labjack outputs? Otherwise uses software timing.
%h.Settings.labjack_timer=1; % must use if there are a large number of pulses (e.g. >4)

%% AUDIO

% duration of stimulus in seconds
%h.Settings.stimdur = 0.3; % modified by oddball settings
% pitch
%h.Settings.f0 = 500; 
%intensity
%h.Settings.inten = 0; % value between 0 and 1, or if decibels value of <=0
%h.Settings.inten_type = 'dB'; % either 'dB' or 'abs'
%h.Settings.df = 0;
% Audio sampling rate 
h.Settings.fs = 96000; % don't change this
% attenuation level in decibels
%h.Settings.atten = 0; 
%h.Settings.attenchan = [1 2]; % apply attenuation (e.g. during thresholding) to these chans
% stimulate responses?
h.Settings.simulate_response = 0;