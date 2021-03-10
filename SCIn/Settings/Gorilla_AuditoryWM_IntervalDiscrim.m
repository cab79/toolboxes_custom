function h = Gorilla_AuditoryWM_IntervalDiscrim(h,opt)

% To use a T7, first install the software from here: https://labjack.com/sites/default/files/software/LabJack-2019-05-20.exe
% https://labjack.com/support/datasheets/t-series/communication/stream-mode/stream-out
% https://labjack.com/support/datasheets/t-series/communication/stream-mode/stream-out/stream-out-description

switch opt
    
    case 'setoptions'
        
    % settings options
    h.SettingsOptions = {'main'};

    
    case 'main'
    % runs random duration CS trains at (mean) 20Hz for 1s each, once per trial.

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
        %h.Settings.nblocks = 1; % must integer-divide each value in h.Settings.cond_rep_init
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
    h.Settings.trialdur = 1; % if 0, consecutive stimuli will occur with no gap. 'Inf' requires participant to respond to move on to next trial.
    % Tactile: number of pulses per trial
    h.Settings.nstim_trial = 1; 
    % Which stims are targets for behavioural responses?
    h.Settings.target_stims = [1];
    
    %% first stimulus
    
    % initial settings
    nstim=4; % minimum frequency of duration changes. Must divide by 2 to produce integer. In this example: on/off/on/off = 4.
    train_dur = 0.24; % 240ms
    
    % create 9 trial types, gaps of 50ms to 210ms in 20ms steps
    gap = 0.05;
    pip = 0.015;
    step = 0.020;
    h.Settings.stim(1).dur = {
        [pip gap pip train_dur-gap-2*pip] % e.g. 15,50,15,160 for on/off/on/off
        [pip gap+1*step pip train_dur-(gap+1*step)-2*pip] 
        [pip gap+2*step pip train_dur-(gap+2*step)-2*pip]
        [pip gap+3*step pip train_dur-(gap+3*step)-2*pip]
        [pip gap+4*step pip train_dur-(gap+4*step)-2*pip]
        [pip gap+5*step pip train_dur-(gap+5*step)-2*pip]
        [pip gap+6*step pip train_dur-(gap+6*step)-2*pip]
        [pip gap+7*step pip train_dur-(gap+7*step)-2*pip]
        [pip gap+8*step pip train_dur-(gap+8*step)-2*pip] % e.g. 15,210,15,0
        };
    ngaps=nstim/2;
    h.Settings.stim(1).gaprandind = []; % random this index in same order as actual gaps. divide by 2 since we are converting from gap indices in whole train (even numbers) to indices of gaps only
    h.Settings.stim(1).stimrandind{1} = []; % first train's gap order is stimrandind
    
    % other settings
    h.Settings.stim(1).durtype = '';% 'oddballvalue' will select which sequence is presented on each trial according to the values in h.Settings.PL.oddballvalue
    h.Settings.stim(1).inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    if h.Settings.labjack
        % use labjack T7
        h.Settings.stim(1).patternvalue = repmat([5 0],1,ngaps); % one per stimdur in each cell; one cell per oddball value
        h.Settings.stim(1).patternmethod = 'duration';% Pattern type method: intensity, pitch, duration. 
        h.Settings.stim(1).control='T7'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt, LJTick-DAQ, T7
        h.Settings.stim(1).inten_type = 'abs'; % either 'dB' or 'abs'
        h.Settings.stim(1).chan = 'DAC0'; h.Settings.stim(1).chanforLJ = 1;
        h.Settings.stim(1).nrchannels = 1; % total number of channels, e.g. on sound card
        h.Settings.stim(1).npulses_train = 0; % Tactile: number of pulses in a train; set to zero to be determined by stimdur
        h.Settings.stim(1).p_freq=[]; % Tactile: within-train frequency (Hz). only needed if not specifying using stimdur
        h.Settings.stim(1).labjack_timer=1; % Use timer to control frequency of labjack outputs? Otherwise uses software timing. must use if there are a large number of pulses (e.g. >4)
        h.Settings.stim(1).wavetype = 'digital'; % sin, square, step or digital.
        h.Settings.stim(1).f0 = 16384/(2*sum(h.Settings.stim(1).dur{1})); % frequency of output (1 / min duration). 16384 bytes is the max buffer size, and each output requires 2 bytes.
        h.Settings.stim(1).maxinten = 200; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    else
        % use audio
        h.Settings.stim(1).patternvalue = repmat([300 0],1,ngaps); % one per stimdur in each cell; one cell per oddball value
        h.Settings.stim(1).patternmethod = 'pitch';% Pattern type method: intensity, pitch, duration. 
        h.Settings.stim(1).f0 = 300; % pitch
        h.Settings.stim(1).inten_type = 'dB'; % either 'dB' or 'abs'
        h.Settings.stim(1).df = 0;
        h.Settings.stim(1).atten = 0; % attenuation level in decibels
        h.Settings.stim(1).attenchan = [1 2]; % apply attenuation (e.g. during thresholding) to these chans
        h.Settings.stim(1).control='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt, LJTick-DAQ, T7
        h.Settings.stim(1).chan = [1 2]; 
        h.Settings.stim(1).nrchannels = 2; % total number of channels, e.g. on sound card
        h.Settings.stim(1).Tukey = 0.25; % Apply Tukey window?
        h.Settings.stim(1).Tukeytype = 1; % 1 = apply to each tone within pattern; 2 = apply to whole pattern
        h.Settings.stim(1).wavetype = 'sin';
        h.Settings.stim(1).maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    end
    
    % within-trial stimulus onset asychrony (SOA) (s)
    h.Settings.wait=[]; % one value per nstim
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, index, pitch, channel
    h.Settings.PL.oddballmethod = 'pattern'; % this will select different patterns according to h.Settings.PL.oddballvalue. Can use same type for pattern only if oddball intensity is adaptive
    h.Settings.PL.oddballvalue = {1}; % values to go into h.Seq.signal. One per oddprob row, or leave blank if determined from GUI
    h.Settings.PL.nonoddballstimvalue = {0};
    h.Settings.PL.oddballtype = 'classical'; % options: 'roving', 'classical'

    %% SEQUENCE
    % Change probablity (CP): each condition is in rows
    h.Settings.PL.oddprob = [
        1 
        ];
    h.Settings.ntrials = [9]; % if 3 stims per trial should multiples of 3 to allow every combination as programmed in oddballvalue. 
    
    %% RESPONSE PARAMETERS
    % record responses during experiment? 0 or 1
    h.Settings.record_response = 0;
    % how to record responses?
    h.Settings.record_response_type = {'all'}; %options: 'all','thistrial','previoustrial'
    %intensity difference
    % buttonpress options: key: keyboard inputs. Blank for no button press
    h.Settings.buttontype='key';
    % range of keyboard presses indicating a recordable response
    h.Settings.buttonopt = {'LeftArrow','RightArrow'}; 
    % how early after start of trial can button press trigger the next trial? Empty if programmed
    % ISI
    h.Settings.response_nexttrialmin = [];
    % when does next trial start after button press? Empty if programmed
    % ISI
    h.Settings.response_nexttrialwait = [];
    
end

function h = setgeneral(h)

if 0 % turn to 1 to use DS8R with EEG
    %% EQUIPMENT CONTROL: DS8R
    % record EEG, NS: netstation, BV: brainvision, 'serial': serial port
    % serial port
    h.Settings.serial = 'COM1';
    h.Settings.record_EEG='labjack_DB15';
    h.Settings.EEGMarkPattern = 0; % mark EEG for every change in stimulus pattern (0 = start of trial only)
    h.Settings.labjack=1; % Use labjack for controlling any equipment?
    h.Settings.labjack_DACport = 4;
    h.Settings.DAC_multiply = 0.01; % multiply DAC output by this (e.g. to get mA on DS8R)end
end

if 0 % turn to 1 to use DS8R without EEG
    %% EQUIPMENT CONTROL: DS8R
    % record EEG, NS: netstation, BV: brainvision, 'serial': serial port
    % serial port
    h.Settings.serial = '';
    h.Settings.record_EEG='';
    h.Settings.EEGMarkPattern = 0; % mark EEG for every change in stimulus pattern (0 = start of trial only)
    h.Settings.labjack=1; % Use labjack for controlling any equipment?
    h.Settings.labjack_DACport = 4;
    h.Settings.DAC_multiply = 0.01; % multiply DAC output by this (e.g. to get mA on DS8R)
end

if 1 % turn to 1 to use audio without EEG
    %% EQUIPMENT CONTROL: AUDIO
    % serial port
    h.Settings.serial = '';
    h.Settings.record_EEG='';
    h.Settings.EEGMarkPattern = 0; % mark EEG for every change in stimulus pattern (0 = start of trial only)
    h.Settings.labjack=0; % Use labjack for controlling any equipment?
    % Audio sampling rate 
    h.Settings.fs = 96000; % don't change this
    % simulate responses?
    h.Settings.simulate_response = 0;
end
