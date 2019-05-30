function h = PressureTest(h,opt)

% for a 10min scan:
% each trial is 10s; total trials is 60; 30 painful, 30 non-painful.

% FILENAME OF SEQUENCE CREATION FUNCTION (without .m)
h.SeqFun = 'CreateDesign';

switch opt
    
    case 'setoptions'
        
    % settings options
    h.SettingsOptions = {'test'};
    
    case 'test'
        h = PressureTest(h,'commonsettings');
    
    case 'commonsettings'    
        %% TRIALS or CONTINUOUS?
        h.Settings.design = 'continuous';
        % if continuous, how many trials ahead should be in the player schedule?
        % (applies to stimulation via soundcard only)
        h.Settings.ntrialsahead = 0;  %0 = all trials
        
        %% EXPERIMENTAL CONDIITIONS
        % name the settings that define orthogonal conditions at a different row
        h.Settings.conds = {'oddprob'};

        %% Output options
        % save sinwave from all trials as part of stim sequence file
        h.Settings.savesinwave = 1;

        %% EQUIPMENT CONTROL
        % record EEG, NS: netstation, BV: brainvision, 'serial': serial port
        h.Settings.record_EEG='';
        %h.Settings.spt1_port = 'COM3'; % only needed for 'serial' EEG triggers
        h.Settings.EEGMarkPattern = 0; % mark EEG for every change in stimulus pattern (0 = start of trial only)
        h.Settings.labjack=1; % Use labjack for controlling any equipment?
        h.Settings.DAC_multiply = 1; % multiply DAC output by this (e.g. to get mA on DS8R)

        %% BLOCKING/RUN OPTIONS
        % 'divide' = equally divide trials by nblocks; 
        % 'cond' = separate block for each condition
        h.Settings.blockopt = '';
        % further options for 'divide':
            % number of blocks (containing multiple conditions)
            %h.Settings.nblocks = 1; % must integer-divide each value in h.Settings.cond_rep_init
            %distribute conditions equally among blocks
            %h.Settings.distblocks = 1;
        % options to start sequence at beginning of every run
        % 'msgbox', 'labjack', 'buttonpress', 'audio' - can have more than one in
        % cell array
        h.Settings.blockstart = {'buttonpress'}; 
        % names of any audiofiles
        h.Settings.audiofile = {''}; 

        %% Condition-independent stimulus parameters - can be superceded by condition-dependent parameters
        % duration of stimulus sequence in seconds
        h.Settings.totdur = []; 
        h.Settings.trialdur = 10; % if 0, consecutive stimuli will occur with no gap
        h.Settings.fs = 100; % don't change this
        h.Settings.wavetype = 'step';
        
        %% stimulus: pressure
        h.Settings.stim(1).dur = 3; 
        h.Settings.stim(1).patternmethod = '';% Pattern type method: intensity, pitch. Not supported: channel, duration
        h.Settings.stim(1).patternvalue = []; % one per stimdur
        h.Settings.stim(1).durtype = 'reg'; % not needed unless 'rand'
        h.Settings.stim(1).inten = []; % value between 2 and 1000mA for Digitimer DS8R
        h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
        h.Settings.stim(1).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
        h.Settings.stim(1).maxinten = 4; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
        h.Settings.stim(1).inten_type = 'abs'; % either 'dB' or 'abs'
        h.Settings.stim(1).control='labjack'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
        h.Settings.stim(1).chan = [1]; h.Settings.stim(1).chanforLJ = 0;
        h.Settings.stim(1).nrchannels = 1; % total number of channels, e.g. on sound card
        h.Settings.stim(1).npulses_train = 1; % Tactile: number of pulses in a train; set to zero to be determined by stimdur
        h.Settings.stim(1).p_freq=[];% Tactile: within-train frequency (Hz)
        h.Settings.stim(1).labjack_timer=0; % Use timer to control frequency of labjack outputs? Otherwise uses software timing. must use if there are a large number of pulses (e.g. >4)
        
        h.Settings.stim(1).alignphase = 0;
        % duration of trial in seconds
        h.Settings.stim(1).patternvaluetarget = 0; % amount to add/substract for targets 
        h.Settings.stim(1).monaural = 0; % Monaural beats instead?
        h.Settings.stim(1).f0 = 1; % 10Hz = alpha. Other options: 1Hz, 25Hz, 40Hz.
        % attenuation level in decibels
        h.Settings.stim(1).atten = 0; 
        

        %% Condition-dependent stimulus parameters
        % Condition method: intensity, pitch, channel
        %sd = h.Settings.stimdur;
        h.Settings.conditionmethod = {};
        h.Settings.conditionvalue = {};% Rows: methods. Columns: each stimtype
        % Oddball method: intensity, pitch, channel
        h.Settings.oddballmethod = 'intensity'; % can use same type for pattern only if oddball intensity is adaptive
        % Odball value: DURATION DEVIANTS
        % i.e. all possible duration options and their probability
        % In general, each row is a different stim type, with columns providing
        % within-stimulus temporal patterns of pitch, intensity or channel.
        % Here the options are chosen as:
            % left column = 1st inten/pitch
            % right column = 2nd inten/pitch
        % and the temporal pattern is defined by fc (from either fpitch or finten)
        h.Settings.oddballvalue = {1 1.5};

        %% SEQUENCE
        h.Settings.oddprob = [0.5 0.5];

        h.Settings.oddballtype = 'classical'; % options: 'roving', 'classical'

        % index of oddball value that are standards
        h.Settings.standardind = 1; % does not apply to 'roving oddball' design
        h.Settings.oddind = 2; % does not apply to 'roving oddball' design

        % keep oddball trials apart by at least sep_odd standards
        h.Settings.sep_odd = [0];%[0 0 0 0 0];
        h.Settings.sep_odd_tol = [1]; % set these to be as high as possible (max 1)
        
        % for sep_odd, which indices of h.Settings.oddballvalue to consider
        % each time? (each list will be separated separately)
        h.Settings.sep_odd_ind = {};

        % for each set, ensure a number of leading standards 
        h.Settings.std_lead = [0];%[0 0 0 0 0];

        % number of minimum sets to randomised together
        h.Settings.n_set = [1];%[1 1 1 1 1]; % 1 = use min set size

        % min number of oddballs within each CP condition
        h.Settings.n_odd = [30];%[300 20 20 20 20]; % overrides h.Settings.totdur
        % min number of oddballs per randomised set, per CP
        h.Settings.n_odd_set = [30];%[300 20 20 20 20]; % overrides h.Settings.totdur
        % randomise sets?
        h.Settings.rand_set = [0];%[0 1 1 1 1]; 
        % randomise within sets?
        %h.Settings.rand_within_set = 1;%[0 1 1 1 1]; 
        h.Settings.condnum = [
        1 2
        ]; 

        %% RESPONSE PARAMETERS
        % record responses during experiment? 0 or 1
        h.Settings.record_response = 1;
        h.Settings.record_response_type = {'all'}; %options: 'all','thistrial','previoustrial'
        % buttonpress options: key: keyboard inputs. Blank for no button press
        h.Settings.buttontype='key';
        % range of keyboard presses indicating a recordable response
        h.Settings.buttonopt = {'LeftArrow','RightArrow'}; 

        
        %% ADAPTIVE: General
        % which ones to run? (i.e. indices of h.Settings.adaptive)
        h.Settings.adaptive_general.adapttypes = [1];
        % alternate or randomise runs over types? Alt must have equal number of
        % runs for each adapttype. Cond = one type per CP block
        h.Settings.adaptive_general.seqtype = 'rand'; % 'alt', 'rand', 'cond' 
        h.Settings.adaptive_general.seqtypecond = []; %if 'cond', associate each CP with an adaptive type
        h.Settings.adaptive_general.seqrandblocksize = 12; % should divide the number of trials in a set
        h.Settings.adaptive_general.selectcond.cp = [1]; % which CP condition to run adaptive on?
        h.Settings.adaptive_general.stim = 1; % which stim to run adaptive on?
        h.Settings.adaptive_general.terminate = ''; % 'block' to terminate within each block only. EMPTY '' TO TURN OFF
        h.Settings.adaptive_general.reestimate = ''; % 'block' to re-estimate with wider prior each block

        %% ADAPTIVE 1
        h.Settings.adaptive(1).type = 'detect';
        h.Settings.adaptive(1).updown = [1 1];
        % how many of each to run?
        h.Settings.adaptive(1).nRuns = 10*12;
        % max number of thresh estimates to average over to get overall
        % estimate (for plotting only - red line)
        h.Settings.adaptive(1).av_thresh = [50,75,100];
        h.Settings.adaptive(1).ci_thresh = 20;
        % number of trials each run
        h.Settings.adaptive(1).trialsperrun = 1;
        % adaptive staircase: meanings of the buttonopt
        h.Settings.adaptive(1).buttonfun = {'LeftArrow','RightArrow'};  
        % adaptive staircase: corresponding signal values that would signify a
        % correct answer
        h.Settings.adaptive(1).signalval = [1 2];
        % number of reversals
        h.Settings.adaptive(1).reversals = [4;8;12];
        % stepsize
        h.Settings.adaptive(1).stepsize = [4;2;1]*2;
        % steptype 0 = multiple/divide by stepsize; steptype 1 = add/subtract
        h.Settings.adaptive(1).steptype = 1;
        % stepdir -1 = level decreases intensity; stepdir 1 = level increases intensity
        h.Settings.adaptive(1).stepdir = -1;
        % starting level of adaptive staircase
        h.Settings.adaptive(1).startinglevel = 1; % should be a value in mA. 
        % adapt to omissions of response (not suitable for 2AFC tasks, so set to 0)
        h.Settings.adaptive(1).omit = 0; % 1 = omission is incorrect; 2 = omission is correct
        % which trials (or oddballs if oddonly selected) to start adaptive procedure if there is an omission?
        h.Settings.adaptive(1).startomit = 0;
        % adapt on every trial or only just before an oddball?
        %h.Settings.adaptive.oddonly = 1;
        % max number of trials after oddball that subject must respond (otherwise counts as omitted response)
        %h.Settings.adaptive.resptrials = 4;
        % number of reversals to average over to calculate threshold.
        h.Settings.adaptive(1).reversalForthresh = 6;
        % use mean from the first X responses of each type (high and low)
        %h.Settings.adaptive(1).getmeanfromresponses = 6;
        % maximum amount to adjust the mean if their responses are very
        % incorrect (should be a small fraction, e.g. 1/5th, of the stimulus intensity)
        %h.Settings.adaptive(1).meanadjustmax = 10;
        % maximum amount - for safety
        h.Settings.adaptive(1).levelmax = 4; % should be value in mA. 
        h.Settings.adaptive(1).levelmin = 0.1;

        % number of trials of plot
        h.Settings.plottrials=0;

end

