function h = Seq_Discrim(h,opt)

switch opt
    
    case 'setoptions'
        
    % settings options
    h.SettingsOptions = {'Discrim_dur'};
    
case 'Discrim_dur'

    % set general options
    h = setgeneral(h);
    
    % FILENAME OF SEQUENCE CREATION FUNCTION (without .m)
    h.SeqFun = 'TriangularSequence';
    
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
    h.Settings.trialdur = inf; % if 0, consecutive stimuli will occur with no gap
    % Tactile: number of pulses per trial
    h.Settings.nstim_trial = 2; % set to zero to be determined by stimdur
    % Tactile: within-trial frequency (Hz) 
    h.Settings.wait=[2 0]; % one value per nstim 
    
    %% first stimulus: audio
    h.Settings.stim(1).patternmethod = 'pitch';% Pattern type method: intensity, pitch. Not supported: channel, duration
    % play with these:
    h.Settings.stim(1).dur = repmat(0.05,1,20); % duration of stimulus in seconds; modified by oddball settings
    h.Settings.stim(1).stimrandind = [];% index of stimdur to randomise. 
    h.Settings.stim(1).patternvalue = repmat([300 0],1,10); % one per stimdur in each cell; one cell per oddball value
    
    h.Settings.stim(1).durtype = 'sequence_rand'; 
    h.Settings.stim(1).inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).f0 = 300; % pitch
    h.Settings.stim(1).inten_type = 'dB'; % either 'dB' or 'abs'
    h.Settings.stim(1).df = 0;
    h.Settings.stim(1).atten = {'inten_diff',1}; % attenuation level in decibels
    h.Settings.stim(1).attenchan = [1 2]; % apply attenuation (e.g. during thresholding) to these chans
    h.Settings.stim(1).control='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    h.Settings.stim(1).chan = [1 2]; 
    h.Settings.stim(1).nrchannels = 2; % total number of channels, e.g. on sound card
    h.Settings.stim(1).Tukey = 0.25; % Apply Tukey window?
    h.Settings.stim(1).Tukeytype = 1; % 1 = apply to each tone within pattern; 2 = apply to whole pattern
    
    % duplicate
    for ns=2:h.Settings.nstim_trial
        h.Settings.stim(ns) = h.Settings.stim(1);
    end
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, index, pitch, channel
    h.Settings.PL.oddballmethod = 'pattern'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.PL.oddballvalue = {[1 2],[3 4]}; % values to go into h.Seq.signal. One per oddprob row, or leave blank if determined from GUI
    h.Settings.PL.oddballtype = 'classical'; % options: 'roving', 'classical'

    %% SEQUENCE
    % Change probablity (CP): each condition is in rows
    h.Settings.PL.oddprob = [
        1/2 1/2
        1/2 1/2
        ];
    % keep oddball trials apart by at least sep_odd standards
    %h.Settings.sep_odd = [0]; % for each CP condition
    % for sep_odd, which indices of h.Settings.oddballvalue to consider
    % each time? (each list will be considered separately)
    %h.Settings.sep_odd_ind = {[1 2]};
    %h.Settings.sep_odd_tol = [1]; % set these to be as high as possible (max 1)
    % for each set, ensure a number of leading standards 
    %h.Settings.std_lead = [0]; % for each CP condition
    % number of sets to randomise together
    %h.Settings.n_set = []; % Leave blank to calculate automatically; or one nunmber per CP condition
    % min number of oddballs within each CP condition
    %h.Settings.n_odd = 10*[12]; % overrides h.Settings.totdur
    % min number of oddballs per randomised set, per CP
    %h.Settings.n_odd_set = 10*[12]; % overrides h.Settings.totdur
    % randomise sets?
    %h.Settings.rand_set = [0]; 
    h.Settings.ntrials = [150 150];
    
    %% RESPONSE PARAMETERS
    % record responses during experiment? 0 or 1
    h.Settings.record_response = 1;
    % how to record responses?
    h.Settings.record_response_type = {'all'}; %options: 'all','thistrial','previoustrial'
    %intensity difference
    % buttonpress options: key: keyboard inputs. Blank for no button press
    h.Settings.buttontype='key';
    % range of keyboard presses indicating a recordable response
    h.Settings.buttonopt = {'LeftArrow','DownArrow','RightArrow'}; 
    % how early after start of trial can button press trigger the next trial? Empty if programmed
    % ISI
    h.Settings.response_nexttrialmin = 2;
    % when does next trial starts after button press? Empty if programmed
    % ISI
    h.Settings.response_nexttrialwait = 0.6:0.2:1.6;

case 'OLD'

    % set general options
    h = setgeneral(h);
    
    % FILENAME OF SEQUENCE CREATION FUNCTION (without .m)
    h.SeqFun = 'TriangularSequence';
    
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
    h.Settings.trialdur = inf; % if 0, consecutive stimuli will occur with no gap
    % Tactile: number of pulses per trial
    h.Settings.nstim_trial = 2; % set to zero to be determined by stimdur
    % Tactile: within-trial frequency (Hz) 
    h.Settings.wait=2*ones(1,h.Settings.nstim_trial); % one value per nstim 
    
    %% first stimulus: audio
    h.Settings.stim(1).patternmethod = 'pitch';% Pattern type method: intensity, pitch. Not supported: channel, duration
    % play with these:
    h.Settings.stim(1).dur = repmat(0.05,1,20); % duration of stimulus in seconds; modified by oddball settings
    h.Settings.stim(1).stimrandind = [2 5 8 11 14 17];% index of stimdur to randomise. 
    h.Settings.stim(1).patternvalue = [315:15:600]; % one per stimdur in each cell; one cell per oddball value
    
    h.Settings.stim(1).durtype = 'sequence_rand'; 
    h.Settings.stim(1).inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).f0 = 0; % pitch
    h.Settings.stim(1).inten_type = 'dB'; % either 'dB' or 'abs'
    h.Settings.stim(1).df = 0;
    h.Settings.stim(1).atten = {'inten_diff',1}; % attenuation level in decibels
    h.Settings.stim(1).attenchan = [1 2]; % apply attenuation (e.g. during thresholding) to these chans
    h.Settings.stim(1).control='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    h.Settings.stim(1).chan = [1 2]; 
    h.Settings.stim(1).nrchannels = 2; % total number of channels, e.g. on sound card
    h.Settings.stim(1).Tukey = 0.25; % Apply Tukey window?
    h.Settings.stim(1).Tukeytype = 2; % 1 = apply to each tone within pattern; 2 = apply to whole pattern
    
    % duplicate
    for ns=2:h.Settings.nstim_trial
        h.Settings.stim(ns) = h.Settings.stim(1);
    end
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, index, pitch, channel
    h.Settings.PL.oddballmethod = 'pattern'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.PL.oddballvalue = {[1 2],[3 4]}; % values to go into h.Seq.signal. One per oddprob row, or leave blank if determined from GUI
    h.Settings.PL.oddballtype = 'classical'; % options: 'roving', 'classical'

    %% SEQUENCE
    % Change probablity (CP): each condition is in rows
    h.Settings.PL.oddprob = [
        1/3 1/3 1/3
        1/3 1/3 1/3
        ];
    % keep oddball trials apart by at least sep_odd standards
    %h.Settings.sep_odd = [0]; % for each CP condition
    % for sep_odd, which indices of h.Settings.oddballvalue to consider
    % each time? (each list will be considered separately)
    %h.Settings.sep_odd_ind = {[1 2]};
    %h.Settings.sep_odd_tol = [1]; % set these to be as high as possible (max 1)
    % for each set, ensure a number of leading standards 
    %h.Settings.std_lead = [0]; % for each CP condition
    % number of sets to randomise together
    %h.Settings.n_set = []; % Leave blank to calculate automatically; or one nunmber per CP condition
    % min number of oddballs within each CP condition
    %h.Settings.n_odd = 10*[12]; % overrides h.Settings.totdur
    % min number of oddballs per randomised set, per CP
    %h.Settings.n_odd_set = 10*[12]; % overrides h.Settings.totdur
    % randomise sets?
    %h.Settings.rand_set = [0]; 
    h.Settings.ntrials = [150 150];
    
    %% RESPONSE PARAMETERS
    % record responses during experiment? 0 or 1
    h.Settings.record_response = 1;
    % how to record responses?
    h.Settings.record_response_type = {'all'}; %options: 'all','thistrial','previoustrial'
    %intensity difference
    % buttonpress options: key: keyboard inputs. Blank for no button press
    h.Settings.buttontype='key';
    % range of keyboard presses indicating a recordable response
    h.Settings.buttonopt = {'LeftArrow','DownArrow','RightArrow'}; 
    % how early after start of trial can button press trigger the next trial? Empty if programmed
    % ISI
    h.Settings.response_nexttrialmin = 5;
    % when does next trial starts after button press? Empty if programmed
    % ISI
    h.Settings.response_nexttrialwait = 0.6:0.2:1.6;
end

function h = setgeneral(h)

%% EQUIPMENT CONTROL
% record EEG, NS: netstation, BV: brainvision, 'serial': serial port
% serial port
h.Settings.serial = '';
h.Settings.record_EEG='';
%h.Settings.EEGport = 'COM3'; % only needed for 'serial' EEG triggers
h.Settings.EEGMarkPattern = 0; % mark EEG for every change in stimulus pattern (0 = start of trial only)
h.Settings.labjack=0; % Use labjack for controlling any equipment?
%h.Settings.stimcontrol='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
%h.Settings.labjack_DACport = 0;
%h.Settings.DAC_multiply = 0; % multiply DAC output by this (e.g. to get mA on DS8R)
%h.Settings.stimchan = [1 2]; 
%h.Settings.stimchanforLJ = 0;
%h.Settings.nrchannels = 2; % total number of channels, e.g. on sound card
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
% simulate responses?
h.Settings.simulate_response = 0;