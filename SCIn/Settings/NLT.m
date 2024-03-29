function h = NLT(h,opt)

switch opt
    
    case 'setoptions'
        
    % settings options
    h.SettingsOptions = {'audio_test','Ascend','Adaptive','RO','ALPL','ALPL_EEG'};
    
    case 'audio_test'

    % set general options
    h = setgeneral(h);
    h.Settings.labjack=0;
    h.Settings.record_EEG='';
    
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
    h.Settings.blockopt = '';
    % further options for 'divide':
        % number of blocks (containing multiple conditions)
        %h.Settings.nblocks = 2; % must integer-divide each value in h.Settings.cond_rep_init
        %distribute conditions equally among blocks
    %    h.Settings.distblocks = 1;
    % options to start sequence at beginning of every run
    % 'msgbox', 'labjack', 'buttonpress', 'audio' - can have more than one in
    % cell array
    h.Settings.blockstart = {}; % audio,labjack,audio
    h.Settings.pauseeachblock = 0; % pause after every block?
    % names of any audiofiles
    h.Settings.audiofile = {};
    
    %% Condition-independent stimulus parameters - can be superceded by condition-dependent parameters
    % duration of stimulus sequence in seconds
    h.Settings.totdur = 0; 
    % duration of trial in seconds
    h.Settings.trialdur = 1.3; % if 0, consecutive stimuli will occur with no gap
    % Tactile: number of pulses per trial
    h.Settings.nstim_trial = 1; % set to zero to be determined by stimdur
    % Tactile: within-trial frequency (Hz) 
    h.Settings.wait=[0]; % one value per nstim 
    
    %% first stimulus: audio
    % duration of stimulus in seconds
    h.Settings.stim(1).dur = [0.15 0.15]; % duration of stimulus in seconds; modified by oddball settings
    h.Settings.stim(1).patternmethod = 'pitch';% Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.stim(1).patternvalue = {[540 500],[500 540]}; % one per stimdur in each cell; one cell per oddball value
    h.Settings.stim(1).durtype = 'reg'; % not needed unless 'rand'
    h.Settings.stim(1).inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).f0 = 500; % pitch
    h.Settings.stim(1).inten_type = 'dB'; % either 'dB' or 'abs'
    h.Settings.stim(1).df = 0;
    h.Settings.stim(1).atten = {'inten_diff',1}; % attenuation level in decibels
    h.Settings.stim(1).attenchan = [1 2]; % apply attenuation (e.g. during thresholding) to these chans
    h.Settings.stim(1).control='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    h.Settings.stim(1).chan = [1 2]; 
    h.Settings.stim(1).nrchannels = 2; % total number of channels, e.g. on sound card
    h.Settings.stim(1).Tukey = 0.25; % Apply Tukey window?
    h.Settings.stim(1).Tukeytype = 2; % 1 = apply to each tone within pattern; 2 = apply to whole pattern

    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
   %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    
    % Oddball method: intensity, index, pitch, channel
    h.Settings.PL.oddballmethod = 'index'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.PL.oddballvalue = {[1 2],[1 2]}; % values to go into h.Seq.signal. One per oddprob row, or leave blank if determined from GUI
    h.Settings.PL.oddballtype = 'classical'; % options: 'roving', 'classical'

     %% SEQUENCES: PL
    % Change probablity (CP): each condition is in rows
    h.Settings.PL.oddprob = [
        % standard (left) vs oddball (right)
        0.5 0.5
        ];
    % index of cols of oddprob that are standards and oddballs. Can be
    % overridden by h.Settings.oddballvalue if using index
    h.Settings.PL.standardind = 1; 
    h.Settings.PL.oddind = 2; 
    % keep oddball trials apart by at least sep_odd standards
    h.Settings.PL.sep_odd = [0];%[0 2 0 2 0 2 0 2 0 2 0 2]; % for each CP condition
    % for sep_odd, which indices of h.Settings.oddballvalue to consider
    % each time? (each list will be considered separately)
    h.Settings.PL.sep_odd_ind = {[1 2]};
    h.Settings.PL.sep_odd_tol = [1]; % set these to be as high as possible (max 1)
    % for each set, ensure a number of leading standards 
    h.Settings.PL.std_lead = [0]; % for each CP condition
    % number of sets to randomise together
    h.Settings.PL.n_set = []; % Leave blank to calculate automatically; or one nunmber per CP condition
    % min number of oddballs within each CP condition
    h.Settings.PL.n_odd = [100]; % overrides h.Settings.totdur
    % min number of oddballs per randomised set, per CP
    h.Settings.PL.n_odd_set = [100]; % overrides h.Settings.totdur
    % randomise sets?
    h.Settings.PL.rand_set = [1]; 
    % condition numbers
    h.Settings.PL.condnum = [
        1 2
        ]; 
    
    
    %% RESPONSE PARAMETERS
    % record responses during experiment? 0 or 1
    h.Settings.record_response = 0;
   
    case 'Ascend'

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
    h.Settings.nstim_trial = 1; % set to zero to be determined by stimdur
    h.Settings.wait=0; % within-trial frequency (Hz); one value per nstim 
    
    %%stimulus: tactile
    h.Settings.stim(1).dur = 0; % duration of stimulus in seconds; modified by oddball settings
    h.Settings.stim(1).patternmethod = '';% Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.stim(1).patternvalue = []; % one per stimdur
    h.Settings.stim(1).durtype = 'reg'; % not needed unless 'rand'
    h.Settings.stim(1).inten = []; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).maxinten = 200; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).control='LJTick-DAQ'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    h.Settings.stim(1).chan = [6]; h.Settings.stim(1).chanforLJ = 1;
    h.Settings.stim(1).nrchannels = 1; % total number of channels, e.g. on sound card
    h.Settings.stim(1).npulses_train = 1; % Tactile: number of pulses in a train; set to zero to be determined by stimdur
    h.Settings.stim(1).p_freq=[];% Tactile: within-train frequency (Hz)
    h.Settings.stim(1).labjack_timer=1; % Use timer to control frequency of labjack outputs? Otherwise uses software timing. must use if there are a large number of pulses (e.g. >4)
    
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
    % how early after start of trial can button press trigger the next trial? Empty if programmed
    % ISI
    h.Settings.response_nexttrialmin = 0.2;
    % when does next trial starts after button press? Empty if programmed
    % ISI
    h.Settings.response_nexttrialwait = [];
    
    %% THRESHOLDING
    % starting level and step size
    h.Settings.threshold.type = 'intensity'; % for intensity
    h.Settings.threshold.startinglevel = 0; % for intensity)
    h.Settings.threshold.step = 2;
    h.Settings.threshold.signalval = [1 2]; % 2 = carrying on increasing; 1 = decrease
    h.Settings.threshold.maxinten = h.Settings.stim(1).maxinten; % 0 is the max
    
case 'Adaptive'

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
    % name the settings that define orthogonal condtions at a different row
    h.Settings.conds = {};
    
    %% Output options
    % save sinwave from all trials as part of stim sequence file
    %h.Settings.savesinwave = 0;
    
    %% BLOCKING/RUN OPTIONS
    % 'divide' = equally divide trials by nblocks; 
    % 'cond' = separate block for each condition
    h.Settings.blockopt = '';
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
    h.Settings.nstim_trial = 1; % set to zero to be determined by stimdur
    h.Settings.wait=0; % within-trial frequency (Hz); one value per nstim 
    
    %%stimulus: tactile
    h.Settings.stim(1).dur = 0; % duration of stimulus in seconds; modified by oddball settings
    h.Settings.stim(1).patternmethod = '';% Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.stim(1).patternvalue = []; % one per stimdur
    h.Settings.stim(1).durtype = 'reg'; % not needed unless 'rand'
    h.Settings.stim(1).inten = []; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).maxinten = 200; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_type = 'abs'; % either 'dB' or 'abs'
    h.Settings.stim(1).control='LJTick-DAQ'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    h.Settings.stim(1).chan = [6]; h.Settings.stim(1).chanforLJ = 1;
    h.Settings.stim(1).nrchannels = 1; % total number of channels, e.g. on sound card
    h.Settings.stim(1).npulses_train = 1; % Tactile: number of pulses in a train; set to zero to be determined by stimdur
    h.Settings.stim(1).p_freq=[];% Tactile: within-train frequency (Hz)
    h.Settings.stim(1).labjack_timer=1; % Use timer to control frequency of labjack outputs? Otherwise uses software timing. must use if there are a large number of pulses (e.g. >4)
    
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, index, pitch, channel
    h.Settings.PL.oddballmethod = 'index'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.PL.oddballvalue = {[1 2],[1 2]}; % values to go into h.Seq.signal. One per oddprob row, or leave blank if determined from GUI
    h.Settings.PL.oddballtype = 'classical'; % options: 'roving', 'classical'
    
    %% SEQUENCE
    % Change probablity (CP): each condition is in rows
    h.Settings.oddprob = [
        % standard (left) vs oddball (right)
        0.5 0.5
        0.5 0.5
        ];
    % index of row of oddprob that are standards and oddballs. Can be
    % overridden by h.Settings.oddballvalue if using index
    h.Settings.standardind = 1; 
    h.Settings.oddind = 2; 
    % keep oddball trials apart by at least sep_odd standards
    h.Settings.sep_odd = [0 0]; % for each CP condition
    % for sep_odd, which indices of h.Settings.oddballvalue to consider
    % each time? (each list will be considered separately)
    h.Settings.sep_odd_ind = {[1 2],[1 2]};
    h.Settings.sep_odd_tol = [1 1]; % set these to be as high as possible (max 1)
    % for each set, ensure a number of leading standards 
    h.Settings.std_lead = [0 0]; % for each CP condition
    % number of sets to randomise together
    h.Settings.n_set = []; % Leave blank to calculate automatically; or one nunmber per CP condition
    % min number of oddballs within each CP condition
    h.Settings.n_odd = 10*[12 12]; % overrides h.Settings.totdur
    % min number of oddballs per randomised set, per CP
    h.Settings.n_odd_set = 10*[12 12]; % overrides h.Settings.totdur
    % randomise sets?
    h.Settings.rand_set = [0 0]; 
    
    
    %% RESPONSE PARAMETERS
    % record responses during experiment? 0 or 1
    h.Settings.record_response = 1;
    % how to record responses?
    h.Settings.record_response_type = {'all'}; %options: 'all','thistrial','previoustrial'
    %intensity difference
    % buttonpress options: key: keyboard inputs. Blank for no button press
    h.Settings.buttontype='key';
    % range of keyboard presses indicating a recordable response
    h.Settings.buttonopt = {'LeftArrow','RightArrow'}; 
    % how early after start of trial can button press trigger the next trial? Empty if programmed
    % ISI
    h.Settings.response_nexttrialmin = 0.2;
    % when does next trial starts after button press? Empty if programmed
    % ISI
    h.Settings.response_nexttrialwait = 0.6:0.2:1.6;
    
    %% THRESHOLDING
    % starting level and step size
    %h.Settings.threshold.startinglevel = 2; % for intensity)
    %h.Settings.threshold.step = 2;
    
    
    %% ADAPTIVE: General
    % which ones to run? (i.e. indices of h.Settings.adaptive)
    h.Settings.adaptive_general.adapttypes = [1 2];
    % alternate or randomise runs over types? Alt must have equal number of
    % runs for each adapttype. Cond = one type per CP block
    h.Settings.adaptive_general.seqtype = 'rand'; % 'alt', 'rand', 'cond' 
    h.Settings.adaptive_general.seqtypecond = [1 2]; %if 'cond', associate each CP with an adaptive type
    h.Settings.adaptive_general.seqrandblocksize = 12; % should divide the number of trials in a set
    h.Settings.adaptive_general.selectcond.cp = [1:2]; % which CP condition to run adaptive on?
    h.Settings.adaptive_general.stim = 1; % which stim to run adaptive on?
    h.Settings.adaptive_general.terminate = 'block'; % terminate within each block only
    h.Settings.adaptive_general.reestimate = ''; % 'block' to re-estimate with wider prior each block
    
    %% ADAPTIVE 1
    h.Settings.adaptive(1).type = 'detect';
    h.Settings.adaptive(1).updown = [1 1];
    % how many of each to run?
    h.Settings.adaptive(1).nRuns = 100*12;
    % max number of thresh estimates to average over to get overall
    % estimate (for plotting only - red line)
    h.Settings.adaptive(1).av_thresh = [];
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
    % starting level of adaptive staircase
    h.Settings.adaptive(1).startinglevel = 0; % should be a value in mA. 
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
    h.Settings.adaptive(1).levelmax = 100; % should be value in mA. 
    h.Settings.adaptive(1).levelmin = 0;
    h.Settings.adaptive(1).maxtrial = inf;
    h.Settings.adaptive(1).slope_stimratio = 4; % number to divide stimulus level by, to calculate slope of ZEST
    %h.Settings.adaptive(1).expected_change = 10; % smaller value increases the precision of the prior for ZEST and reduces step size of changes in estimates
    h.Settings.adaptive(1).ignoretrials = 0;
    
    %% ADAPTIVE 2
    h.Settings.adaptive(2).type = 'discrim';
    h.Settings.adaptive(2).updown = [1 2];
    % how many of each to run?
    h.Settings.adaptive(2).nRuns = 100*12;
    % max number of thresh estimates to average over to get overall estimate
    h.Settings.adaptive(2).av_thresh = [];
    h.Settings.adaptive(2).ci_thresh = 20;
    % number of trials each run
    h.Settings.adaptive(2).trialsperrun = 1;
    % adaptive staircase: meanings of the buttonopt
    h.Settings.adaptive(2).buttonfun = {'LeftArrow','RightArrow'}; 
    % adaptive staircase: corresponding signal values that would signify a
    % correct answer
    h.Settings.adaptive(2).signalval = [1 2];
    % reversals
    h.Settings.adaptive(2).reversals = [4;8;12];
    % stepsize
    h.Settings.adaptive(2).stepsize = [2;sqrt(2);sqrt(sqrt(2))];
    % steptype 0 = multiple/divide by stepsize; steptype 1 = add/subtract
    h.Settings.adaptive(2).steptype = 0;
    % stepdir -1 = level decreases intensity; stepdir 1 = level increases intensity
    h.Settings.adaptive(2).stepdir = -1;
    % starting level of adaptive staircase
    h.Settings.adaptive(2).startinglevel = 4; % should be a DIFFERENCE value in mA. Keep small as it will increase naturally over time.
    % adapt to omissions of response (not suitable for 2AFC tasks, so set to 0)
    h.Settings.adaptive(2).omit = 0; % 1 = omission is incorrect; 2 = omission is correct
    % which trials (or oddballs if oddonly selected) to start adaptive procedure if there is an omission?
    h.Settings.adaptive(2).startomit = 0;
    % adapt on every trial or only just before an oddball?
    %h.Settings.adaptive.oddonly = 1;
    % max number of trials after oddball that subject must respond (otherwise counts as omitted response)
    %h.Settings.adaptive.resptrials = 4;
    % number of reversals to average over to calculate threshold.
    h.Settings.adaptive(2).reversalForthresh = 6;
    % use mean from the first X responses of each type (high and low)
    %h.Settings.adaptive(2).getmeanfromresponses = 6;
    % maximum amount to adjust the mean if their responses are very
    % incorrect (should be a small fraction, e.g. 1/5th, of the stimulus intensity)
    %h.Settings.adaptive(2).meanadjustmax = 10;
    % maximum amount of the difference value (should be a small fraction, e.g. 1/5th, of the stimulus intensity)
    h.Settings.adaptive(2).levelmax = 100; % should be a DIFFERENCE value in mA.
    h.Settings.adaptive(2).levelmin = 0;
    h.Settings.adaptive(2).maxtrial = inf;
    %h.Settings.adaptive(2).expected_change = 5; % smaller value increases the precision of the prior for ZEST and reduces step size of changes in estimates
    h.Settings.adaptive(2).ignoretrials = 0;
    h.Settings.adaptive(2).eta_divide = 2;
    h.Settings.adaptive(2).slope_stimratio = 1; % number to divide stimulus level by, to calculate slope of ZEST. Should be smaller for thresholding, larger for adjustments during an expt
    
    
    case 'RO' % roving oddball (EEG)

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
    % name the settings that define orthogonal condtions at a different row
    h.Settings.conds = {'oddprob'};
    
    %% Output options
    % save sinwave from all trials as part of stim sequence file
    %h.Settings.savesinwave = 0;
    
    %% BLOCKING/RUN OPTIONS
    % 'divide' = equally divide trials by nblocks; 
    % 'cond' = separate block for each condition
    h.Settings.blockopt = 'divide';
    % further options for 'divide':
        % number of blocks (containing multiple conditions)
        h.Settings.nblocks = 3; % must integer-divide each value in h.Settings.cond_rep_init
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
    h.Settings.trialdur = 0.4; % if 0, consecutive stimuli will occur with no gap
    % Tactile: number of pulses per trial
    h.Settings.nstim_trial = 1; % set to zero to be determined by stimdur
    % Tactile: within-trial frequency (Hz) 
    h.Settings.wait=[0]; % one value per nstim 
    
    %% stimulus: tactile
    % duration of stimulus in seconds
    h.Settings.stim(1).dur = 0; % modified by oddball settings
    h.Settings.stim(1).patternmethod = '';% Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.stim(1).patternvalue = []; % one per stimdur
    h.Settings.stim(1).durtype = 'reg'; % not needed unless 'rand'
    h.Settings.stim(1).inten = []; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).maxinten = 200; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_type = 'abs'; % either 'dB' or 'abs'
    h.Settings.stim(1).control='LJTick-DAQ'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    h.Settings.stim(1).chan = [6]; h.Settings.stim(2).chanforLJ = 1;
    h.Settings.stim(1).nrchannels = 1; % total number of channels, e.g. on sound card
    h.Settings.stim(1).npulses_train = 1; % Tactile: number of pulses in a train; set to zero to be determined by stimdur
    h.Settings.stim(1).p_freq=[];% Tactile: within-train frequency (Hz)
    h.Settings.stim(1).labjack_timer=1; % Use timer to control frequency of labjack outputs? Otherwise uses software timing. must use if there are a large number of pulses (e.g. >4)
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, index, pitch, channel
    h.Settings.oddballmethod = 'index'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.oddballvalue = {[1 2], [3 4], [1 2], [3 4]}; % values to go into h.Seq.signal. One per oddprob row, or leave blank if determined from GUI
    h.Settings.oddballtype = 'roving'; % options: 'roving', 'classical'

    %% SEQUENCE
    % Change probablity (CP): each condition is in rows
    h.Settings.PL.oddprob = [
        0.8 0.2
        0.8 0.2
        0.5 0.5
        0.5 0.5
        ];

% index of row of oddprob that are standards and oddballs. Can be
    % overridden by h.Settings.oddballvalue if using index
    h.Settings.PL.standardind = 1; 
    h.Settings.PL.oddind = 2; 
    % keep oddball trials apart by at least sep_odd standards
    h.Settings.PL.sep_odd = [2 2 1 1];%[0 2 0 2 0 2 0 2 0 2 0 2]; % for each CP condition
    h.Settings.PL.sep_odd_tol = [1 0.8]; % set these to be as high as possible (max 1)
    % for sep_odd, which indices of h.Settings.oddballvalue to consider
    % each time? (each list will be considered separately)
    h.Settings.PL.sep_odd_ind = {[1 2],[1 2],[1 2],[1 2]};
    % for each set, ensure a number of leading standards 
    h.Settings.PL.std_lead = [0 0 0 0]; % for each CP condition
    % number of sets to randomise together
    h.Settings.PL.n_set = [9 9 9 9]; % Leave blank to calculate automatically; or one nunmber per CP condition
    % min number of oddballs within each CP condition
    h.Settings.PL.n_odd = [180 180 180 180]; % overrides h.Settings.totdur
    % min number of oddballs per randomised set, per CP
    h.Settings.PL.n_odd_set = [20 20 20 20]; % overrides h.Settings.totdur
    % randomise sets?
    h.Settings.PL.rand_set = [1 1 1 1]; 
    % condition numbers
    h.Settings.PL.condnum = [
        1 2
        3 4
        5 6
        7 8
        ]; 
    
    %% RESPONSE PARAMETERS
    % record responses during experiment? 0 or 1
    h.Settings.record_response = 0;
%     % how to record responses?
%     h.Settings.record_response_type = {'all'}; %options: 'all','thistrial','previoustrial'
%     %intensity difference
%     % buttonpress options: key: keyboard inputs. Blank for no button press
%     h.Settings.buttontype='key';
%     % range of keyboard presses indicating a recordable response
%     h.Settings.buttonopt = {'LeftArrow','RightArrow'}; 
%     % how early after start of trial can button press trigger the next trial? Empty if programmed
%     % ISI
%     h.Settings.response_nexttrialmin = 0.2;
%     % when does next trial starts after button press? Empty if programmed
%     % ISI
%     h.Settings.response_nexttrialwait = 0.6:0.2:1.6;
    
    
    case 'AL'

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
    h.Settings.blockopt = 'cond';
    % further options for 'divide':
        % number of blocks (containing multiple conditions)
        %h.Settings.nblocks = 2; % must integer-divide each value in h.Settings.cond_rep_init
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
    h.Settings.trialdur = 10; % if 0, consecutive stimuli will occur with no gap
    % Tactile: number of pulses per trial
    h.Settings.nstim_trial = 2; % set to zero to be determined by stimdur
    % Tactile: within-trial frequency (Hz) 
    h.Settings.wait=[1 0]; % one value per nstim 
    
    %% first stimulus: audio
    % duration of stimulus in seconds
    h.Settings.stim(1).dur = [0.15 0.15]; % duration of stimulus in seconds; modified by oddball settings
    h.Settings.stim(1).patternmethod = 'pitch';% Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.stim(1).patternvalue = {[540 500],[460 500]}; % one per stimdur in each cell; one cell per oddball value
    h.Settings.stim(1).durtype = 'reg'; % not needed unless 'rand'
    h.Settings.stim(1).inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).f0 = 500; % pitch
    h.Settings.stim(1).inten_type = 'dB'; % either 'dB' or 'abs'
    h.Settings.stim(1).df = 0;
    h.Settings.stim(1).atten = {'inten_diff',1}; % attenuation level in decibels
    h.Settings.stim(1).attenchan = [1 2]; % apply attenuation (e.g. during thresholding) to these chans
    h.Settings.stim(1).control='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    h.Settings.stim(1).chan = [1 2]; 
    h.Settings.stim(1).nrchannels = 2; % total number of channels, e.g. on sound card
    h.Settings.stim(1).Tukey = 0.25; % Apply Tukey window?
    h.Settings.stim(1).Tukeytype = 2; % 1 = apply to each tone within pattern; 2 = apply to whole pattern

    %% second stimulus: tactile
    % duration of stimulus in seconds
    h.Settings.stim(2).dur = 0; % modified by oddball settings
    h.Settings.stim(2).patternmethod = '';% Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.stim(2).patternvalue = []; % one per stimdur
    h.Settings.stim(2).durtype = 'reg'; % not needed unless 'rand'
    h.Settings.stim(2).inten = []; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(2).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(2).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(2).maxinten = 200; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(2).inten_type = 'abs'; % either 'dB' or 'abs'
    h.Settings.stim(2).control='LJTick-DAQ'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    h.Settings.stim(2).chan = [6]; h.Settings.stim(2).chanforLJ = 1;
    h.Settings.stim(2).nrchannels = 1; % total number of channels, e.g. on sound card
    h.Settings.stim(2).npulses_train = 1; % Tactile: number of pulses in a train; set to zero to be determined by stimdur
    h.Settings.stim(2).p_freq=[];% Tactile: within-train frequency (Hz)
    h.Settings.stim(2).labjack_timer=1; % Use timer to control frequency of labjack outputs? Otherwise uses software timing. must use if there are a large number of pulses (e.g. >4)
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
   %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, index, pitch, channel
    h.Settings.oddballmethod = 'index'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.oddballvalue = {[1 2], [1 2],  [1 2], [2 1], [1 2],  [1 2], [1 2], [2 1], [1 2],  [1 2], [1 2], [2 1]}; % values to go into h.Seq.signal. One per oddprob row, or leave blank if determined from GUI
    h.Settings.oddballtype = 'classical'; % options: 'roving', 'classical'

    %% SEQUENCE
    % Change probablity (CP): each condition is in rows
    h.Settings.AL.oddprob = [
        % standard (left) vs oddball (right)
%         0.8 0.2
%         0.5 0.5
%         0.8 0.2
%         0.8 0.2
%         0.5 0.5
%         0.8 0.2
%         0.8 0.2
%         0.5 0.5
%         0.8 0.2
        0.5 0.5
        0.8 0.2
        0.8 0.2
        0.8 0.2
        0.8 0.2
        ];
%     % index of row of oddprob that are standards and oddballs. Can be
%     % overridden by h.Settings.oddballvalue if using index
%     h.Settings.standardind = 1; 
%     h.Settings.oddind = 2; 
%     % keep oddball trials apart by at least sep_odd standards
%     h.Settings.sep_odd = [2 0 2 2 0 2 2 0 2]; % for each CP condition
%     % for sep_odd, which indices of h.Settings.oddballvalue to consider
%     % each time? (each list will be considered separately)
%     h.Settings.sep_odd_ind = {[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2]};
%     % for each set, ensure a number of leading standards 
%     h.Settings.std_lead = [0 0 0 0 0 0 0 0 0]; % for each CP condition
%     % number of sets to randomise together
%     h.Settings.n_set = []; % Leave blank to calculate automatically; or one nunmber per CP condition
%     % min number of oddballs within each CP condition
%     h.Settings.n_odd = [4, 4, 4,  8, 8, 8,  12, 12, 12]; % overrides h.Settings.totdur
%     % min number of oddballs per randomised set, per CP
%     h.Settings.n_odd_set = h.Settings.n_odd; % overrides h.Settings.totdur
%     % randomise sets?
%     h.Settings.rand_set = [1 1 1 1 1 1 1 1 1]; 
%     % condition numbers
%     h.Settings.condnum = [
%         1 2
%         3 4
%         5 6
%         1 2
%         3 4
%         5 6
%         1 2
%         3 4
%         5 6
%         ]; 
% index of row of oddprob that are standards and oddballs. Can be
    % overridden by h.Settings.oddballvalue if using index
    h.Settings.AL.standardind = 1; 
    h.Settings.AL.oddind = 2; 
    % keep oddball trials apart by at least sep_odd standards
    h.Settings.AL.sep_odd = [0 2 2 2 2];%[0 2 0 2 0 2 0 2 0 2 0 2]; % for each CP condition
    % for sep_odd, which indices of h.Settings.oddballvalue to consider
    % each time? (each list will be considered separately)
    h.Settings.AL.sep_odd_ind = {[1 2],[1 2],[1 2],[1 2],[1 2]};
    h.Settings.AL.sep_odd_tol = [1 1 1 1 1]; % set these to be as high as possible (max 1)
    % for each set, ensure a number of leading standards 
    h.Settings.AL.std_lead = [0 0 0 0 0]; % for each CP condition
    % number of sets to randomise together
    h.Settings.AL.n_set = []; % Leave blank to calculate automatically; or one nunmber per CP condition
    % min number of oddballs within each CP condition
    h.Settings.AL.n_odd = [12 12 12 12 12]; % overrides h.Settings.totdur
    % min number of oddballs per randomised set, per CP
    h.Settings.AL.n_odd_set = h.Settings.AL.n_odd; % overrides h.Settings.totdur
    % randomise sets?
    h.Settings.AL.rand_set = [0 0 0 0 0]; 
    % condition numbers
    h.Settings.AL.condnum = [
        3 4
        5 6
        1 2
        5 6
        1 2
        ]; 
    
    %% associative learning experiments
    h.Settings.AL.pairing = {
      % cue 1  cue 2
        [1], [2]; %pairing option 1 (e.g. most porbable in some blocks)
        [2], [1]; %pairing option 2 (e.g. least probable)
        }; % more than one number per pairing/cue: will be randomly assigned to each cue
    % if there is more than one number per pairing/cue, what is their
    % percent representation?
    h.Settings.AL.probstim = {1,1,1,1,1,1}; % e.g. low intensity diff vs. high intensity diff
    %h.Settings.assoc.probstim = {[0.2 0.8],[1 0],[0.5 0.5],[0.5 0.5],[0.2 0.8],[1 0]}; % e.g. low intensity diff vs. high intensity diff
    % for each condnum, which pair to use?
    h.Settings.AL.pair = [1 2 1 2 2 1];
    % for each stimtype (unique value) within h.Settings.assoc.pairing, 
    % what inten_diff multiplier to use?
    h.Settings.AL.intenstim = [1 -1];
    h.Settings.AL.stimnums = [1 2]; % which two stimulus numbers relate to the two levels of h.Seq.signal?
    h.Settings.AL.stimpart = 1; % which part of the discriminative stim to run adaptive on?
    
    %% RESPONSE PARAMETERS
    % record responses during experiment? 0 or 1
    h.Settings.record_response = 1;
    % how to record responses?
    h.Settings.record_response_type = {'all'}; %options: 'all','thistrial','previoustrial'
    %intensity difference
    % buttonpress options: key: keyboard inputs. Blank for no button press
    h.Settings.buttontype='key';
    % range of keyboard presses indicating a recordable response
    h.Settings.buttonopt = {'LeftArrow','RightArrow'}; 
    % how early after start of trial can button press trigger the next trial? Empty if programmed
    % ISI
    h.Settings.response_nexttrialmin = 0.2;
    % when does next trial starts after button press? Empty if programmed
    % ISI
    h.Settings.response_nexttrialwait = 0.6:0.2:1.6;
    
    %% THRESHOLDING
    % starting level and step size
    %h.Settings.threshold.startinglevel = 2; % for intensity)
    %h.Settings.threshold.step = 2;
    
    %% ADAPTIVE: General
    % which ones to run? (i.e. indices of h.Settings.adaptive)
    h.Settings.adaptive_general.adapttypes = [1 2];
    % alternate or randomise runs over types? Alt must have equal number of
    % runs for each adapttype. Cond = one type per CP block
    h.Settings.adaptive_general.seqtype = 'rand'; % 'alt', 'rand', 'cond' 
    h.Settings.adaptive_general.seqtypecond = []; %if 'cond', associate each CP with an adaptive type
    h.Settings.adaptive_general.seqrandblocksize = 12; % should divide the number of trials in a set
    h.Settings.adaptive_general.selectcond.cp = [1:5]; % which CP condition to run adaptive on?
    h.Settings.adaptive_general.stim = 2; % which stim to run adaptive on?
    h.Settings.adaptive_general.terminate = ''; % 'block' to terminate within each block only. EMPTY '' TO TURN OFF
    h.Settings.adaptive_general.reestimate = ''; % 'block' to re-estimate with wider prior each block
    
    %% ADAPTIVE 1
    h.Settings.adaptive(1).type = 'detect';
    h.Settings.adaptive(1).updown = [1 1];
    % how many of each to run?
    h.Settings.adaptive(1).nRuns = 100*12;
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
    h.Settings.adaptive(1).startinglevel = 0; % should be a value in mA. 
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
    h.Settings.adaptive(1).levelmax = 100; % should be value in mA. 
    h.Settings.adaptive(1).levelmin = 2;
    
    %% ADAPTIVE 2
    h.Settings.adaptive(2).type = 'discrim';
    h.Settings.adaptive(2).updown = [1 2];
    % how many of each to run?
    h.Settings.adaptive(2).nRuns = 100*12;
    % max number of thresh estimates to average over to get overall estimate
    h.Settings.adaptive(2).av_thresh = [50,75,100];
    h.Settings.adaptive(2).ci_thresh = 20;
    % number of trials each run
    h.Settings.adaptive(2).trialsperrun = 1;
    % adaptive staircase: meanings of the buttonopt
    h.Settings.adaptive(2).buttonfun = {'LeftArrow','RightArrow'}; 
    % adaptive staircase: corresponding signal values that would signify a
    % correct answer
    h.Settings.adaptive(2).signalval = [2 1];
    % reversals
    h.Settings.adaptive(2).reversals = [4;8;12];
    % stepsize
    h.Settings.adaptive(2).stepsize = [2;sqrt(2);sqrt(sqrt(2))];
    % steptype 0 = multiple/divide by stepsize; steptype 1 = add/subtract
    h.Settings.adaptive(2).steptype = 0;
    % stepdir -1 = level decreases intensity; stepdir 1 = level increases intensity
    h.Settings.adaptive(2).stepdir = -1;
    % starting level of adaptive staircase
    h.Settings.adaptive(2).startinglevel = 4; % should be a DIFFERENCE value in mA. Keep small as it will increase naturally over time.
    % adapt to omissions of response (not suitable for 2AFC tasks, so set to 0)
    h.Settings.adaptive(2).omit = 0; % 1 = omission is incorrect; 2 = omission is correct
    % which trials (or oddballs if oddonly selected) to start adaptive procedure if there is an omission?
    h.Settings.adaptive(2).startomit = 0;
    % adapt on every trial or only just before an oddball?
    %h.Settings.adaptive.oddonly = 1;
    % max number of trials after oddball that subject must respond (otherwise counts as omitted response)
    %h.Settings.adaptive.resptrials = 4;
    % number of reversals to average over to calculate threshold.
    h.Settings.adaptive(2).reversalForthresh = 6;
    % use mean from the first X responses of each type (high and low)
    %h.Settings.adaptive(2).getmeanfromresponses = 6;
    % maximum amount to adjust the mean if their responses are very
    % incorrect (should be a small fraction, e.g. 1/5th, of the stimulus intensity)
    %h.Settings.adaptive(2).meanadjustmax = 10;
    % maximum amount of the difference value (should be a small fraction, e.g. 1/5th, of the stimulus intensity)
    h.Settings.adaptive(2).levelmax = 50; % should be a DIFFERENCE value in mA.
    h.Settings.adaptive(2).levelmin = 0.1;
    
   case 'ALPL'

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
    h.Settings.blockopt = 'cond';
    % further options for 'divide':
        % number of blocks (containing multiple conditions)
        %h.Settings.nblocks = 2; % must integer-divide each value in h.Settings.cond_rep_init
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
    h.Settings.wait=[1 0]; % one value per nstim 
    
    %% first stimulus: audio
    % duration of stimulus in seconds
    h.Settings.stim(1).dur = [0.15 0.15]; % duration of stimulus in seconds; modified by oddball settings
    h.Settings.stim(1).patternmethod = 'pitch';% Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.stim(1).patternvalue = {[540 500],[500 540]}; % one per stimdur in each cell; one cell per oddball value
    h.Settings.stim(1).durtype = 'reg'; % not needed unless 'rand'
    h.Settings.stim(1).inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).f0 = 500; % pitch
    h.Settings.stim(1).inten_type = 'dB'; % either 'dB' or 'abs'
    h.Settings.stim(1).df = 0;
    h.Settings.stim(1).atten = {'inten_diff',1}; % attenuation level in decibels
    h.Settings.stim(1).attenchan = [1 2]; % apply attenuation (e.g. during thresholding) to these chans
    h.Settings.stim(1).control='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    h.Settings.stim(1).chan = [1 2]; 
    h.Settings.stim(1).nrchannels = 2; % total number of channels, e.g. on sound card
    h.Settings.stim(1).Tukey = 0.25; % Apply Tukey window?
    h.Settings.stim(1).Tukeytype = 2; % 1 = apply to each tone within pattern; 2 = apply to whole pattern

    %% second stimulus: tactile
    % duration of stimulus in seconds
    h.Settings.stim(2).dur = 0; % modified by oddball settings
    h.Settings.stim(2).patternmethod = '';% Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.stim(2).patternvalue = []; % one per stimdur
    h.Settings.stim(2).durtype = 'reg'; % not needed unless 'rand'
    h.Settings.stim(2).inten = []; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(2).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(2).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(2).maxinten = 200; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(2).inten_type = 'abs'; % either 'dB' or 'abs'
    h.Settings.stim(2).control='LJTick-DAQ'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    h.Settings.stim(2).chan = [6]; h.Settings.stim(2).chanforLJ = 1;
    h.Settings.stim(2).nrchannels = 1; % total number of channels, e.g. on sound card
    h.Settings.stim(2).npulses_train = 1; % Tactile: number of pulses in a train; set to zero to be determined by stimdur
    h.Settings.stim(2).p_freq=[];% Tactile: within-train frequency (Hz)
    h.Settings.stim(2).labjack_timer=1; % Use timer to control frequency of labjack outputs? Otherwise uses software timing. must use if there are a large number of pulses (e.g. >4)
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
   %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    
    % Oddball method: intensity, index, pitch, channel
    h.Settings.PL.oddballmethod = 'index'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.PL.oddballvalue = {[1 2], [1 2], [1 2], [2 1], [1 2], [2 1], [1 2], [1 2]}; % values to go into h.Seq.signal. One per oddprob row, or leave blank if determined from GUI
    h.Settings.PL.oddballtype = 'classical'; % options: 'roving', 'classical'

     %% SEQUENCES: PL
    % Change probablity (CP): each condition is in rows
    h.Settings.PL.oddprob = [
        % standard (left) vs oddball (right)
        0.5 0.5
        0.8 0.2
        0.5 0.5
        0.8 0.2
        0.5 0.5
        0.8 0.2
        0.5 0.5
        0.8 0.2
        ];
    % index of cols of oddprob that are standards and oddballs. Can be
    % overridden by h.Settings.oddballvalue if using index
    h.Settings.PL.standardind = 1; 
    h.Settings.PL.oddind = 2; 
    % keep oddball trials apart by at least sep_odd standards
    h.Settings.PL.sep_odd = [0 2 0 2 0 2 0 2];%[0 2 0 2 0 2 0 2 0 2 0 2]; % for each CP condition
    % for sep_odd, which indices of h.Settings.oddballvalue to consider
    % each time? (each list will be considered separately)
    h.Settings.PL.sep_odd_ind = {[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2]};
    h.Settings.PL.sep_odd_tol = [1 1 1 1 1 1 1 1]; % set these to be as high as possible (max 1)
    % for each set, ensure a number of leading standards 
    h.Settings.PL.std_lead = [0 0 0 0 0 0 0 0]; % for each CP condition
    % number of sets to randomise together
    h.Settings.PL.n_set = []; % Leave blank to calculate automatically; or one nunmber per CP condition
    % min number of oddballs within each CP condition
    h.Settings.PL.n_odd = [25 10 25 10 25 10 25 10]; % overrides h.Settings.totdur
    % min number of oddballs per randomised set, per CP
    h.Settings.PL.n_odd_set = [25 10 25 10 25 10 25 10]; % overrides h.Settings.totdur
    % randomise sets?
    h.Settings.PL.rand_set = [0 0 0 0 0 0 0 0]; 
    % condition numbers
    h.Settings.PL.condnum = [
        5 6
        3 4 % low intensity more likely
        5 6
        1 2 % high intensity more likely
        5 6
        1 2 % high intensity more likely
        5 6
        3 4 % low intensity more likely
        ]; 
    
    %% SEQUENCE
    
    % Oddball method: intensity, index, pitch, channel
    h.Settings.AL.oddballmethod = 'index'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.AL.oddballvalue = {[2 1], [2 1], [1 2], [1 2], [2 1], [2 1], [1 2], [1 2]}; % values to go into h.Seq.signal. One per oddprob row, or leave blank if determined from GUI
    h.Settings.AL.oddballtype = 'classical'; % options: 'roving', 'classical'
    % Change probablity (CP): each condition is in rows
    h.Settings.AL.oddprob = [
        % standard (left) vs oddball (right)
        0.8 0.2
        0.8 0.2
        0.8 0.2
        0.8 0.2
        0.8 0.2
        0.8 0.2
        0.8 0.2
        0.8 0.2
        ];
%     % index of row of oddprob that are standards and oddballs. Can be
%     % overridden by h.Settings.oddballvalue if using index
%     h.Settings.standardind = 1; 
%     h.Settings.oddind = 2; 
%     % keep oddball trials apart by at least sep_odd standards
%     h.Settings.sep_odd = [2 0 2 2 0 2 2 0 2]; % for each CP condition
%     % for sep_odd, which indices of h.Settings.oddballvalue to consider
%     % each time? (each list will be considered separately)
%     h.Settings.sep_odd_ind = {[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2]};
%     % for each set, ensure a number of leading standards 
%     h.Settings.std_lead = [0 0 0 0 0 0 0 0 0]; % for each CP condition
%     % number of sets to randomise together
%     h.Settings.n_set = []; % Leave blank to calculate automatically; or one nunmber per CP condition
%     % min number of oddballs within each CP condition
%     h.Settings.n_odd = [4, 4, 4,  8, 8, 8,  12, 12, 12]; % overrides h.Settings.totdur
%     % min number of oddballs per randomised set, per CP
%     h.Settings.n_odd_set = h.Settings.n_odd; % overrides h.Settings.totdur
%     % randomise sets?
%     h.Settings.rand_set = [1 1 1 1 1 1 1 1 1]; 
%     % condition numbers
%     h.Settings.condnum = [
%         1 2
%         3 4
%         5 6
%         1 2
%         3 4
%         5 6
%         1 2
%         3 4
%         5 6
%         ]; 
% index of row of oddprob that are standards and oddballs. Can be
    % overridden by h.Settings.oddballvalue if using index
    h.Settings.AL.standardind = 1; 
    h.Settings.AL.oddind = 2; 
    % keep oddball trials apart by at least sep_odd standards
    h.Settings.AL.sep_odd = [2 2 2 2 2 2 2 2];%[0 2 0 2 0 2 0 2 0 2 0 2]; % for each CP condition
    % for sep_odd, which indices of h.Settings.oddballvalue to consider
    % each time? (each list will be considered separately)
    h.Settings.AL.sep_odd_ind = {[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2]};
    h.Settings.AL.sep_odd_tol = [1 1 1 1 1 1 1 1]; % set these to be as high as possible (max 1)
    % for each set, ensure a number of leading standards 
    h.Settings.AL.std_lead = [0 0 0 0 0 0 0 0]; % for each CP condition
    % number of sets to randomise together
    h.Settings.AL.n_set = []; % Leave blank to calculate automatically; or one nunmber per CP condition
    % min number of oddballs within each CP condition
    h.Settings.AL.n_odd = [10 10 10 10 10 10 10 10]; % overrides h.Settings.totdur
    % min number of oddballs per randomised set, per CP
    h.Settings.AL.n_odd_set = [10 10 10 10 10 10 10 10]; % overrides h.Settings.totdur
    % randomise sets?
    h.Settings.AL.rand_set = [0 0 0 0 0 0 0 0]; 
    % condition numbers
    h.Settings.AL.condnum = [
        
        % block pair 1
        3 4 % pairing 2, AL
        3 4 % pairing 2, ALPL, PL low, AL prior pairing consistent
        
        % block pair 2
        1 2 % pairing 1, AL
        1 2 % pairing 1, ALPL, PL high, AL prior pairing consistent
        
        % block pair 3
        3 4 % pairing 2, AL
        3 4 % pairing 1, ALPL, PL high, AL prior pairing consistent
        
        % block pair 4
        1 2 % pairing 1, AL
        1 2 % pairing 2, ALPL, PL low, AL prior pairing consistent
        
        ]; 
    
    
    %% associative learning experiments
    h.Settings.AL.pairing = {
      % cue 1  cue 2
        [1], [2]; %pairing option 1 (e.g. most porbable in some blocks)
        [2], [1]; %pairing option 2 (e.g. least probable)
        }; % more than one number per pairing/cue: will be randomly assigned to each cue
    % if there is more than one number per pairing/cue, what is their
    % percent representation?
    h.Settings.AL.probstim = {1,1,1,1}; % e.g. low intensity diff vs. high intensity diff
    %h.Settings.assoc.probstim = {[0.2 0.8],[1 0],[0.5 0.5],[0.5 0.5],[0.2 0.8],[1 0]}; % e.g. low intensity diff vs. high intensity diff
    % for each condnum, which pair to use?
    h.Settings.AL.pair = [1 1 2 2];
    % for each stimtype (unique value) within h.Settings.assoc.pairing, 
    % what inten_diff multiplier to use?
    h.Settings.AL.intenstim = [1 -1];
    h.Settings.AL.stimnums = [1 2]; % which two stimulus numbers relate to the two levels of h.Seq.signal?
    h.Settings.AL.stimpart = 1; % which part of the discriminative stim to run adaptive on?
    
    %% RESPONSE PARAMETERS
    % record responses during experiment? 0 or 1
    h.Settings.record_response = 1;
    % how to record responses?
    h.Settings.record_response_type = {'all'}; %options: 'all','thistrial','previoustrial'
    %intensity difference
    % buttonpress options: key: keyboard inputs. Blank for no button press
    h.Settings.buttontype='key';
    % range of keyboard presses indicating a recordable response
    h.Settings.buttonopt = {'LeftArrow','RightArrow'}; 
    % how early after start of trial can button press trigger the next trial? Empty if programmed
    % ISI
    h.Settings.response_nexttrialmin = 0.2;
    % when does next trial starts after button press? Empty if programmed
    % ISI
    h.Settings.response_nexttrialwait = 0.6:0.2:1.6;
    
    %% THRESHOLDING
    % starting level and step size
    %h.Settings.threshold.startinglevel = 2; % for intensity)
    %h.Settings.threshold.step = 2;
    
    %% ADAPTIVE: General
    % which ones to run? (i.e. indices of h.Settings.adaptive)
    h.Settings.adaptive_general.adapttypes = [1 2];
    % alternate or randomise runs over types? Alt must have equal number of
    % runs for each adapttype. Cond = one type per CP block
    h.Settings.adaptive_general.seqtype = 'rand'; % 'alt', 'rand', 'cond' 
    h.Settings.adaptive_general.seqtypecond = []; %if 'cond', associate each CP with an adaptive type
    h.Settings.adaptive_general.seqrandblocksize = 12; % should divide the number of trials in a set
    h.Settings.adaptive_general.selectcond.cp = [1:2:7]; % which CP condition to run adaptive on? E.g. AL-only conds.
    h.Settings.adaptive_general.stim = 2; % which stim to run adaptive on?
    h.Settings.adaptive_general.terminate = ''; % terminate within each block only
    h.Settings.adaptive_general.reestimate = 'block'; % 'block' to re-estimate with wider prior each block
    
    %% ADAPTIVE 1
    h.Settings.adaptive(1).type = 'detect';
    h.Settings.adaptive(1).subtype = 'ratio';
    h.Settings.adaptive(1).updown = [1 1];
    % how many of each to run?
    h.Settings.adaptive(1).nRuns = 100*12;
    % max number of thresh estimates to average over to get overall
    % estimate (for plotting only - red line)
    h.Settings.adaptive(1).av_thresh = [];
    h.Settings.adaptive(1).ci_thresh = 20;
    % number of trials each run
    h.Settings.adaptive(1).trialsperrun = 1;
    % adaptive staircase: meanings of the buttonopt
    h.Settings.adaptive(1).buttonfun = {'LeftArrow','RightArrow'};  
    % adaptive staircase: corresponding signal values that would signify a
    % correct answer
    h.Settings.adaptive(1).signalval = [1 2];
    % starting level of adaptive staircase
    h.Settings.adaptive(1).startinglevel = 0; % should be a value in mA. 
    % adapt to omissions of response (not suitable for 2AFC tasks, so set to 0)
    h.Settings.adaptive(1).omit = 0; % 1 = omission is incorrect; 2 = omission is correct
    % which trials (or oddballs if oddonly selected) to start adaptive procedure if there is an omission?
    h.Settings.adaptive(1).startomit = 0;
    % maximum amount - for safety
    h.Settings.adaptive(1).levelmax = 100; % should be value in mA. 
    h.Settings.adaptive(1).levelmin = 0;
    %h.Settings.adaptive(1).expected_change = 5; % smaller value increases the precision of the prior for ZEST and reduces step size of changes in estimates
    h.Settings.adaptive(1).ignoretrials = 5;
    h.Settings.adaptive(1).accuracy_ratio = 3; % for ratio subtype
    h.Settings.adaptive(1).ratio_trials = 4; % number of trials per intensity to use for calculating ratio
    h.Settings.adaptive(1).min_trial_percond_ratio = 2;
    h.Settings.adaptive(1).slope_stimratio = 32; % number to divide stimulus level by, to calculate slope of ZEST
    
    
    %% ADAPTIVE 2
    h.Settings.adaptive(2).type = 'discrim';
    h.Settings.adaptive(2).updown = [1 2];
    % how many of each to run?
    h.Settings.adaptive(2).nRuns = 100*12;
    % max number of thresh estimates to average over to get overall estimate
    h.Settings.adaptive(2).av_thresh = [];
    h.Settings.adaptive(2).ci_thresh = 20;
    % number of trials each run
    h.Settings.adaptive(2).trialsperrun = 1;
    % adaptive staircase: meanings of the buttonopt
    h.Settings.adaptive(2).buttonfun = {'LeftArrow','RightArrow'}; 
    % adaptive staircase: corresponding signal values that would signify a
    % correct answer
    h.Settings.adaptive(2).signalval = [1 2];
    % reversals
    h.Settings.adaptive(2).reversals = [4;8;12];
    % stepsize
    h.Settings.adaptive(2).stepsize = [2;sqrt(2);sqrt(sqrt(2))];
    % steptype 0 = multiple/divide by stepsize; steptype 1 = add/subtract
    h.Settings.adaptive(2).steptype = 0;
    % stepdir -1 = level decreases intensity; stepdir 1 = level increases intensity
    h.Settings.adaptive(2).stepdir = 1;
    % starting level of adaptive staircase
    h.Settings.adaptive(2).startinglevel = 4; % should be a DIFFERENCE value in mA. Keep small as it will increase naturally over time.
    % adapt to omissions of response (not suitable for 2AFC tasks, so set to 0)
    h.Settings.adaptive(2).omit = 0; % 1 = omission is incorrect; 2 = omission is correct
    % which trials (or oddballs if oddonly selected) to start adaptive procedure if there is an omission?
    h.Settings.adaptive(2).startomit = 0;
    % adapt on every trial or only just before an oddball?
    %h.Settings.adaptive.oddonly = 1;
    % max number of trials after oddball that subject must respond (otherwise counts as omitted response)
    %h.Settings.adaptive.resptrials = 4;
    % number of reversals to average over to calculate threshold.
    h.Settings.adaptive(2).reversalForthresh = 6;
    % use mean from the first X responses of each type (high and low)
    %h.Settings.adaptive(2).getmeanfromresponses = 6;
    % maximum amount to adjust the mean if their responses are very
    % incorrect (should be a small fraction, e.g. 1/5th, of the stimulus intensity)
    %h.Settings.adaptive(2).meanadjustmax = 10;
    % maximum amount of the difference value (should be a small fraction, e.g. 1/5th, of the stimulus intensity)
    h.Settings.adaptive(2).levelmax = 100; % should be a DIFFERENCE value in mA.
    h.Settings.adaptive(2).levelmin = 0;
    h.Settings.adaptive(2).maxtrial = inf;
    %h.Settings.adaptive(2).expected_change = 2; % smaller value increases the precision of the prior for ZEST and reduces step size of changes in estimates
    h.Settings.adaptive(2).ignoretrials = 5;
    h.Settings.adaptive(2).eta_divide = 2;
    h.Settings.adaptive(2).slope_stimratio = 2; % number to divide stimulus level by, to calculate slope of ZEST. Should be smaller for thresholding, larger for adjustments during an expt
    
   case 'ALPL_EEG'

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
    % '' = no blocking
    h.Settings.blockopt = 'cond';
    % further options for 'divide':
        % number of blocks (containing multiple conditions)
        %h.Settings.nblocks = 2; % must integer-divide each value in h.Settings.cond_rep_init
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
    h.Settings.trialdur = 2.5:0.5:3.5; % 0.3s auditory cue + 1s anticipation + 0.6s min for ERP
    % Tactile: number of pulses per trial
    h.Settings.nstim_trial = 2; % set to zero to be determined by stimdur
    % Tactile: within-trial frequency (Hz) 
    h.Settings.wait=[1 0]; % one value per nstim 
    
    %% first stimulus: audio
    % duration of stimulus in seconds
    h.Settings.stim(1).dur = [0.15 0.15]; % duration of stimulus in seconds; modified by oddball settings
    h.Settings.stim(1).patternmethod = 'pitch';% Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.stim(1).patternvalue = {[540 500],[460 500]}; % one per stimdur in each cell; one cell per oddball value
    h.Settings.stim(1).durtype = 'reg'; % not needed unless 'rand'
    h.Settings.stim(1).inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).f0 = 500; % pitch
    h.Settings.stim(1).inten_type = 'dB'; % either 'dB' or 'abs'
    h.Settings.stim(1).df = 0;
    h.Settings.stim(1).atten = {'inten_diff',1}; % attenuation level in decibels
    h.Settings.stim(1).attenchan = [1 2]; % apply attenuation (e.g. during thresholding) to these chans
    h.Settings.stim(1).control='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    h.Settings.stim(1).chan = [1 2]; 
    h.Settings.stim(1).nrchannels = 2; % total number of channels, e.g. on sound card
    h.Settings.stim(1).Tukey = 0.25; % Apply Tukey window?
    h.Settings.stim(1).Tukeytype = 2; % 1 = apply to each tone within pattern; 2 = apply to whole pattern

    %% second stimulus: tactile
    % duration of stimulus in seconds
    h.Settings.stim(2).dur = 0; % modified by oddball settings
    h.Settings.stim(2).patternmethod = '';% Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.stim(2).patternvalue = []; % one per stimdur
    h.Settings.stim(2).durtype = 'reg'; % not needed unless 'rand'
    h.Settings.stim(2).inten = []; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(2).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(2).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(2).maxinten = 200; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(2).inten_type = 'abs'; % either 'dB' or 'abs'
    h.Settings.stim(2).control='LJTick-DAQ'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    h.Settings.stim(2).chan = [6]; h.Settings.stim(2).chanforLJ = 1;
    h.Settings.stim(2).nrchannels = 1; % total number of channels, e.g. on sound card
    h.Settings.stim(2).npulses_train = 1; % Tactile: number of pulses in a train; set to zero to be determined by stimdur
    h.Settings.stim(2).p_freq=[];% Tactile: within-train frequency (Hz)
    h.Settings.stim(2).labjack_timer=1; % Use timer to control frequency of labjack outputs? Otherwise uses software timing. must use if there are a large number of pulses (e.g. >4)
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
   %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    
    % Oddball method: intensity, index, pitch, channel
    h.Settings.PL.oddballmethod = 'index'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.PL.oddballvalue = {[1 2], [1 2], [1 2], [2 1], [1 2], [2 1], [1 2], [1 2]}; % values to go into h.Seq.signal. One per oddprob row, or leave blank if determined from GUI
    h.Settings.PL.oddballtype = 'classical'; % options: 'roving', 'classical'

     %% SEQUENCES: PL
    % Change probablity (CP): each condition is in rows
    h.Settings.PL.oddprob = [
        % standard (left) vs oddball (right)
        0.5 0.5
        0.8 0.2
        0.5 0.5
        0.8 0.2
        0.5 0.5
        0.8 0.2
        0.5 0.5
        0.8 0.2
        ];
    % index of cols of oddprob that are standards and oddballs. Can be
    % overridden by h.Settings.oddballvalue if using index
    h.Settings.PL.standardind = 1; 
    h.Settings.PL.oddind = 2; 
    % keep oddball trials apart by at least sep_odd standards
    h.Settings.PL.sep_odd = [0 2 0 2 0 2 0 2];%[0 2 0 2 0 2 0 2 0 2 0 2]; % for each CP condition
    % for sep_odd, which indices of h.Settings.oddballvalue to consider
    % each time? (each list will be considered separately)
    h.Settings.PL.sep_odd_ind = {[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2]};
    h.Settings.PL.sep_odd_tol = [1 1 1 1 1 1 1 1]; % set these to be as high as possible (max 1)
    % for each set, ensure a number of leading standards 
    h.Settings.PL.std_lead = [0 0 0 0 0 0 0 0]; % for each CP condition
    % number of sets to randomise together
    h.Settings.PL.n_set = []; % Leave blank to calculate automatically; or one nunmber per CP condition
    % min number of oddballs within each CP condition
    h.Settings.PL.n_odd = [25 10 25 10 25 10 25 10]; % overrides h.Settings.totdur
    % min number of oddballs per randomised set, per CP
    h.Settings.PL.n_odd_set = [25 10 25 10 25 10 25 10]; % overrides h.Settings.totdur
    % randomise sets?
    h.Settings.PL.rand_set = [0 0 0 0 0 0 0 0]; 
    % condition numbers
    h.Settings.PL.condnum = [
        5 6
        3 4 % low intensity more likely
        5 6
        1 2 % high intensity more likely
        5 6
        1 2 % high intensity more likely
        5 6
        3 4 % low intensity more likely
        ]; 
    
    %% SEQUENCE
    
    % Oddball method: intensity, index, pitch, channel
    h.Settings.AL.oddballmethod = 'index'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.AL.oddballvalue = {[2 1], [2 1], [1 2], [1 2], [2 1], [2 1], [1 2], [1 2]}; % values to go into h.Seq.signal. One per oddprob row, or leave blank if determined from GUI
    h.Settings.AL.oddballtype = 'classical'; % options: 'roving', 'classical'
    % Change probablity (CP): each condition is in rows
    h.Settings.AL.oddprob = [
        % standard (left) vs oddball (right)
        0.8 0.2
        0.8 0.2
        0.8 0.2
        0.8 0.2
        0.8 0.2
        0.8 0.2
        0.8 0.2
        0.8 0.2
        ];
%     % index of row of oddprob that are standards and oddballs. Can be
%     % overridden by h.Settings.oddballvalue if using index
%     h.Settings.standardind = 1; 
%     h.Settings.oddind = 2; 
%     % keep oddball trials apart by at least sep_odd standards
%     h.Settings.sep_odd = [2 0 2 2 0 2 2 0 2]; % for each CP condition
%     % for sep_odd, which indices of h.Settings.oddballvalue to consider
%     % each time? (each list will be considered separately)
%     h.Settings.sep_odd_ind = {[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2]};
%     % for each set, ensure a number of leading standards 
%     h.Settings.std_lead = [0 0 0 0 0 0 0 0 0]; % for each CP condition
%     % number of sets to randomise together
%     h.Settings.n_set = []; % Leave blank to calculate automatically; or one nunmber per CP condition
%     % min number of oddballs within each CP condition
%     h.Settings.n_odd = [4, 4, 4,  8, 8, 8,  12, 12, 12]; % overrides h.Settings.totdur
%     % min number of oddballs per randomised set, per CP
%     h.Settings.n_odd_set = h.Settings.n_odd; % overrides h.Settings.totdur
%     % randomise sets?
%     h.Settings.rand_set = [1 1 1 1 1 1 1 1 1]; 
%     % condition numbers
%     h.Settings.condnum = [
%         1 2
%         3 4
%         5 6
%         1 2
%         3 4
%         5 6
%         1 2
%         3 4
%         5 6
%         ]; 
% index of row of oddprob that are standards and oddballs. Can be
    % overridden by h.Settings.oddballvalue if using index
    h.Settings.AL.standardind = 1; 
    h.Settings.AL.oddind = 2; 
    % keep oddball trials apart by at least sep_odd standards
    h.Settings.AL.sep_odd = [2 2 2 2 2 2 2 2];%[0 2 0 2 0 2 0 2 0 2 0 2]; % for each CP condition
    % for sep_odd, which indices of h.Settings.oddballvalue to consider
    % each time? (each list will be considered separately)
    h.Settings.AL.sep_odd_ind = {[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2]};
    h.Settings.AL.sep_odd_tol = [1 1 1 1 1 1 1 1]; % set these to be as high as possible (max 1)
    % for each set, ensure a number of leading standards 
    h.Settings.AL.std_lead = [0 0 0 0 0 0 0 0]; % for each CP condition
    % number of sets to randomise together
    h.Settings.AL.n_set = []; % Leave blank to calculate automatically; or one nunmber per CP condition
    % min number of oddballs within each CP condition
    h.Settings.AL.n_odd = [10 10 10 10 10 10 10 10]; % overrides h.Settings.totdur
    % min number of oddballs per randomised set, per CP
    h.Settings.AL.n_odd_set = [10 10 10 10 10 10 10 10]; % overrides h.Settings.totdur
    % randomise sets?
    h.Settings.AL.rand_set = [0 0 0 0 0 0 0 0]; 
    % condition numbers
    h.Settings.AL.condnum = [
        
        % block pair 1
        3 4 % pairing 2, AL
        3 4 % pairing 2, ALPL, PL low, AL prior pairing consistent
        
        % block pair 2
        1 2 % pairing 1, AL
        1 2 % pairing 1, ALPL, PL high, AL prior pairing consistent
        
        % block pair 3
        3 4 % pairing 2, AL
        3 4 % pairing 1, ALPL, PL high, AL prior pairing consistent
        
        % block pair 4
        1 2 % pairing 1, AL
        1 2 % pairing 2, ALPL, PL low, AL prior pairing consistent
        
        ]; 
    
    
    %% associative learning experiments
    h.Settings.AL.pairing = {
      % cue 1  cue 2
        [1], [2]; %pairing option 1 (e.g. most porbable in some blocks)
        [2], [1]; %pairing option 2 (e.g. least probable)
        }; % more than one number per pairing/cue: will be randomly assigned to each cue
    % if there is more than one number per pairing/cue, what is their
    % percent representation?
    h.Settings.AL.probstim = {1,1,1,1}; % e.g. low intensity diff vs. high intensity diff
    %h.Settings.assoc.probstim = {[0.2 0.8],[1 0],[0.5 0.5],[0.5 0.5],[0.2 0.8],[1 0]}; % e.g. low intensity diff vs. high intensity diff
    % for each condnum, which pair to use?
    h.Settings.AL.pair = [1 1 2 2];
    % for each stimtype (unique value) within h.Settings.assoc.pairing, 
    % what inten_diff multiplier to use?
    h.Settings.AL.intenstim = [1 -1];
    h.Settings.AL.stimnums = [1 2]; % which two stimulus numbers relate to the two levels of h.Seq.signal?
    h.Settings.AL.stimpart = 1; % which part of the discriminative stim to run adaptive on?
    
    %% RESPONSE PARAMETERS
    % record responses during experiment? 0 or 1
    h.Settings.record_response = 0;
    
    %% THRESHOLDING
    % starting level and step size
    %h.Settings.threshold.startinglevel = 2; % for intensity)
    %h.Settings.threshold.step = 2;
    
   
    
    
end

function h = setgeneral(h)

%% EQUIPMENT CONTROL
% record EEG, NS: netstation, BV: brainvision, 'serial': serial port
% serial port
h.Settings.serial = 'COM1';
h.Settings.record_EEG='labjack_DB15';
%h.Settings.EEGport = 'COM3'; % only needed for 'serial' EEG triggers
h.Settings.EEGMarkPattern = 0; % mark EEG for every change in stimulus pattern (0 = start of trial only)
h.Settings.labjack=1; % Use labjack for controlling any equipment?
h.Settings.labjack_DACport = 4;
h.Settings.DAC_multiply = 0.01; % multiply DAC output by this (e.g. to get mA on DS8R)

% Audio sampling rate 
h.Settings.fs = 96000; % don't change this


    % simulate responses?
    h.Settings.simulate_response = 0;