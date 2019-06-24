function h = Seq_Discrim(h,opt)

switch opt
    
    case 'setoptions'
        
    % settings options
    h.SettingsOptions = {'Discrim_dur','Discrim_isochron_cues','Discrim_isochron','Discrim_random_cues','Discrim_random','Discrim_isochron_freqcond'};

   
case 'Discrim_isochron_freqcond'

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
    h.Settings.nstim_trial = 2; %
    % Which stims that are targets for behavioural responses?
    h.Settings.target_stims = [1 2];
    
    %% first stimulus: audio
    h.Settings.stim(1).patternmethod = 'pitch';% Pattern type method: intensity, pitch. Not supported: channel, duration
    % play with these:
    h.Settings.stim(1).dur = {repmat(1/20,1,20),repmat(1/24,1,24),repmat(1/28,1,28)}; % duration of stimulus in seconds; modified by oddball settings. Value set to zero
    h.Settings.stim(1).stimrandind = [];% index of stimdur to randomise/adapt. 
    h.Settings.stim(1).patternvalue = {repmat([300 0],1,20/2),repmat([300 0],1,24/2),repmat([300 0],1,28/2)}; % one per stimdur in each cell; one cell per oddball value
    h.Settings.stim(1).durtype = 'oddballvalue';%'oddballvalue','sequence_rand'; 
    h.Settings.stim(1).inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).f0 = 300; % pitch
    h.Settings.stim(1).inten_type = 'dB'; % either 'dB' or 'abs'
    h.Settings.stim(1).df = 0;
    h.Settings.stim(1).atten = 0; % attenuation level in decibels
    h.Settings.stim(1).attenchan = [1 2]; % apply attenuation (e.g. during thresholding) to these chans
    h.Settings.stim(1).control='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    h.Settings.stim(1).chan = [1 2]; 
    h.Settings.stim(1).nrchannels = 2; % total number of channels, e.g. on sound card
    h.Settings.stim(1).Tukey = 0.25; % Apply Tukey window?
    h.Settings.stim(1).Tukeytype = 1; % 1 = apply to each tone within pattern; 2 = apply to whole pattern
    
    %% first stimulus: audio pip
    pip.patternmethod = '';% Pattern type method: intensity, pitch. Not supported: channel, duration
    % play with these:
    pip.dur = 0.05; % duration of stimulus in seconds; modified by oddball settings. Value set to zero
    pip.stimrandind = [];% index of stimdur to randomise/adapt. 
    pip.patternvalue = []; % one per stimdur in each cell; one cell per oddball value
    pip.durtype = ''; 
    pip.inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    pip.inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    pip.inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    pip.maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    pip.f0 = 600; % pitch
    pip.inten_type = 'dB'; % either 'dB' or 'abs'
    pip.df = 0;
    pip.atten = 0; % attenuation level in decibels
    pip.attenchan = [1 2]; % apply attenuation (e.g. during thresholding) to these chans
    pip.control='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    pip.chan = [1 2]; 
    pip.nrchannels = 2; % total number of channels, e.g. on sound card
    pip.Tukey = 0.25; % Apply Tukey window?
    pip.Tukeytype = 1; % 1 = apply to each tone within pattern; 2 = apply to whole pattern
    
    % duplicate
    for ns=[1 2]
        h.Settings.stim(ns) = h.Settings.stim(1);
    end
    %h.Settings.stim(2) = pip;
    
    % within-trial frequency (s)
    h.Settings.wait=[2 0]; % one value per nstim 
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, index, pitch, channel
    h.Settings.PL.oddballmethod = 'duration'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.PL.oddballvalue = {[1 2],[2 1],[3 4],[4 3],[5 6],[6 5]}; % values to go into h.Seq.signal. One per oddprob row, or leave blank if determined from GUI
    h.Settings.PL.nonoddballstimvalue = {0 0 0 0 0 0};
    h.Settings.PL.oddballtype = 'classical'; % options: 'roving', 'classical'

    %% SEQUENCE
    % Change probablity (CP): each condition is in rows
    h.Settings.PL.oddprob = [
        1/2 1/2
        1/2 1/2
        1/2 1/2
        1/2 1/2
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
    h.Settings.ntrials = [20 20 20 20 20 20];
    
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
    h.Settings.response_nexttrialmin = 2;
    % when does next trial starts after button press? Empty if programmed
    % ISI
    h.Settings.response_nexttrialwait = 0.6:0.2:1.6;
    
    
    %% ADAPTIVE: General
    % which ones to run? (i.e. indices of h.Settings.adaptive)
    h.Settings.adaptive_general.adapttypes = [1];
    % alternate or randomise runs over types? Alt must have equal number of
    % runs for each adapttype. Cond = one type per CP block
    h.Settings.adaptive_general.seqtype = 'rand'; % 'alt', 'rand', 'cond' 
    h.Settings.adaptive_general.seqtypecond = [1 1]; %if 'cond', associate each CP with an adaptive type
    h.Settings.adaptive_general.seqrandblocksize = 12; % should divide the number of trials in a set
    h.Settings.adaptive_general.selectcond.cp = [2]; % which CP condition to run adaptive on?
    h.Settings.adaptive_general.stim = h.Settings.target_stims; % which stims to apply changes to?
    h.Settings.adaptive_general.stim_judge = h.Settings.target_stims; % which stims do participants judge?
    h.Settings.adaptive_general.terminate = ''; % terminate within each block only
    h.Settings.adaptive_general.reestimate = ''; % 'block' to re-estimate with wider prior each block
    
    %% ADAPTIVE 2
    h.Settings.adaptive(1).type = 'discrim';
    h.Settings.adaptive(1).updown = [1 2];
    % how many of each to run?
    h.Settings.adaptive(1).nRuns = 100*12;
    % max number of thresh estimates to average over to get overall estimate
    h.Settings.adaptive(1).av_thresh = [];
    h.Settings.adaptive(1).ci_thresh = 20;
    % number of trials each run
    h.Settings.adaptive(1).trialsperrun = 1;
    % adaptive staircase: meanings of the buttonopt
    h.Settings.adaptive(1).buttonfun = {'LeftArrow','RightArrow'}; 
    % adaptive staircase: corresponding signal values that would signify a
    % correct answer
    h.Settings.adaptive(1).signalval = [1 2]; % first: correct response. Second: incorrect
    h.Settings.adaptive(1).signaltarget = {[1 2 3 4 5 6],[1 2 1 2 1 2]}; % translate actual signal values to signalval indices
    h.Settings.adaptive(1).response_type = '2AFC';%'samediff'; % discrimination is same/difference
    % reversals
    h.Settings.adaptive(1).reversals = [4;8;12];
    % stepsize
    h.Settings.adaptive(1).stepsize = [2;sqrt(2);sqrt(sqrt(2))];
    % steptype 0 = multiple/divide by stepsize; steptype 1 = add/subtract
    h.Settings.adaptive(1).steptype = 0;
    % stepdir -1 = level decreases intensity; stepdir 1 = level increases intensity
    h.Settings.adaptive(1).stepdir = -1;
    % starting level of adaptive staircase
    h.Settings.adaptive(1).startinglevel = 1/(25*10); % should be a DIFFERENCE value. Keep small as it will increase naturally over time.
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
    % maximum amount of the difference value (should be a small fraction, e.g. 1/5th, of the stimulus intensity)
    h.Settings.adaptive(1).levelmax = 0.05; % should be a DIFFERENCE value.
    h.Settings.adaptive(1).levelmin = 0;
    h.Settings.adaptive(1).maxtrial = inf;
    %h.Settings.adaptive(1).expected_change = 5; % smaller value increases the precision of the prior for ZEST and reduces step size of changes in estimates
    h.Settings.adaptive(1).ignoretrials = 0;
    h.Settings.adaptive(1).eta_divide = 2;
    h.Settings.adaptive(1).slope_stimratio = 1; % number to divide stimulus level by, to calculate slope of ZEST. Should be smaller for thresholding, larger for adjustments during an expt
    
    
case 'Discrim_random_cues'

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
    % number of pulses per trial
    h.Settings.nstim_trial = 5; %
    % Which stims that are targets for behavioural responses?
    h.Settings.target_stims = [4 5];
    
    %% first stimulus: audio
    h.Settings.stim(1).patternmethod = 'pitch';% Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.stim(1).gap_ind=2:2:16;
    h.Settings.stim(1).ngap_rand_oddball=5;
    h.Settings.stim(1).dur = [0.05*ones(1,8) 0.01:0.02:0.15]; % duration of stimulus in seconds; modified by oddball settings. Value set to zero
    h.Settings.stim(1).dur([h.Settings.stim(1).gap_ind-1 h.Settings.stim(1).gap_ind]) = h.Settings.stim(1).dur;
    h.Settings.stim(1).stimrandind = h.Settings.stim(1).gap_ind(randperm(length(h.Settings.stim(1).gap_ind)));% index of stimdur to randomise/adapt. 
    h.Settings.stim(1).stimrandind_oddball = h.Settings.stim(1).stimrandind(1:h.Settings.stim(1).ngap_rand_oddball); % indices to further randomise as oddballs
    h.Settings.stim(1).patternvalue = repmat([300 0],1,8); % one per stimdur in each cell; one cell per oddball value
    h.Settings.stim(1).durtype = '';%'sequence_rand'; 
    h.Settings.stim(1).inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).f0 = 300; % pitch
    h.Settings.stim(1).inten_type = 'dB'; % either 'dB' or 'abs'
    h.Settings.stim(1).df = 0;
    h.Settings.stim(1).atten = 0; % attenuation level in decibels
    h.Settings.stim(1).attenchan = [1 2]; % apply attenuation (e.g. during thresholding) to these chans
    h.Settings.stim(1).control='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    h.Settings.stim(1).chan = [1 2]; 
    h.Settings.stim(1).nrchannels = 2; % total number of channels, e.g. on sound card
    h.Settings.stim(1).Tukey = 0.25; % Apply Tukey window?
    h.Settings.stim(1).Tukeytype = 1; % 1 = apply to each tone within pattern; 2 = apply to whole pattern
    
    %% first stimulus: audio pip
    pip.patternmethod = '';% Pattern type method: intensity, pitch. Not supported: channel, duration
    % play with these:
    pip.gap_ind=[];
    pip.ngap_rand_oddball=0;
    pip.dur = 0.05; % duration of stimulus in seconds; modified by oddball settings. Value set to zero
    pip.stimrandind = [];% index of stimdur to randomise/adapt. 
    pip.stimrandind_oddball = []; % indices to further randomise as oddballs
    pip.patternvalue = []; % one per stimdur in each cell; one cell per oddball value
    pip.durtype = ''; 
    pip.inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    pip.inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    pip.inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    pip.maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    pip.f0 = 600; % pitch
    pip.inten_type = 'dB'; % either 'dB' or 'abs'
    pip.df = 0;
    pip.atten = 0; % attenuation level in decibels
    pip.attenchan = [1 2]; % apply attenuation (e.g. during thresholding) to these chans
    pip.control='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    pip.chan = [1 2]; 
    pip.nrchannels = 2; % total number of channels, e.g. on sound card
    pip.Tukey = 0.25; % Apply Tukey window?
    pip.Tukeytype = 1; % 1 = apply to each tone within pattern; 2 = apply to whole pattern
    
    % duplicate
    for ns=[1 2 4 5]
        h.Settings.stim(ns) = h.Settings.stim(1);
    end
    h.Settings.stim(3) = pip;
    
    % within-trial frequency (s)
    h.Settings.wait=[sum(h.Settings.stim(1).dur)+1 sum(h.Settings.stim(2).dur)+0.5 sum(h.Settings.stim(3).dur)+0.5 sum(h.Settings.stim(4).dur)+1 0]; % one value per nstim 
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, index, pitch, channel
    h.Settings.PL.oddballmethod = 'duration_shuffle'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.PL.oddballvalue = {[1 1],[1 2],[2 1],[2 2],[1 1],[1 2],[2 1],[2 2]}; % values to go into h.Seq.signal. One per oddprob row, or leave blank if determined from GUI
    h.Settings.PL.nonoddballstimvalue = {1,1,1,1,2,2,2,2};
    h.Settings.PL.oddballtype = 'classical'; % options: 'roving', 'classical'

    %% SEQUENCE
    % Change probablity (CP): each condition is in rows
    h.Settings.PL.oddprob = [
        1/2 1/2
        1/2 1/2
        1/2 1/2
        1/2 1/2
        1/2 1/2
        1/2 1/2
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
    h.Settings.ntrials = [10 10 10 10 10 10 10 10];
    
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
    h.Settings.response_nexttrialmin = 2;
    % when does next trial starts after button press? Empty if programmed
    % ISI
    h.Settings.response_nexttrialwait = 0.6:0.2:1.6;
    
    
    %% ADAPTIVE: General
    % which ones to run? (i.e. indices of h.Settings.adaptive)
    h.Settings.adaptive_general.adapttypes = [1];
    % alternate or randomise runs over types? Alt must have equal number of
    % runs for each adapttype. Cond = one type per CP block
    h.Settings.adaptive_general.seqtype = 'rand'; % 'alt', 'rand', 'cond' 
    h.Settings.adaptive_general.seqtypecond = [1 1]; %if 'cond', associate each CP with an adaptive type
    h.Settings.adaptive_general.seqrandblocksize = 12; % should divide the number of trials in a set
    h.Settings.adaptive_general.selectcond.cp = [2]; % which CP condition to run adaptive on?
    h.Settings.adaptive_general.stim = [1 2 4 5]; % which stims to apply changes to?
    h.Settings.adaptive_general.stim_judge = h.Settings.target_stims; % which stims to participants judge?
    h.Settings.adaptive_general.terminate = ''; % terminate within each block only
    h.Settings.adaptive_general.reestimate = ''; % 'block' to re-estimate with wider prior each block
    
    %% ADAPTIVE 2
    h.Settings.adaptive(1).type = 'discrim';
    h.Settings.adaptive(1).updown = [1 2];
    % how many of each to run?
    h.Settings.adaptive(1).nRuns = 100*12;
    % max number of thresh estimates to average over to get overall estimate
    h.Settings.adaptive(1).av_thresh = [];
    h.Settings.adaptive(1).ci_thresh = 20;
    % number of trials each run
    h.Settings.adaptive(1).trialsperrun = 1;
    % adaptive staircase: meanings of the buttonopt
    h.Settings.adaptive(1).buttonfun = {'LeftArrow','RightArrow'}; 
    % adaptive staircase: corresponding signal values that would signify a
    % correct answer
    h.Settings.adaptive(1).signalval = [1 2];
    h.Settings.adaptive(1).response_type = 'samediff'; % discrimination is same/difference
    % reversals
    h.Settings.adaptive(1).reversals = [4;8;12];
    % stepsize
    h.Settings.adaptive(1).stepsize = [2;sqrt(2);sqrt(sqrt(2))];
    % steptype 0 = multiple/divide by stepsize; steptype 1 = add/subtract
    h.Settings.adaptive(1).steptype = 0;
    % stepdir -1 = level decreases intensity; stepdir 1 = level increases intensity
    h.Settings.adaptive(1).stepdir = -1;
    % starting level of adaptive staircase
    h.Settings.adaptive(1).startinglevel = h.Settings.stim(h.Settings.adaptive_general.stim(1)).ngap_rand_oddball; % should be a DIFFERENCE value. Keep small as it will increase naturally over time.
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
    % maximum amount of the difference value (should be a small fraction, e.g. 1/5th, of the stimulus intensity)
    h.Settings.adaptive(1).levelmax = length(h.Settings.stim(h.Settings.adaptive_general.stim(1)).gap_ind); % should be a DIFFERENCE value.
    h.Settings.adaptive(1).levelmin = 0;
    h.Settings.adaptive(1).maxtrial = inf;
    %h.Settings.adaptive(1).expected_change = 5; % smaller value increases the precision of the prior for ZEST and reduces step size of changes in estimates
    h.Settings.adaptive(1).ignoretrials = 5;
    h.Settings.adaptive(1).eta_divide = 2;
    h.Settings.adaptive(1).slope_stimratio = 1; % number to divide stimulus level by, to calculate slope of ZEST. Should be smaller for thresholding, larger for adjustments during an expt
    

case 'Discrim_random'

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
    % number of pulses per trial
    h.Settings.nstim_trial = 2; %
    % Which stims that are targets for behavioural responses?
    h.Settings.target_stims = [1 2];
    
    %% first stimulus: audio
    h.Settings.stim(1).patternmethod = 'pitch';% Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.stim(1).gap_ind=2:2:16;
    h.Settings.stim(1).ngap_rand_oddball=5;
    h.Settings.stim(1).dur = [0.05*ones(1,8) 0.01:0.02:0.15]; % duration of stimulus in seconds; modified by oddball settings. Value set to zero
    h.Settings.stim(1).dur([h.Settings.stim(1).gap_ind-1 h.Settings.stim(1).gap_ind]) = h.Settings.stim(1).dur;
    h.Settings.stim(1).stimrandind = h.Settings.stim(1).gap_ind(randperm(length(h.Settings.stim(1).gap_ind)));% index of stimdur to randomise/adapt. 
    h.Settings.stim(1).stimrandind_oddball = h.Settings.stim(1).stimrandind(1:h.Settings.stim(1).ngap_rand_oddball); % indices to further randomise as oddballs
    h.Settings.stim(1).patternvalue = repmat([300 0],1,8); % one per stimdur in each cell; one cell per oddball value
    h.Settings.stim(1).durtype = '';%'sequence_rand'; % randomises patternmethod during sequence setup
    h.Settings.stim(1).inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).f0 = 300; % pitch
    h.Settings.stim(1).inten_type = 'dB'; % either 'dB' or 'abs'
    h.Settings.stim(1).df = 0;
    h.Settings.stim(1).atten = 0; % attenuation level in decibels
    h.Settings.stim(1).attenchan = [1 2]; % apply attenuation (e.g. during thresholding) to these chans
    h.Settings.stim(1).control='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    h.Settings.stim(1).chan = [1 2]; 
    h.Settings.stim(1).nrchannels = 2; % total number of channels, e.g. on sound card
    h.Settings.stim(1).Tukey = 0.25; % Apply Tukey window?
    h.Settings.stim(1).Tukeytype = 1; % 1 = apply to each tone within pattern; 2 = apply to whole pattern
    
    %% first stimulus: audio pip
    pip.patternmethod = '';% Pattern type method: intensity, pitch. Not supported: channel, duration
    % play with these:
    pip.gap_ind=[];
    pip.ngap_rand_oddball=0;
    pip.dur = 0.05; % duration of stimulus in seconds; modified by oddball settings. Value set to zero
    pip.stimrandind = [];% index of stimdur to randomise/adapt. 
    pip.stimrandind_oddball = []; % indices to further randomise as oddballs
    pip.patternvalue = []; % one per stimdur in each cell; one cell per oddball value
    pip.durtype = ''; 
    pip.inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    pip.inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    pip.inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    pip.maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    pip.f0 = 600; % pitch
    pip.inten_type = 'dB'; % either 'dB' or 'abs'
    pip.df = 0;
    pip.atten = 0; % attenuation level in decibels
    pip.attenchan = [1 2]; % apply attenuation (e.g. during thresholding) to these chans
    pip.control='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    pip.chan = [1 2]; 
    pip.nrchannels = 2; % total number of channels, e.g. on sound card
    pip.Tukey = 0.25; % Apply Tukey window?
    pip.Tukeytype = 1; % 1 = apply to each tone within pattern; 2 = apply to whole pattern
    
    % duplicate
    for ns=[1 2]
        h.Settings.stim(ns) = h.Settings.stim(1);
    end
    %h.Settings.stim(2) = pip;
    
    % within-trial frequency (s)
    h.Settings.wait=[sum(h.Settings.stim(1).dur)+1 0]; % one value per nstim 
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, index, pitch, channel
    h.Settings.PL.oddballmethod = 'duration_shuffle'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.PL.oddballvalue = {[1 1],[1 2],[2 1],[2 2]}; % values to go into h.Seq.signal. One per oddprob row, or leave blank if determined from GUI
    h.Settings.PL.nonoddballstimvalue = {0 0 0 0};
    h.Settings.PL.oddballtype = 'classical'; % options: 'roving', 'classical'

    %% SEQUENCE
    % Change probablity (CP): each condition is in rows
    h.Settings.PL.oddprob = [
        1/2 1/2
        1/2 1/2
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
    h.Settings.ntrials = [10 10 10 10 10 10 10 10];
    
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
    h.Settings.response_nexttrialmin = 2;
    % when does next trial starts after button press? Empty if programmed
    % ISI
    h.Settings.response_nexttrialwait = 0.6:0.2:1.6;
    
    
    %% ADAPTIVE: General
    % which ones to run? (i.e. indices of h.Settings.adaptive)
    h.Settings.adaptive_general.adapttypes = [1];
    % alternate or randomise runs over types? Alt must have equal number of
    % runs for each adapttype. Cond = one type per CP block
    h.Settings.adaptive_general.seqtype = 'rand'; % 'alt', 'rand', 'cond' 
    h.Settings.adaptive_general.seqtypecond = [1 1]; %if 'cond', associate each CP with an adaptive type
    h.Settings.adaptive_general.seqrandblocksize = 12; % should divide the number of trials in a set
    h.Settings.adaptive_general.selectcond.cp = [2]; % which CP condition to run adaptive on?
    h.Settings.adaptive_general.stim = h.Settings.target_stims; % which stims to apply changes to?
    h.Settings.adaptive_general.stim_judge = h.Settings.target_stims; % which stims to participants judge?
    h.Settings.adaptive_general.terminate = 'block'; % terminate within each block only
    h.Settings.adaptive_general.reestimate = ''; % 'block' to re-estimate with wider prior each block
    
    %% ADAPTIVE 2
    h.Settings.adaptive(1).type = 'discrim';
    h.Settings.adaptive(1).updown = [1 2];
    % how many of each to run?
    h.Settings.adaptive(1).nRuns = 100*12;
    % max number of thresh estimates to average over to get overall estimate
    h.Settings.adaptive(1).av_thresh = [];
    h.Settings.adaptive(1).ci_thresh = 20;
    % number of trials each run
    h.Settings.adaptive(1).trialsperrun = 1;
    % adaptive staircase: meanings of the buttonopt
    h.Settings.adaptive(1).buttonfun = {'LeftArrow','RightArrow'}; 
    % adaptive staircase: corresponding signal values that would signify a
    % correct answer
    h.Settings.adaptive(1).signalval = [1 2];
    h.Settings.adaptive(1).response_type = 'samediff'; % discrimination is same/difference
    % reversals
    h.Settings.adaptive(1).reversals = [4;8;12];
    % stepsize
    h.Settings.adaptive(1).stepsize = [2;sqrt(2);sqrt(sqrt(2))];
    % steptype 0 = multiple/divide by stepsize; steptype 1 = add/subtract
    h.Settings.adaptive(1).steptype = 0;
    % stepdir -1 = level decreases intensity; stepdir 1 = level increases intensity
    h.Settings.adaptive(1).stepdir = -1;
    % starting level of adaptive staircase
    h.Settings.adaptive(1).startinglevel = h.Settings.stim(h.Settings.adaptive_general.stim(1)).ngap_rand_oddball; % should be a DIFFERENCE value. Keep small as it will increase naturally over time.
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
    % maximum amount of the difference value (should be a small fraction, e.g. 1/5th, of the stimulus intensity)
    h.Settings.adaptive(1).levelmax = length(h.Settings.stim(h.Settings.adaptive_general.stim(1)).gap_ind); % should be a DIFFERENCE value.
    h.Settings.adaptive(1).levelmin = 0;
    h.Settings.adaptive(1).maxtrial = inf;
    %h.Settings.adaptive(1).expected_change = 5; % smaller value increases the precision of the prior for ZEST and reduces step size of changes in estimates
    h.Settings.adaptive(1).ignoretrials = 0;
    h.Settings.adaptive(1).eta_divide = 2;
    h.Settings.adaptive(1).slope_stimratio = 1; % number to divide stimulus level by, to calculate slope of ZEST. Should be smaller for thresholding, larger for adjustments during an expt
    
case 'Discrim_isochron_cues'

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
    h.Settings.nstim_trial = 5; %
    % Which stims that are targets for behavioural responses?
    h.Settings.target_stims = [4 5];
    
    %% first stimulus: audio
    h.Settings.stim(1).patternmethod = 'pitch';% Pattern type method: intensity, pitch. Not supported: channel, duration
    % play with these:
    h.Settings.stim(1).dur = repmat(0.05,1,20); % duration of stimulus in seconds; modified by oddball settings. Value set to zero
    h.Settings.stim(1).stimrandind = [];% index of stimdur to randomise/adapt. 
    h.Settings.stim(1).patternvalue = repmat([300 0],1,10); % one per stimdur in each cell; one cell per oddball value
    h.Settings.stim(1).durtype = '';% 'sequence_rand'; 
    h.Settings.stim(1).inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).f0 = 300; % pitch
    h.Settings.stim(1).inten_type = 'dB'; % either 'dB' or 'abs'
    h.Settings.stim(1).df = 0;
    h.Settings.stim(1).atten = 0; % attenuation level in decibels
    h.Settings.stim(1).attenchan = [1 2]; % apply attenuation (e.g. during thresholding) to these chans
    h.Settings.stim(1).control='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    h.Settings.stim(1).chan = [1 2]; 
    h.Settings.stim(1).nrchannels = 2; % total number of channels, e.g. on sound card
    h.Settings.stim(1).Tukey = 0.25; % Apply Tukey window?
    h.Settings.stim(1).Tukeytype = 1; % 1 = apply to each tone within pattern; 2 = apply to whole pattern
    
    %% first stimulus: audio pip
    pip.patternmethod = '';% Pattern type method: intensity, pitch. Not supported: channel, duration
    % play with these:
    pip.dur = 0.05; % duration of stimulus in seconds; modified by oddball settings. Value set to zero
    pip.stimrandind = [];% index of stimdur to randomise/adapt. 
    pip.patternvalue = []; % one per stimdur in each cell; one cell per oddball value
    pip.durtype = ''; 
    pip.inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    pip.inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    pip.inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    pip.maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    pip.f0 = 600; % pitch
    pip.inten_type = 'dB'; % either 'dB' or 'abs'
    pip.df = 0;
    pip.atten = 0; % attenuation level in decibels
    pip.attenchan = [1 2]; % apply attenuation (e.g. during thresholding) to these chans
    pip.control='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    pip.chan = [1 2]; 
    pip.nrchannels = 2; % total number of channels, e.g. on sound card
    pip.Tukey = 0.25; % Apply Tukey window?
    pip.Tukeytype = 1; % 1 = apply to each tone within pattern; 2 = apply to whole pattern
    
    % duplicate
    for ns=[1 2 4 5]
        h.Settings.stim(ns) = h.Settings.stim(1);
    end
    h.Settings.stim(3) = pip;
    
    % within-trial frequency (s)
    h.Settings.wait=[sum(h.Settings.stim(1).dur)+1 sum(h.Settings.stim(2).dur)+0.5 sum(h.Settings.stim(3).dur)+0.5 sum(h.Settings.stim(4).dur)+1 0]; % one value per nstim 
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, index, pitch, channel
    h.Settings.PL.oddballmethod = 'duration'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.PL.oddballvalue = {[1 1],[1 2],[2 1],[2 2],[1 1],[1 2],[2 1],[2 2]}; % values to go into h.Seq.signal. One per oddprob row, or leave blank if determined from GUI
    h.Settings.PL.nonoddballstimvalue = {1,1,1,1,2,2,2,2};
    h.Settings.PL.oddballtype = 'classical'; % options: 'roving', 'classical' 

    %% SEQUENCE
    % Change probablity (CP): each condition is in rows
    h.Settings.PL.oddprob = [
        1/2 1/2
        1/2 1/2
        1/2 1/2
        1/2 1/2
        1/2 1/2
        1/2 1/2
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
    h.Settings.ntrials = [10 10 10 10 10 10 10 10];
    
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
    h.Settings.response_nexttrialmin = 2;
    % when does next trial starts after button press? Empty if programmed
    % ISI
    h.Settings.response_nexttrialwait = 0.6:0.2:1.6;
    
    
    %% ADAPTIVE: General
    % which ones to run? (i.e. indices of h.Settings.adaptive)
    h.Settings.adaptive_general.adapttypes = [1];
    % alternate or randomise runs over types? Alt must have equal number of
    % runs for each adapttype. Cond = one type per CP block
    h.Settings.adaptive_general.seqtype = 'rand'; % 'alt', 'rand', 'cond' 
    h.Settings.adaptive_general.seqtypecond = [1 1]; %if 'cond', associate each CP with an adaptive type
    h.Settings.adaptive_general.seqrandblocksize = 12; % should divide the number of trials in a set
    h.Settings.adaptive_general.selectcond.cp = [2]; % which CP condition to run adaptive on?
    h.Settings.adaptive_general.stim = h.Settings.target_stims; % which stims to apply changes to?
    h.Settings.adaptive_general.stim_judge = h.Settings.target_stims; % which stims to participants judge?
    h.Settings.adaptive_general.terminate = ''; % terminate within each block only
    h.Settings.adaptive_general.reestimate = ''; % 'block' to re-estimate with wider prior each block
    
    %% ADAPTIVE 2
    h.Settings.adaptive(1).type = 'discrim';
    h.Settings.adaptive(1).updown = [1 2];
    % how many of each to run?
    h.Settings.adaptive(1).nRuns = 100*12;
    % max number of thresh estimates to average over to get overall estimate
    h.Settings.adaptive(1).av_thresh = [];
    h.Settings.adaptive(1).ci_thresh = 20;
    % number of trials each run
    h.Settings.adaptive(1).trialsperrun = 1;
    % adaptive staircase: meanings of the buttonopt
    h.Settings.adaptive(1).buttonfun = {'LeftArrow','RightArrow'}; 
    % adaptive staircase: corresponding signal values that would signify a
    % correct answer
    h.Settings.adaptive(1).signalval = [1 2];
    h.Settings.adaptive(1).response_type = 'samediff'; % discrimination is same/difference
    % reversals
    h.Settings.adaptive(1).reversals = [4;8;12];
    % stepsize
    h.Settings.adaptive(1).stepsize = [2;sqrt(2);sqrt(sqrt(2))];
    % steptype 0 = multiple/divide by stepsize; steptype 1 = add/subtract
    h.Settings.adaptive(1).steptype = 0;
    % stepdir -1 = level decreases intensity; stepdir 1 = level increases intensity
    h.Settings.adaptive(1).stepdir = -1;
    % starting level of adaptive staircase
    h.Settings.adaptive(1).startinglevel = 0.005; % should be a DIFFERENCE value. Keep small as it will increase naturally over time.
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
    % maximum amount of the difference value (should be a small fraction, e.g. 1/5th, of the stimulus intensity)
    h.Settings.adaptive(1).levelmax = 0.05; % should be a DIFFERENCE value.
    h.Settings.adaptive(1).levelmin = 0;
    h.Settings.adaptive(1).maxtrial = inf;
    %h.Settings.adaptive(1).expected_change = 5; % smaller value increases the precision of the prior for ZEST and reduces step size of changes in estimates
    h.Settings.adaptive(1).ignoretrials = 5;
    h.Settings.adaptive(1).eta_divide = 2;
    h.Settings.adaptive(1).slope_stimratio = 1; % number to divide stimulus level by, to calculate slope of ZEST. Should be smaller for thresholding, larger for adjustments during an expt
    
case 'Discrim_isochron'

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
    h.Settings.nstim_trial = 2; %
    % Which stims that are targets for behavioural responses?
    h.Settings.target_stims = [1 2];
    
    %% first stimulus: audio
    h.Settings.stim(1).patternmethod = 'pitch';% Pattern type method: intensity, pitch. Not supported: channel, duration
    % play with these:
    h.Settings.stim(1).dur = repmat(0.05,1,20); % duration of stimulus in seconds; modified by oddball settings. Value set to zero
    h.Settings.stim(1).stimrandind = [];% index of stimdur to randomise/adapt. 
    h.Settings.stim(1).patternvalue = repmat([300 0],1,10); % one per stimdur in each cell; one cell per oddball value
    h.Settings.stim(1).durtype = '';%'sequence_rand'; 
    h.Settings.stim(1).inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).f0 = 300; % pitch
    h.Settings.stim(1).inten_type = 'dB'; % either 'dB' or 'abs'
    h.Settings.stim(1).df = 0;
    h.Settings.stim(1).atten = 0; % attenuation level in decibels
    h.Settings.stim(1).attenchan = [1 2]; % apply attenuation (e.g. during thresholding) to these chans
    h.Settings.stim(1).control='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    h.Settings.stim(1).chan = [1 2]; 
    h.Settings.stim(1).nrchannels = 2; % total number of channels, e.g. on sound card
    h.Settings.stim(1).Tukey = 0.25; % Apply Tukey window?
    h.Settings.stim(1).Tukeytype = 1; % 1 = apply to each tone within pattern; 2 = apply to whole pattern
    
    %% first stimulus: audio pip
    pip.patternmethod = '';% Pattern type method: intensity, pitch. Not supported: channel, duration
    % play with these:
    pip.dur = 0.05; % duration of stimulus in seconds; modified by oddball settings. Value set to zero
    pip.stimrandind = [];% index of stimdur to randomise/adapt. 
    pip.patternvalue = []; % one per stimdur in each cell; one cell per oddball value
    pip.durtype = ''; 
    pip.inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    pip.inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    pip.inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    pip.maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    pip.f0 = 600; % pitch
    pip.inten_type = 'dB'; % either 'dB' or 'abs'
    pip.df = 0;
    pip.atten = 0; % attenuation level in decibels
    pip.attenchan = [1 2]; % apply attenuation (e.g. during thresholding) to these chans
    pip.control='PsychPortAudio'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt
    pip.chan = [1 2]; 
    pip.nrchannels = 2; % total number of channels, e.g. on sound card
    pip.Tukey = 0.25; % Apply Tukey window?
    pip.Tukeytype = 1; % 1 = apply to each tone within pattern; 2 = apply to whole pattern
    
    % duplicate
    for ns=[1 2]
        h.Settings.stim(ns) = h.Settings.stim(1);
    end
    %h.Settings.stim(2) = pip;
    
    % within-trial frequency (s)
    h.Settings.wait=[sum(h.Settings.stim(1).dur)+1 0]; % one value per nstim 
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, index, pitch, channel
    h.Settings.PL.oddballmethod = 'duration'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.PL.oddballvalue = {[1 1],[1 2],[2 1],[2 2]}; % values to go into h.Seq.signal. One per oddprob row, or leave blank if determined from GUI
    h.Settings.PL.nonoddballstimvalue = {0 0 0 0};
    h.Settings.PL.oddballtype = 'classical'; % options: 'roving', 'classical'

    %% SEQUENCE
    % Change probablity (CP): each condition is in rows
    h.Settings.PL.oddprob = [
        1/2 1/2
        1/2 1/2
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
    h.Settings.ntrials = [10 10 10 10 10 10 10 10];
    
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
    h.Settings.response_nexttrialmin = 2;
    % when does next trial starts after button press? Empty if programmed
    % ISI
    h.Settings.response_nexttrialwait = 0.6:0.2:1.6;
    
    
    %% ADAPTIVE: General
    % which ones to run? (i.e. indices of h.Settings.adaptive)
    h.Settings.adaptive_general.adapttypes = [1];
    % alternate or randomise runs over types? Alt must have equal number of
    % runs for each adapttype. Cond = one type per CP block
    h.Settings.adaptive_general.seqtype = 'rand'; % 'alt', 'rand', 'cond' 
    h.Settings.adaptive_general.seqtypecond = [1 1]; %if 'cond', associate each CP with an adaptive type
    h.Settings.adaptive_general.seqrandblocksize = 12; % should divide the number of trials in a set
    h.Settings.adaptive_general.selectcond.cp = [2]; % which CP condition to run adaptive on?
    h.Settings.adaptive_general.stim = h.Settings.target_stims; % which stims to apply changes to?
    h.Settings.adaptive_general.stim_judge = h.Settings.target_stims; % which stims to participants judge?
    h.Settings.adaptive_general.terminate = 'block'; % terminate within each block only
    h.Settings.adaptive_general.reestimate = ''; % 'block' to re-estimate with wider prior each block
    
    %% ADAPTIVE 2
    h.Settings.adaptive(1).type = 'discrim';
    h.Settings.adaptive(1).updown = [1 2];
    % how many of each to run?
    h.Settings.adaptive(1).nRuns = 100*12;
    % max number of thresh estimates to average over to get overall estimate
    h.Settings.adaptive(1).av_thresh = [];
    h.Settings.adaptive(1).ci_thresh = 20;
    % number of trials each run
    h.Settings.adaptive(1).trialsperrun = 1;
    % adaptive staircase: meanings of the buttonopt
    h.Settings.adaptive(1).buttonfun = {'LeftArrow','RightArrow'}; 
    % adaptive staircase: corresponding signal values that would signify a
    % correct answer
    h.Settings.adaptive(1).signalval = [1 2];
    h.Settings.adaptive(1).response_type = 'samediff'; % discrimination is same/difference
    % reversals
    h.Settings.adaptive(1).reversals = [4;8;12];
    % stepsize
    h.Settings.adaptive(1).stepsize = [2;sqrt(2);sqrt(sqrt(2))];
    % steptype 0 = multiple/divide by stepsize; steptype 1 = add/subtract
    h.Settings.adaptive(1).steptype = 0;
    % stepdir -1 = level decreases intensity; stepdir 1 = level increases intensity
    h.Settings.adaptive(1).stepdir = -1;
    % starting level of adaptive staircase
    h.Settings.adaptive(1).startinglevel = 0.005; % should be a DIFFERENCE value. Keep small as it will increase naturally over time.
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
    % maximum amount of the difference value (should be a small fraction, e.g. 1/5th, of the stimulus intensity)
    h.Settings.adaptive(1).levelmax = 0.05; % should be a DIFFERENCE value.
    h.Settings.adaptive(1).levelmin = 0;
    h.Settings.adaptive(1).maxtrial = inf;
    %h.Settings.adaptive(1).expected_change = 5; % smaller value increases the precision of the prior for ZEST and reduces step size of changes in estimates
    h.Settings.adaptive(1).ignoretrials = 0;
    h.Settings.adaptive(1).eta_divide = 2;
    h.Settings.adaptive(1).slope_stimratio = 1; % number to divide stimulus level by, to calculate slope of ZEST. Should be smaller for thresholding, larger for adjustments during an expt
    
     
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
    h.Settings.stim(1).dur = repmat(0.05,1,20); % duration of stimulus in seconds; modified by oddball settings. Value set to zero
    h.Settings.stim(1).stimrandind = [];% index of stimdur to randomise/adapt. 
    h.Settings.stim(1).patternvalue = repmat([300 0],1,10); % one per stimdur in each cell; one cell per oddball value
    
    h.Settings.stim(1).durtype = 'sequence_rand'; 
    h.Settings.stim(1).inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).f0 = 300; % pitch
    h.Settings.stim(1).inten_type = 'dB'; % either 'dB' or 'abs'
    h.Settings.stim(1).df = 0;
    h.Settings.stim(1).atten = 0; % attenuation level in decibels
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
    h.Settings.PL.oddballmethod = 'duration'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.PL.oddballvalue = {[1 1],[1 2]}; % values to go into h.Seq.signal. One per oddprob row, or leave blank if determined from GUI
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
    h.Settings.buttonopt = {'LeftArrow','RightArrow'}; 
    % how early after start of trial can button press trigger the next trial? Empty if programmed
    % ISI
    h.Settings.response_nexttrialmin = 2;
    % when does next trial starts after button press? Empty if programmed
    % ISI
    h.Settings.response_nexttrialwait = 0.6:0.2:1.6;
     
    %% ADAPTIVE: General
    % which ones to run? (i.e. indices of h.Settings.adaptive)
    h.Settings.adaptive_general.adapttypes = [1];
    % alternate or randomise runs over types? Alt must have equal number of
    % runs for each adapttype. Cond = one type per CP block
    h.Settings.adaptive_general.seqtype = 'rand'; % 'alt', 'rand', 'cond' 
    h.Settings.adaptive_general.seqtypecond = [1 1]; %if 'cond', associate each CP with an adaptive type
    h.Settings.adaptive_general.seqrandblocksize = 12; % should divide the number of trials in a set
    h.Settings.adaptive_general.selectcond.cp = [2]; % which CP condition to run adaptive on?
    h.Settings.adaptive_general.stim = [1 2]; % which stim to run adaptive on?
    h.Settings.adaptive_general.terminate = 'block'; % terminate within each block only
    h.Settings.adaptive_general.reestimate = ''; % 'block' to re-estimate with wider prior each block
    
    %% ADAPTIVE 2
    h.Settings.adaptive(1).type = 'discrim';
    h.Settings.adaptive(1).updown = [1 2];
    % how many of each to run?
    h.Settings.adaptive(1).nRuns = 100*12;
    % max number of thresh estimates to average over to get overall estimate
    h.Settings.adaptive(1).av_thresh = [];
    h.Settings.adaptive(1).ci_thresh = 20;
    % number of trials each run
    h.Settings.adaptive(1).trialsperrun = 1;
    % adaptive staircase: meanings of the buttonopt
    h.Settings.adaptive(1).buttonfun = {'LeftArrow','RightArrow'}; 
    % adaptive staircase: corresponding signal values that would signify a
    % correct answer
    h.Settings.adaptive(1).signalval = [1 2];
    h.Settings.adaptive(1).response_type = 'samediff'; % discrimination is same/difference
    % reversals
    h.Settings.adaptive(1).reversals = [4;8;12];
    % stepsize
    h.Settings.adaptive(1).stepsize = [2;sqrt(2);sqrt(sqrt(2))];
    % steptype 0 = multiple/divide by stepsize; steptype 1 = add/subtract
    h.Settings.adaptive(1).steptype = 0;
    % stepdir -1 = level decreases intensity; stepdir 1 = level increases intensity
    h.Settings.adaptive(1).stepdir = -1;
    % starting level of adaptive staircase
    h.Settings.adaptive(1).startinglevel = 0.005; % should be a DIFFERENCE value. Keep small as it will increase naturally over time.
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
    % maximum amount of the difference value (should be a small fraction, e.g. 1/5th, of the stimulus intensity)
    h.Settings.adaptive(1).levelmax = 0.05; % should be a DIFFERENCE value.
    h.Settings.adaptive(1).levelmin = 0;
    h.Settings.adaptive(1).maxtrial = inf;
    %h.Settings.adaptive(1).expected_change = 5; % smaller value increases the precision of the prior for ZEST and reduces step size of changes in estimates
    h.Settings.adaptive(1).ignoretrials = 0;
    h.Settings.adaptive(1).eta_divide = 2;
    h.Settings.adaptive(1).slope_stimratio = 1; % number to divide stimulus level by, to calculate slope of ZEST. Should be smaller for thresholding, larger for adjustments during an expt
    

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
    h.Settings.wait=[2 0]; % one value per nstim 
    
    %% first stimulus: audio
    h.Settings.stim(1).patternmethod = 'pitch';% Pattern type method: intensity, pitch. Not supported: channel, duration
    % play with these:
    h.Settings.stim(1).dur = repmat(0,1,20); % duration of stimulus in seconds; modified by oddball settings. Value set to zero
    h.Settings.stim(1).stimrandind = [];% index of stimdur to randomise/adapt. 
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