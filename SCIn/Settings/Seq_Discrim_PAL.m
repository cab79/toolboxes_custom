function h = Seq_Discrim_PAL(h,opt)

% To use a T7, first install the software from here: https://labjack.com/sites/default/files/software/LabJack-2019-05-20.exe
% https://labjack.com/support/datasheets/t-series/communication/stream-mode/stream-out
% https://labjack.com/support/datasheets/t-series/communication/stream-mode/stream-out/stream-out-description

switch opt
    
    case 'setoptions'
        
    % settings options
    h.SettingsOptions = {'US','GS','GSseries','Ascend','Adaptive_pain','Adaptive_gaprand'};
    
    case 'US'
    % runs isochronic US stimuli at 20Hz for 1s each, once per trial.

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
    h.Settings.trialdur = inf; % if 0, consecutive stimuli will occur with no gap. 'Inf' requires participant to respond to move on to next trial.
    % Tactile: number of pulses per trial
    h.Settings.nstim_trial = 1; %
    % Which stims are targets for behavioural responses?
    h.Settings.target_stims = [1];
    
    %% first stimulus
    nstim=20; % minimum frequency of duration changes. Must divide by 2 to produce integer.
    train_dur = 5;
    % don't change these
    min_gap = 0.001; % minimum gap in ms. standard DS8R cannot be triggered more frequently than every 1ms. Use 0.001 for DS8R.
    min_trigger = max(0.0001,1/(16384/(2*train_dur))); % DS8R can detect down to 0.01ms but 0.1ms is sufficient for most purposes. Also need to take the sampling frequency into account.
    % inital settings/calculation to create gap series
    stimsum = nstim * min_trigger; % duration all the stimuli add up to in a train
    gapsum  = train_dur - stimsum; % duration in s that gap series must sum to
    ngaps = nstim-1; % number of gaps within series
    gaps = repmat(gapsum/ngaps,1,ngaps); % equally space gaps for US
    % check gaps are ok
    if min(gaps)<min_gap
        error('check gap minimum')
    end
    % create gaps and randomise
    gap_ind=2:2:ngaps*2; % specifies the gaps so they can be randomised later
    stim_ind=1:2:nstim*2; % specifies the gaps so they can be randomised later
    dur(stim_ind) = (stimsum/nstim)*ones(1,nstim); % duration of stimulus in seconds. 0.0225 pips with gaps ranging from 0.01 to 0.1. Should sum to 1s for a 1s train.
    dur(gap_ind) = gaps; 
    h.Settings.stim(1).dur = dur;
    % plot gap durations
    figure
    plotdur = h.Settings.stim(1).dur(gap_ind);
    bar(plotdur);
    title('gap durations')
    ylabel('seconds')
    % other settings
    h.Settings.stim(1).patternvalue = [repmat([5 0],1,ngaps) 5]; % one per stimdur in each cell; one cell per oddball value
    h.Settings.stim(1).durtype = '';%'oddballvalue','sequence_rand'; 
    h.Settings.stim(1).inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    if h.Settings.labjack
        % use labjack T7
        h.Settings.stim(1).patternvalue = [repmat([5 0],1,ngaps) 5]; % one per stimdur in each cell; one cell per oddball value
        h.Settings.stim(1).patternmethod = 'duration';% Pattern type method: intensity, pitch, duration. 
        h.Settings.stim(1).control='T7'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt, LJTick-DAQ, T7
        h.Settings.stim(1).inten_type = 'abs'; % either 'dB' or 'abs'
        h.Settings.stim(1).chan = 'DAC0'; h.Settings.stim(1).chanforLJ = 1;
        h.Settings.stim(1).nrchannels = 1; % total number of channels, e.g. on sound card
        h.Settings.stim(1).npulses_train = 0; % Tactile: number of pulses in a train; set to zero to be determined by stimdur
        h.Settings.stim(1).p_freq=[]; % Tactile: within-train frequency (Hz). only needed if not specifying using stimdur
        h.Settings.stim(1).labjack_timer=1; % Use timer to control frequency of labjack outputs? Otherwise uses software timing. must use if there are a large number of pulses (e.g. >4)
        h.Settings.stim(1).wavetype = 'digital'; % sin, square, step or digital.
        h.Settings.stim(1).f0 = 16384/(2*sum(h.Settings.stim(1).dur)); % frequency of output (1 / min duration). 16384 bytes is the max buffer size, and each output requires 2 bytes.
        h.Settings.stim(1).maxinten = 200; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    else
        % use audio
        h.Settings.stim(1).patternvalue = [repmat([300 0],1,ngaps) 300]; % one per stimdur in each cell; one cell per oddball value
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
    
    % within-trial frequency (s)
    h.Settings.wait=[]; % one value per nstim 
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, index, pitch, channel
    h.Settings.PL.oddballmethod = ''; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.PL.oddballvalue = {[1]}; % values to go into h.Seq.signal. One per oddprob row, or leave blank if determined from GUI
    h.Settings.PL.nonoddballstimvalue = {0};
    h.Settings.PL.oddballtype = 'classical'; % options: 'roving', 'classical'

    %% SEQUENCE
    % Change probablity (CP): each condition is in rows
    h.Settings.PL.oddprob = [
       1
        ];
    h.Settings.ntrials = [10];
    
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
    h.Settings.response_nexttrialmin = 1;
    % when does next trial start after button press? Empty if programmed
    % ISI
    h.Settings.response_nexttrialwait = 2;
    
    
    case 'GS'
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
    h.Settings.trialdur = inf; % if 0, consecutive stimuli will occur with no gap. 'Inf' requires participant to respond to move on to next trial.
    % Tactile: number of pulses per trial
    h.Settings.nstim_trial = 1; %
    % Which stims are targets for behavioural responses?
    h.Settings.target_stims = [1];
    
    %% first stimulus
    nstim=20; % minimum frequency of duration changes. Must divide by 2 to produce integer.
    train_dur = 5;
    % don't change these
    min_gap = 0.001; % minimum gap in ms. standard DS8R cannot be triggered more frequently than every 1ms. Use 0.001 for DS8R.
    min_trigger = max(0.0001,1/(16384/(2*train_dur))); % DS8R can detect down to 0.01ms but 0.1ms is sufficient for most purposes. Also need to take the sampling frequency into account.
    % inital settings/calculation to create gap series
    stimsum = nstim * min_trigger; % duration all the stimuli add up to in a train
    gapsum  = train_dur - stimsum; % duration in s that gap series must sum to
    ngaps = nstim-1; % number of gaps within series
    
    % this bit differs from the US code
    nx = 0; for i = 1:ngaps; nx = nx+i; end % triangular number - multipler of x to get train_dur
    x=gapsum/nx;
    gaps = x:x:(x*ngaps);
    
    % check gaps are ok
    if min(gaps)<min_gap
        error('check gap minimum')
    end
    % create gaps and randomise
    gap_ind=2:2:ngaps*2; % specifies the gaps so they can be randomised later
    stim_ind=1:2:nstim*2; % specifies the gaps so they can be randomised later
    dur(stim_ind) = (stimsum/nstim)*ones(1,nstim); % duration of stimulus in seconds. 0.0225 pips with gaps ranging from 0.01 to 0.1. Should sum to 1s for a 1s train.
    dur(gap_ind) = gaps; 
    h.Settings.stim(1).dur = dur;
    
    % this bit differs from the US code
    stimrandind = gap_ind(randperm(length(gap_ind))); % randomises the gaps
    h.Settings.stim(1).dur(sort(stimrandind)) = dur(stimrandind);
    
    % plot gap durations
    figure
    plotdur = h.Settings.stim(1).dur(sort(stimrandind));
    bar(plotdur);
    title('gap durations')
    ylabel('seconds')
    % other settings
    h.Settings.stim(1).durtype = '';%'oddballvalue','sequence_rand'; 
    h.Settings.stim(1).inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    if h.Settings.labjack
        % use labjack T7
        h.Settings.stim(1).patternvalue = [repmat([5 0],1,ngaps) 5]; % one per stimdur in each cell; one cell per oddball value
        h.Settings.stim(1).patternmethod = 'duration';% Pattern type method: intensity, pitch, duration. 
        h.Settings.stim(1).control='T7'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt, LJTick-DAQ, T7
        h.Settings.stim(1).inten_type = 'abs'; % either 'dB' or 'abs'
        h.Settings.stim(1).chan = 'DAC0'; h.Settings.stim(1).chanforLJ = 1;
        h.Settings.stim(1).nrchannels = 1; % total number of channels, e.g. on sound card
        h.Settings.stim(1).npulses_train = 0; % Tactile: number of pulses in a train; set to zero to be determined by stimdur
        h.Settings.stim(1).p_freq=[]; % Tactile: within-train frequency (Hz). only needed if not specifying using stimdur
        h.Settings.stim(1).labjack_timer=1; % Use timer to control frequency of labjack outputs? Otherwise uses software timing. must use if there are a large number of pulses (e.g. >4)
        h.Settings.stim(1).wavetype = 'digital'; % sin, square, step or digital.
        h.Settings.stim(1).f0 = 16384/(2*sum(h.Settings.stim(1).dur)); % frequency of output (1 / min duration). 16384 bytes is the max buffer size, and each output requires 2 bytes.
        h.Settings.stim(1).maxinten = 200; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    else
        % use audio
        h.Settings.stim(1).patternvalue = [repmat([300 0],1,ngaps) 300]; % one per stimdur in each cell; one cell per oddball value
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
    
    % within-trial frequency (s)
    h.Settings.wait=[]; % one value per nstim 
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, index, pitch, channel
    h.Settings.PL.oddballmethod = ''; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.PL.oddballvalue = {[1]}; % values to go into h.Seq.signal. One per oddprob row, or leave blank if determined from GUI
    h.Settings.PL.nonoddballstimvalue = {0};
    h.Settings.PL.oddballtype = 'classical'; % options: 'roving', 'classical'

    %% SEQUENCE
    % Change probablity (CP): each condition is in rows
    h.Settings.PL.oddprob = [
       1
        ];
    h.Settings.ntrials = [10];
    
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
    h.Settings.response_nexttrialmin = 1;
    % when does next trial start after button press? Empty if programmed
    % ISI
    h.Settings.response_nexttrialwait = 2;
    
    
    case 'GSseries'
    % runs random duration CS trains at (mean) 20Hz for 1s each, once per trial.

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
    h.Settings.trialdur = inf; % if 0, consecutive stimuli will occur with no gap. 'Inf' requires participant to respond to move on to next trial.
    % Tactile: number of pulses per trial
    h.Settings.nstim_trial = 3; 
    % Which stims are targets for behavioural responses?
    h.Settings.target_stims = [1 2 3];
    
    %% first stimulus
    
    % initial settings
    nstim=20; % minimum frequency of duration changes. Must divide by 2 to produce integer.
    train_dur = 1;
    
    % don't change these
    min_gap = 0.001; % minimum gap in ms. standard DS8R cannot be triggered more frequently than every 1ms. Use 0.001 for DS8R.
    min_trigger = max(0.0001,1/(16384/(2*train_dur))); % DS8R can detect down to 0.01ms but 0.1ms is sufficient for most purposes. Also need to take the sampling frequency into account.
    % min_trigger = 0.05; % AUDIO VERSION
    
    % inital settings/calculation to create gap series
    stimsum = nstim * min_trigger; % duration all the stimuli add up to in a train
    gapsum  = train_dur - stimsum; % duration in s that gap series must sum to
    ngaps = nstim-1; % number of gaps within series
    nx = 0; for i = 1:ngaps; nx = nx+i; end % triangular number - multipler of x to get train_dur
    x=gapsum/nx;
    gaps = x:x:(x*ngaps);
    % check gaps are ok
    if min(gaps)<min_gap
        error('check gap minimum')
    end
    % create gaps and randomise
    gap_ind=2:2:ngaps*2; % specifies the gaps so they can be randomised later
    stim_ind=1:2:nstim*2; % specifies the gaps so they can be randomised later
    dur(stim_ind) = (stimsum/nstim)*ones(1,nstim); % duration of stimulus in seconds. 0.0225 pips with gaps ranging from 0.01 to 0.1. Should sum to 1s for a 1s train.
    dur(gap_ind) = gaps; 
    
    % create 7 versions in total
    h.Settings.stim(1).dur = {dur,dur,dur,dur,dur,dur,dur};
    % define sq
    sq = [1 3 5 6 4 2]; % create sequence sq. Values from 1:6 since we are creating 6 new stimulus trains. Numbers are indices of original gap order (ascending)
    sq=[sq,sq,sq,7]; % ensure sq adds up to the total number of gaps (e.g. 19) using 7 padding. 
    
    % create initial randomised train (e.g. for GS1)
    method = 2; 
    if method==1 % stimrandind randomisation
        stimrandind = gap_ind(randperm(length(gap_ind)));% randomises the gaps
    elseif method==2 % stimrandind randomisation using sq to evenly distribute
        stimrandind=nan(1,length(sq));
        for nr = 1:7
            old_ind = find(sq==nr);
            new_ind = old_ind(randperm(length(old_ind)));
            stimrandind(old_ind) = gap_ind(new_ind);% randomises the gaps
        end
    end
    h.Settings.stim(1).gaprandind = sq(stimrandind/2); % random this index in same order as actual gaps. divide by 2 since we are converting from gap indices in whole train (even numbers) to indices of gaps only
    
    % sequential swapping of gaps to create 7 trains:
    swap_order = [2 3 1]; % fixed order to swap the indices for each new train.
    h.Settings.stim(1).stimrandind{1} = stimrandind; % first train's gap order is stimrandind
    h.Settings.stim(1).dur{1}(sort(stimrandind)) = dur(h.Settings.stim(1).stimrandind{1});
    for nr = 2:7 % subsequent trains are a modification of the previous one
        old_ind = find(h.Settings.stim(1).gaprandind==(nr-1));
        new_ind = old_ind(swap_order);
        h.Settings.stim(1).stimrandind{nr} = h.Settings.stim(1).stimrandind{nr-1};
        h.Settings.stim(1).stimrandind{nr}(new_ind) = h.Settings.stim(1).stimrandind{nr-1}(old_ind);
        h.Settings.stim(1).dur{nr}(sort(stimrandind)) = dur(h.Settings.stim(1).stimrandind{nr});
    end
    
    % plot gap durations
    figure
    for nr = 1:7 
        subplot(7,1,nr)
        plotdur = h.Settings.stim(1).dur{nr}(sort(stimrandind));
        b=bar(plotdur);
        ylabel('seconds')
        b.FaceColor = 'flat';
        if nr>1
            old_ind = find(h.Settings.stim(1).gaprandind==(nr-1));
            new_ind = old_ind(swap_order);
            for bi = new_ind
                b.CData(bi,:) = [1 0 0];
            end
        end
    end
    % other settings
    h.Settings.stim(1).durtype = 'oddballvalue';% 'oddballvalue' will select which sequence is presented on each trial according to the values in h.Settings.PL.oddballvalue
    h.Settings.stim(1).inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    if h.Settings.labjack
        % use labjack T7
        h.Settings.stim(1).patternvalue = [repmat([5 0],1,ngaps) 5]; % one per stimdur in each cell; one cell per oddball value
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
        h.Settings.stim(1).patternvalue = [repmat([300 0],1,ngaps) 300]; % one per stimdur in each cell; one cell per oddball value
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
    
    % duplicate so that there are h.Settings.nstim_trial stimuli per trial
    h.Settings.stim(2) = h.Settings.stim(1);
    h.Settings.stim(3) = h.Settings.stim(1);
    
    % within-trial stimulus onset asychrony (SOA) (s)
    h.Settings.wait=[6 6 6]; % one value per nstim
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, index, pitch, channel
    h.Settings.PL.oddballmethod = 'pattern'; % this will select different patterns according to h.Settings.PL.oddballvalue. Can use same type for pattern only if oddball intensity is adaptive
    h.Settings.PL.oddballvalue = {[1 2 2],[2 1 2],[2 2 1],[6 7 7],[7 6 7],[7 7 6]}; % values to go into h.Seq.signal. One per oddprob row, or leave blank if determined from GUI
    h.Settings.PL.nonoddballstimvalue = {0 0 0 0 0 0};
    h.Settings.PL.oddballtype = 'classical'; % options: 'roving', 'classical'

    %% SEQUENCE
    % Change probablity (CP): each condition is in rows
    h.Settings.PL.oddprob = [
        1 % for triangular sequence, each row here just has to add up to one. Simplest just to put 1 for each row.
        1
        1
        1
        1
        1
        ];
    h.Settings.ntrials = [9 9 9 9 9 9]; % if 3 stims per trial should multiples of 3 to allow every combination as programmed in oddballvalue. 
    
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
    h.Settings.response_nexttrialmin = 1;
    % when does next trial start after button press? Empty if programmed
    % ISI
    h.Settings.response_nexttrialwait = 2;
    
    
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
    h.Settings.nstim_trial = 1; % set to zero to be determined by stimdur
    h.Settings.wait=0; % within-trial frequency (Hz); one value per nstim 
    
    %% stimulus
    nstim=20; % minimum frequency of duration changes. Must divide by 2 to produce integer.
    train_dur = 1;
    h.Settings.trialdur = 2; % should be longer than train_dur
    % don't change these
    min_gap = 0.001; % minimum gap in ms. standard DS8R cannot be triggered more frequently than every 1ms. Use 0.001 for DS8R.
    min_trigger = max(0.0001,1/(16384/(2*train_dur))); % DS8R can detect down to 0.01ms but 0.1ms is sufficient for most purposes. Also need to take the sampling frequency into account.
    % inital settings/calculation to create gap series
    stimsum = nstim * min_trigger; % duration all the stimuli add up to in a train
    gapsum  = train_dur - stimsum; % duration in s that gap series must sum to
    ngaps = nstim-1; % number of gaps within series
    gaps = repmat(gapsum/ngaps,1,ngaps); % equally space gaps for US
    % check gaps are ok
    if min(gaps)<min_gap
        error('check gap minimum')
    end
    % create gaps and randomise
    gap_ind=2:2:ngaps*2; % specifies the gaps so they can be randomised later
    stim_ind=1:2:nstim*2; % specifies the gaps so they can be randomised later
    dur(stim_ind) = (stimsum/nstim)*ones(1,nstim); % duration of stimulus in seconds. 0.0225 pips with gaps ranging from 0.01 to 0.1. Should sum to 1s for a 1s train.
    dur(gap_ind) = gaps; 
    h.Settings.stim(1).dur = dur;
    % plot gap durations
    figure
    plotdur = h.Settings.stim(1).dur(gap_ind);
    bar(plotdur);
    title('gap durations')
    ylabel('seconds')
    % other settings
    h.Settings.stim(1).patternvalue = [repmat([5 0],1,ngaps) 5]; % one per stimdur in each cell; one cell per oddball value
    h.Settings.stim(1).durtype = '';%'oddballvalue','sequence_rand'; 
    h.Settings.stim(1).inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    if h.Settings.labjack
        % use labjack T7
        h.Settings.stim(1).patternvalue = [repmat([5 0],1,ngaps) 5]; % one per stimdur in each cell; one cell per oddball value
        h.Settings.stim(1).patternmethod = 'duration';% Pattern type method: intensity, pitch, duration. 
        h.Settings.stim(1).control='T7'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt, LJTick-DAQ, T7
        h.Settings.stim(1).inten_type = 'abs'; % either 'dB' or 'abs'
        h.Settings.stim(1).chan = 'DAC0'; h.Settings.stim(1).chanforLJ = 1;
        h.Settings.stim(1).nrchannels = 1; % total number of channels, e.g. on sound card
        h.Settings.stim(1).npulses_train = 0; % Tactile: number of pulses in a train; set to zero to be determined by stimdur
        h.Settings.stim(1).p_freq=[]; % Tactile: within-train frequency (Hz). only needed if not specifying using stimdur
        h.Settings.stim(1).labjack_timer=1; % Use timer to control frequency of labjack outputs? Otherwise uses software timing. must use if there are a large number of pulses (e.g. >4)
        h.Settings.stim(1).wavetype = 'digital'; % sin, square, step or digital.
        h.Settings.stim(1).f0 = 16384/(2*sum(h.Settings.stim(1).dur)); % frequency of output (1 / min duration). 16384 bytes is the max buffer size, and each output requires 2 bytes.
        h.Settings.stim(1).maxinten = 200; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    else
        % use audio
        h.Settings.stim(1).patternvalue = [repmat([300 0],1,ngaps) 300]; % one per stimdur in each cell; one cell per oddball value
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
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, pitch, channel
    h.Settings.PL.oddballmethod = ''; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.PL.oddballvalue = [];
    h.Settings.PL.oddballtype = ''; % options: 'roving', 'classical'

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
    
    case 'Adaptive_pain'

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
    
    %% stimulus
    nstim=20; % minimum frequency of duration changes. Must divide by 2 to produce integer.
    train_dur = 1;
    % don't change these
    min_gap = 0.001; % minimum gap in ms. standard DS8R cannot be triggered more frequently than every 1ms. Use 0.001 for DS8R.
    min_trigger = max(0.0001,1/(16384/(2*train_dur))); % DS8R can detect down to 0.01ms but 0.1ms is sufficient for most purposes. Also need to take the sampling frequency into account.
    % inital settings/calculation to create gap series
    stimsum = nstim * min_trigger; % duration all the stimuli add up to in a train
    gapsum  = train_dur - stimsum; % duration in s that gap series must sum to
    ngaps = nstim-1; % number of gaps within series
    gaps = repmat(gapsum/ngaps,1,ngaps); % equally space gaps for US
    % check gaps are ok
    if min(gaps)<min_gap
        error('check gap minimum')
    end
    % create gaps and randomise
    gap_ind=2:2:ngaps*2; % specifies the gaps so they can be randomised later
    stim_ind=1:2:nstim*2; % specifies the gaps so they can be randomised later
    dur(stim_ind) = (stimsum/nstim)*ones(1,nstim); % duration of stimulus in seconds. 0.0225 pips with gaps ranging from 0.01 to 0.1. Should sum to 1s for a 1s train.
    dur(gap_ind) = gaps; 
    h.Settings.stim(1).dur = dur;
    % plot gap durations
    figure
    plotdur = h.Settings.stim(1).dur(gap_ind);
    bar(plotdur);
    title('gap durations')
    ylabel('seconds')
    % other settings
    h.Settings.stim(1).patternvalue = [repmat([5 0],1,ngaps) 5]; % one per stimdur in each cell; one cell per oddball value
    h.Settings.stim(1).durtype = '';%'oddballvalue','sequence_rand'; 
    h.Settings.stim(1).inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    if h.Settings.labjack
        % use labjack T7
        h.Settings.stim(1).patternvalue = [repmat([5 0],1,ngaps) 5]; % one per stimdur in each cell; one cell per oddball value
        h.Settings.stim(1).patternmethod = 'duration';% Pattern type method: intensity, pitch, duration. 
        h.Settings.stim(1).control='T7'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt, LJTick-DAQ, T7
        h.Settings.stim(1).inten_type = 'abs'; % either 'dB' or 'abs'
        h.Settings.stim(1).chan = 'DAC0'; h.Settings.stim(1).chanforLJ = 1;
        h.Settings.stim(1).nrchannels = 1; % total number of channels, e.g. on sound card
        h.Settings.stim(1).npulses_train = 0; % Tactile: number of pulses in a train; set to zero to be determined by stimdur
        h.Settings.stim(1).p_freq=[]; % Tactile: within-train frequency (Hz). only needed if not specifying using stimdur
        h.Settings.stim(1).labjack_timer=1; % Use timer to control frequency of labjack outputs? Otherwise uses software timing. must use if there are a large number of pulses (e.g. >4)
        h.Settings.stim(1).wavetype = 'digital'; % sin, square, step or digital.
        h.Settings.stim(1).f0 = 16384/(2*sum(h.Settings.stim(1).dur)); % frequency of output (1 / min duration). 16384 bytes is the max buffer size, and each output requires 2 bytes.
        h.Settings.stim(1).maxinten = 200; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    else
        % use audio
        h.Settings.stim(1).patternvalue = [repmat([300 0],1,ngaps) 300]; % one per stimdur in each cell; one cell per oddball value
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
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, index, pitch, channel
    h.Settings.PL.oddballmethod = 'index'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.PL.oddballvalue = {[1 2]}; % values to go into h.Seq.signal. One per oddprob row, or leave blank if determined from GUI
    h.Settings.PL.oddballtype = 'classical'; % options: 'roving', 'classical'
    
    %% SEQUENCE
    % Change probablity (CP): each condition is in rows
    h.Settings.PL.oddprob = [
        % standard (left) vs oddball (right)
        0.5 0.5
        ];
    % index of row of oddprob that are standards and oddballs. Can be
    % overridden by h.Settings.oddballvalue if using index
    h.Settings.PL.standardind = 1; 
    h.Settings.PL.oddind = 2; 
    % keep oddball trials apart by at least sep_odd standards
    h.Settings.PL.sep_odd = [0]; % for each CP condition
    % for sep_odd, which indices of h.Settings.oddballvalue to consider
    % each time? (each list will be considered separately)
    h.Settings.PL.sep_odd_ind = {[1 2]};
    h.Settings.PL.sep_odd_tol = [1]; % set these to be as high as possible (max 1)
    % for each set, ensure a number of leading standards 
    h.Settings.PL.std_lead = [0]; % for each CP condition
    % number of sets to randomise together
    h.Settings.PL.n_set = []; % Leave blank to calculate automatically; or one nunmber per CP condition
    % min number of oddballs within each CP condition
    h.Settings.PL.n_odd = 10*[12]; % overrides h.Settings.totdur
    % min number of oddballs per randomised set, per CP
    h.Settings.PL.n_odd_set = 10*[12]; % overrides h.Settings.totdur
    % randomise sets?
    h.Settings.PL.rand_set = [0]; 
    
    
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
    h.Settings.adaptive_general.adapttypes = [1];
    % alternate or randomise runs over types? Alt must have equal number of
    % runs for each adapttype. Cond = one type per CP block
    h.Settings.adaptive_general.seqtype = 'rand'; % 'alt', 'rand', 'cond' 
    h.Settings.adaptive_general.seqtypecond = [1]; %if 'cond', associate each CP with an adaptive type
    h.Settings.adaptive_general.seqrandblocksize = 12; % should divide the number of trials in a set
    h.Settings.adaptive_general.selectcond.cp = [1]; % which CP condition to run adaptive on?
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
    
    
    case 'Adaptive_gaprand'

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
    h.Settings.nstim_trial = 2; % set to zero to be determined by stimdur
    % Which stims that are targets for behavioural responses?
    h.Settings.target_stims = [1 2];
    
    %% first stimulus
    nstim=20; % minimum frequency of duration changes. Must divide by 2 to produce integer.
    train_dur = 1;
    % don't change these
    min_gap = 0.001; % minimum gap in ms. standard DS8R cannot be triggered more frequently than every 1ms. Use 0.001 for DS8R.
    min_trigger = max(0.0001,1/(16384/(2*train_dur))); % DS8R can detect down to 0.01ms but 0.1ms is sufficient for most purposes. Also need to take the sampling frequency into account.
    % inital settings/calculation to create gap series
    stimsum = nstim * min_trigger; % duration all the stimuli add up to in a train
    gapsum  = train_dur - stimsum; % duration in s that gap series must sum to
    ngaps = nstim-1; % number of gaps within series
    % this bit differs from the US code
    nx = 0; for i = 1:ngaps; nx = nx+i; end % triangular number - multipler of x to get train_dur
    x=gapsum/nx;
    gaps = x:x:(x*ngaps);
    % check gaps are ok
    if min(gaps)<min_gap
        error('check gap minimum')
    end
    % create gaps and randomise
    gap_ind=2:2:ngaps*2; % specifies the gaps so they can be randomised later
    stim_ind=1:2:nstim*2; % specifies the gaps so they can be randomised later
    dur(stim_ind) = (stimsum/nstim)*ones(1,nstim); % duration of stimulus in seconds. 0.0225 pips with gaps ranging from 0.01 to 0.1. Should sum to 1s for a 1s train.
    dur(gap_ind) = gaps; 
    h.Settings.stim(1).dur = dur;
    % this bit differs from the US code
    stimrandind = gap_ind(randperm(length(gap_ind))); % randomises the gaps
    h.Settings.stim(1).dur(sort(stimrandind)) = dur(stimrandind);
    % plot gap durations
    figure
    plotdur = h.Settings.stim(1).dur(sort(stimrandind));
    bar(plotdur);
    title('gap durations')
    ylabel('seconds')
    
    % these settings are unique to adaptive procedure for gap
    % randomisations
    h.Settings.stim(1).ngap_rand_oddball=4; % initial number of gaps to randomise (which update using adaptive procedure during experiment)
    h.Settings.stim(1).stimrandind_durind = [1 2]; % select which values of signal to apply stimrandind to
    h.Settings.stim(1).stimrandind_oddball = [13 10 7 4 1 8 11 14]*2; % indices to randomise as oddballs
    
    % other settings
    h.Settings.stim(1).durtype = 'oddballvalue';%'oddballvalue','sequence_rand'; 
    h.Settings.stim(1).inten = 0; % value between 2 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).inten_diff_max = []; % value between 0 and 1000mA for Digitimer DS8R
    h.Settings.stim(1).maxinten = 0; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    if h.Settings.labjack
        % use labjack T7
        h.Settings.stim(1).patternvalue = [repmat([5 0],1,ngaps) 5]; % one per stimdur in each cell; one cell per oddball value
        h.Settings.stim(1).patternmethod = 'duration';% Pattern type method: intensity, pitch, duration. 
        h.Settings.stim(1).control='T7'; % How to control stimulator? Options: PsychPortAudio, audioplayer, labjack, spt, LJTick-DAQ, T7
        h.Settings.stim(1).inten_type = 'abs'; % either 'dB' or 'abs'
        h.Settings.stim(1).chan = 'DAC0'; h.Settings.stim(1).chanforLJ = 1;
        h.Settings.stim(1).nrchannels = 1; % total number of channels, e.g. on sound card
        h.Settings.stim(1).npulses_train = 0; % Tactile: number of pulses in a train; set to zero to be determined by stimdur
        h.Settings.stim(1).p_freq=[]; % Tactile: within-train frequency (Hz). only needed if not specifying using stimdur
        h.Settings.stim(1).labjack_timer=1; % Use timer to control frequency of labjack outputs? Otherwise uses software timing. must use if there are a large number of pulses (e.g. >4)
        h.Settings.stim(1).wavetype = 'digital'; % sin, square, step or digital.
        h.Settings.stim(1).f0 = 16384/(2*sum(h.Settings.stim(1).dur)); % frequency of output (1 / min duration). 16384 bytes is the max buffer size, and each output requires 2 bytes.
        h.Settings.stim(1).maxinten = 200; % max output value for safety purposes. Value between 2 and 1000mA for Digitimer DS8R
    else
        % use audio
        h.Settings.stim(1).patternvalue = [repmat([300 0],1,ngaps) 300]; % one per stimdur in each cell; one cell per oddball value
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
    
    % duplicate
    for ns=[1 2]
        h.Settings.stim(ns) = h.Settings.stim(1);
    end
    
    
    % This bit was previously in Triangular sequence function but I removed
    % to put it here since it's experiment-specific
    % It simply duplicates the pattern so that one of them can have the
    % gaps randomised and the other one kept the same.
    for s = 1:length(h.Settings.stim)
        if ~iscell(h.Settings.stim(s).patternvalue)
            h.Settings.stim(s).patternvalue = {h.Settings.stim(s).patternvalue};
            h.Settings.stim(s).dur = {h.Settings.stim(s).dur};
        end
        patternvalue = h.Settings.stim(s).patternvalue;
        dur = h.Settings.stim(s).dur;
        for i = 1:length(patternvalue)
            h.Settings.stim(s).patternvalue{i*2-1} = patternvalue{i};
            h.Settings.stim(s).patternvalue{i*2} = patternvalue{i};
            h.Settings.stim(s).dur{i*2-1} = dur{i};
            h.Settings.stim(s).dur{i*2} = dur{i};
        end
    end
    
    % within-trial frequency (s)
    h.Settings.wait=[sum(h.Settings.stim(1).dur{1})+1 0]; % one value per nstim
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    
    % Oddball method
    h.Settings.PL.oddballmethod = 'duration_shuffle'; % can use same type for pattern only if oddball intensity is adaptive
    h.Settings.PL.oddballvalue = {[1 1],[1 2],[2 1],[2 2]}; % values to go into h.Seq.signal. One per oddprob row, or leave blank if determined from GUI
    h.Settings.PL.nonoddballstimvalue = {0 0 0 0};
    h.Settings.PL.oddballtype = 'classical'; % options: 'roving', 'classical'
    
    %% SEQUENCE
    h.Settings.PL.oddprob = [
        1/2 1/2
        1/2 1/2
        1/2 1/2
        1/2 1/2
        ];
%     % index of row of oddprob that are standards and oddballs. Can be
%     % overridden by h.Settings.oddballvalue if using index
%     h.Settings.PL.standardind = 1; 
%     h.Settings.PL.oddind = 2; 
%     % keep oddball trials apart by at least sep_odd standards
%     h.Settings.PL.sep_odd = [0]; % for each CP condition
%     % for sep_odd, which indices of h.Settings.oddballvalue to consider
%     % each time? (each list will be considered separately)
%     h.Settings.PL.sep_odd_ind = {[1 2]};
%     h.Settings.PL.sep_odd_tol = [1]; % set these to be as high as possible (max 1)
%     % for each set, ensure a number of leading standards 
%     h.Settings.PL.std_lead = [0]; % for each CP condition
%     % number of sets to randomise together
%     h.Settings.PL.n_set = []; % Leave blank to calculate automatically; or one nunmber per CP condition
%     % min number of oddballs within each CP condition
%     h.Settings.PL.n_odd = 10*[12]; % overrides h.Settings.totdur
%     % min number of oddballs per randomised set, per CP
%     h.Settings.PL.n_odd_set = 10*[12]; % overrides h.Settings.totdur
%     % randomise sets?
%     h.Settings.PL.rand_set = [0]; 
    h.Settings.ntrials = [10 10 10 10];
    h.Settings.randomise_conds_evenly = 1;
    
    
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
    h.Settings.adaptive_general.seqtypecond = [1]; %if 'cond', associate each CP with an adaptive type
    h.Settings.adaptive_general.seqrandblocksize = 12; % should divide the number of trials in a set
    h.Settings.adaptive_general.selectcond.cp = [1]; % which CP condition to run adaptive on?
    h.Settings.adaptive_general.stim = h.Settings.target_stims; % which stims to apply changes to?
    h.Settings.adaptive_general.stim_judge = h.Settings.target_stims; % which stims to participants judge?
    h.Settings.adaptive_general.terminate = 'block'; % terminate within each block only
    h.Settings.adaptive_general.reestimate = ''; % 'block' to re-estimate with wider prior each block
    
    %% ADAPTIVE
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
    h.Settings.adaptive(1).signalcorrect = []; % correct response
    h.Settings.adaptive(1).signaltarget = {[1 2],[1 2]}; % translate actual signal values to signalval indices
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
    h.Settings.adaptive(1).levelmax = length(gap_ind); % should be a DIFFERENCE value.
    h.Settings.adaptive(1).levelmin = 0;
    h.Settings.adaptive(1).maxtrial = inf;
    %h.Settings.adaptive(1).expected_change = 5; % smaller value increases the precision of the prior for ZEST and reduces step size of changes in estimates
    h.Settings.adaptive(1).ignoretrials = 0;
    h.Settings.adaptive(1).eta_divide = 2;
    h.Settings.adaptive(1).slope_stimratio = 1; % number to divide stimulus level by, to calculate slope of ZEST. Should be smaller for thresholding, larger for adjustments during an expt
    
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

if 1 % turn to 1 to use DS8R without EEG
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

if 0 % turn to 1 to use audio without EEG
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
