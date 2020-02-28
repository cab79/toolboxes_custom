function h = Seq_Discrim_PAL(h,opt)

% To use a T7, first install the software from here: https://labjack.com/sites/default/files/software/LabJack-2019-05-20.exe
% https://labjack.com/support/datasheets/t-series/communication/stream-mode/stream-out
% https://labjack.com/support/datasheets/t-series/communication/stream-mode/stream-out/stream-out-description

switch opt
    
    case 'setoptions'
        
    % settings options
    h.SettingsOptions = {'US','GS','GSseries'};
    
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
    min_gap = 0.001; % minimum gap in ms. standard DS8R cannot be triggered more frequently than every 1ms. 
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
    min_gap = 0.001; % minimum gap in ms. standard DS8R cannot be triggered more frequently than every 1ms. 
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
    nstim=20; % minimum frequency of duration changes. Must divide by 2 to produce integer.
    train_dur = 5;
    % don't change these:
    min_gap = 0.001; % minimum gap in ms. standard DS8R cannot be triggered more frequently than every 1ms. Use 0.001 for DS8R.
%     min_trigger = max(0.0001,1/(16384/(2*train_dur))); % DS8R can detect down to 0.01ms but 0.1ms is sufficient for most purposes. Also need to take the sampling frequency into account.
    min_trigger = 0.05; % AUDIO VERSION. 0.05 = 50ms.
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
    h.Settings.DAC_multiply = 0.01; % multiply DAC output by this (e.g. to get mA on DS8R)
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
