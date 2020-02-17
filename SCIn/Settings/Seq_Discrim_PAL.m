function h = Seq_Discrim_PAL(h,opt)

switch opt
    
    case 'setoptions'
        
    % settings options
    h.SettingsOptions = {'US','CS','GSseries'};
    
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
    
    %% first stimulus: audio
    h.Settings.stim(1).patternmethod = 'pitch';% Pattern type method: intensity, pitch. Not supported: channel, duration
    % play with these:
    h.Settings.stim(1).dur = repmat(1/40,1,40); % duration of stimulus in seconds; modified by oddball settings. Value set to zero
    h.Settings.stim(1).stimrandind = [];% index of stimdur to randomise/adapt. 
    h.Settings.stim(1).patternvalue = repmat([300 0],1,40/2); % one per stimdur in each cell; one cell per oddball value
    h.Settings.stim(1).durtype = '';%'oddballvalue','sequence_rand'; 
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
    
    % within-trial frequency (s)
    h.Settings.wait=[]; % one value per nstim 
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, index, pitch, channel
    h.Settings.PL.oddballmethod = 'duration'; % can use same type for pattern only if oddball intensity is adaptive
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
    
    %% first stimulus: audio
    % inital settings/calculation to create gap series
    gapsum  = 0.5; % duration in s that gap series must sum to
    ngaps = 20; % number of gaps within series
    nx = 0; for i = 1:ngaps; nx = nx+i; end % triangular number - multipler of x to get train_dur
    x=gapsum/nx;
    gaps = x:x:(x*nincr);
    % main settings
    h.Settings.stim(1).patternmethod = 'pitch';% Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.stim(1).gap_ind=2:2:ngaps*2; % specifies the gaps so they can be randomised later
    h.Settings.stim(1).dur = [0.025*ones(1,ngaps) gaps]; % duration of stimulus in seconds. 0.0225 pips with gaps ranging from 0.01 to 0.1. Should sum to 1s for a 1s train.
    h.Settings.stim(1).dur([h.Settings.stim(1).gap_ind-1 h.Settings.stim(1).gap_ind]) = h.Settings.stim(1).dur; % adds stimuli to odd numbers and gaps to even numbers
    h.Settings.stim(1).stimrandind = h.Settings.stim(1).gap_ind(randperm(length(h.Settings.stim(1).gap_ind)));% randomises the gaps
    h.Settings.stim(1).patternvalue = repmat([300 0],1,ngaps); % one per stimdur in each cell; one cell per oddball value
    h.Settings.stim(1).durtype = '';%'oddballvalue','sequence_rand'; 
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
    
    % within-trial frequency (s)
    h.Settings.wait=[]; % one value per nstim 
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, index, pitch, channel
    h.Settings.PL.oddballmethod = 'duration'; % can use same type for pattern only if oddball intensity is adaptive
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
    
    %% first stimulus: audio
    
    % inital settings/calculation to create gap series
    gapsum  = 0.5; % duration in s that gap series must sum to
    ngaps = 20; % number of gaps within series
    nx = 0; for i = 1:ngaps; nx = nx+i; end % triangular number - multipler of x to get train_dur
    x=gapsum/nx;
    gaps = x:x:(x*ngaps);
    % main settings
    h.Settings.stim(1).patternmethod = 'pitch';% Pattern type method: intensity, pitch. Not supported: channel, duration
    h.Settings.stim(1).gap_ind=2:2:ngaps*2; % specifies the gaps so they can be randomised later
    % create basic form of the stimulus train with randomised gaps
    dur = [0.025*ones(1,ngaps) gaps]; % duration of stimulus in seconds. 0.0225 pips with gaps ranging from 0.01 to 0.1. Should sum to 1s for a 1s train.
    dur([h.Settings.stim(1).gap_ind-1 h.Settings.stim(1).gap_ind]) = dur; % adds stimuli to odd numbers and gaps to even numbers
    % create 7 versions in total
    h.Settings.stim(1).dur = {dur,dur,dur,dur,dur,dur,dur}; 
    sq = [1 3 5 6 4 2]; % create sequence sq. Values from 1:6 since we are creating 6 new stimulus trains. Numbers are indices of original gap order (ascending)
    sq=[sq,sq,sq,7,7]; % ensure sq adds up to the total number of gaps (e.g. 20) using 7 padding. 
    
    method = 2; 
    if method==1 % stimrandind randomisation
        stimrandind = h.Settings.stim(1).gap_ind(randperm(length(h.Settings.stim(1).gap_ind)));% randomises the gaps
    elseif method==2 % stimrandind randomisation using sq to evenly distribute
        stimrandind=nan(1,length(sq));
        for nr = 1:7
            old_ind = find(sq==nr);
            new_ind = old_ind(randperm(length(old_ind)));
            stimrandind(old_ind) = h.Settings.stim(1).gap_ind(new_ind);% randomises the gaps
        end
    end
    
    h.Settings.stim(1).gaprandind = sq(stimrandind/2); % random this index in same order as actual gaps. divide by 2 since we are converting from gap indices in whole train (even numbers) to indices of gaps only
    swap_order = [2 3 1]; % fixed order to swap the indices for each new train.
    h.Settings.stim(1).stimrandind{1} = stimrandind; % first train's gap order is stimrandind
    for nr = 2:7 % subsequent trains are a modification of the previous one
        old_ind = find(h.Settings.stim(1).gaprandind==(nr-1));
        new_ind = old_ind(swap_order);
        h.Settings.stim(1).stimrandind{nr} = h.Settings.stim(1).stimrandind{nr-1};
        h.Settings.stim(1).stimrandind{nr}(new_ind) = h.Settings.stim(1).stimrandind{nr-1}(old_ind);
    end
    % plot gap durations
    figure
    for nr = 1:7 
        subplot(7,1,nr)
        plotdur = dur(h.Settings.stim(1).stimrandind{nr});
        b=bar(plotdur);
        b.FaceColor = 'flat';
        if nr>1
            old_ind = find(h.Settings.stim(1).gaprandind==(nr-1));
            new_ind = old_ind(swap_order);
            for bi = new_ind
                b.CData(bi,:) = [1 0 0];
            end
        end
    end
    % the rest
    h.Settings.stim(1).patternvalue = repmat([300 0],1,ngaps); % one per stimdur in each cell; one cell per oddball value
    h.Settings.stim(1).durtype = '';%'oddballvalue','sequence_rand'; 
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
    
    % within-trial frequency (s)
    h.Settings.wait=[]; % one value per nstim 
    
    %% CHANGING STIMULUS INTENSITY EVERY X PULSES
    % REFER TO "TIMER STOP": https://labjack.com/support/ud/df-lj-app-guide/10.5
    
    %% Condition-dependent stimulus parameters
    % Condition method: intensity, pitch, channel
    h.Settings.conditionmethod = {};
    h.Settings.conditionvalue = [];% Rows: methods. Columns: each stimtype
    % Oddball method: intensity, index, pitch, channel
    h.Settings.PL.oddballmethod = 'duration'; % can use same type for pattern only if oddball intensity is adaptive
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
