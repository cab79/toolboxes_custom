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
    h.Settings.trialdur = 2; % if 0, consecutive stimuli will occur with no gap. 'Inf' requires participant to respond to move on to next trial.
    % Tactile: number of pulses per trial
    h.Settings.nstim_trial = 1; 
    % Which stims are targets for behavioural responses?
    h.Settings.target_stims = [1];
    
    %% first stimulus
    
%     % OPTION 1: stimulus intervals. create 9 trial types in 50ms steps
%     nstim=4; % minimum frequency of duration changes. Must divide by 2 to produce integer. In this example: on/off/on/off = 4.
%     train_dur = 0.64; % 640ms
%     gap = 0.05;
%     pip = 0.030;
%     step = 0.050;
%     h.Settings.stim(1).dur = {
%         [pip gap pip train_dur-gap-2*pip] % e.g. 15,50,15,180 for on/off/on/off
%         [pip gap+1*step pip train_dur-(gap+1*step)-2*pip] 
%         [pip gap+2*step pip train_dur-(gap+2*step)-2*pip]
%         [pip gap+3*step pip train_dur-(gap+3*step)-2*pip]
%         [pip gap+4*step pip train_dur-(gap+4*step)-2*pip]
%         [pip gap+5*step pip train_dur-(gap+5*step)-2*pip]
%         [pip gap+6*step pip train_dur-(gap+6*step)-2*pip]
%         [pip gap+7*step pip train_dur-(gap+7*step)-2*pip]
%         [pip gap+8*step pip train_dur-(gap+8*step)-2*pip] % e.g. 15,210,15,0
%         };
%     ngaps=nstim/2;
%     h.Settings.stim(1).patternvalue = repmat([400 0],1,ngaps);
    
%     % OPTION 2: contant filled stimuli. create 9 trial types in 50ms steps
%     nstim=2; % minimum frequency of duration changes. Must divide by 2 to produce integer. In this example: on/off/on/off = 4.
%     train_dur = 0.64; % 640ms
%     gap = 0.05;
%     %pip = 0.030;
%     step = 0.050;
%     h.Settings.stim(1).dur = {
%         [gap train_dur-gap] % e.g. 15,50,15,180 for on/off/on/off
%         [gap+1*step train_dur-(gap+1*step)] 
%         [gap+2*step train_dur-(gap+2*step)]
%         [gap+3*step train_dur-(gap+3*step)]
%         [gap+4*step train_dur-(gap+4*step)]
%         [gap+5*step train_dur-(gap+5*step)]
%         [gap+6*step train_dur-(gap+6*step)]
%         [gap+7*step train_dur-(gap+7*step)]
%         [gap+8*step train_dur-(gap+8*step)] % e.g. 15,210,15,0
%         };
%     ngaps=1;
%     h.Settings.stim(1).patternvalue = repmat([400 0],1,ngaps);
    
%     % OPTION 3: contant filled stimuli. create 9 trial types in 50ms steps
%     nstim=2; % minimum frequency of duration changes. Must divide by 2 to produce integer. In this example: on/off/on/off = 4.
%     train_dur = 2; % 640ms
%     gap = 0.05;
%     %pip = 0.030;
%     step = 1.4;
%     h.Settings.stim(1).dur = {
%         [gap train_dur-gap] % e.g. 15,50,15,180 for on/off/on/off
%         [gap*step train_dur-(gap*step)] 
%         [gap*step^2 train_dur-(gap*step^2)]
%         [gap*step^3 train_dur-(gap*step^3)]
%         [gap*step^4 train_dur-(gap*step^4)]
%         [gap*step^5 train_dur-(gap*step^5)]
%         [gap*step^6 train_dur-(gap*step^6)]
%         [gap*step^7 train_dur-(gap*step^7)]
%         [gap*step^8 train_dur-(gap*step^8)] % e.g. 15,210,15,0
%         };
%     ngaps=1;
%     h.Settings.stim(1).patternvalue = repmat([400 0],1,ngaps);
    
%     % OPTION 4: stimulus intervals - multiplier. create 9 trial types in 50ms steps
%     nstim=2; % minimum frequency of duration changes. Must divide by 2 to produce integer. In this example: on/off/on/off = 4.
%     train_dur = 2; % 640ms
%     gap = 0.05;
%     pip = 0.030;
%     step = 1.4;
%     h.Settings.stim(1).dur = {
%         [pip gap pip train_dur-gap-2*pip] % e.g. 15,50,15,180 for on/off/on/off
%         [pip gap*step pip train_dur-(gap*step)-2*pip] 
%         [pip gap*step^2 pip train_dur-(gap*step^2)-2*pip]
%         [pip gap*step^3 pip train_dur-(gap*step^3)-2*pip]
%         [pip gap*step^4 pip train_dur-(gap*step^4)-2*pip]
%         [pip gap*step^5 pip train_dur-(gap*step^5)-2*pip]
%         [pip gap*step^6 pip train_dur-(gap*step^6)-2*pip]
%         [pip gap*step^7 pip train_dur-(gap*step^7)-2*pip]
%         [pip gap*step^8 pip train_dur-(gap*step^8)-2*pip] % e.g. 15,210,15,0
%         };
%     ngaps=2;
%     h.Settings.stim(1).patternvalue = repmat([400 0],1,ngaps);
    
%     % OPTION 5: frequencies
%     train_dur = 1; % seconds
%     freq1 = 10; % first frequency, Hz
%     step = 1.2; % multiplier of freq
%     h.Settings.stim(1).patternvalue = {
%         repmat([400 0],1,freq1)
%         repmat([400 0],1,round(freq1*step)) 
%         repmat([400 0],1,round(freq1*step^2)) 
%         repmat([400 0],1,round(freq1*step^3)) 
%         repmat([400 0],1,round(freq1*step^4)) 
%         repmat([400 0],1,round(freq1*step^5)) 
%         repmat([400 0],1,round(freq1*step^6)) 
%         repmat([400 0],1,round(freq1*step^7)) 
%         repmat([400 0],1,round(freq1*step^8)) 
%         %repmat([400 0],1,round(freq1*step^9)) 
%         };
%     h.Settings.stim(1).dur = {
%         repmat([1/(freq1*2)],1,length(h.Settings.stim(1).patternvalue{1}))
%         repmat([1/(freq1*step*2)],1,length(h.Settings.stim(1).patternvalue{2})) 
%         repmat([1/(freq1*step^2*2)],1,length(h.Settings.stim(1).patternvalue{3})) 
%         repmat([1/(freq1*step^3*2)],1,length(h.Settings.stim(1).patternvalue{4})) 
%         repmat([1/(freq1*step^4*2)],1,length(h.Settings.stim(1).patternvalue{5})) 
%         repmat([1/(freq1*step^5*2)],1,length(h.Settings.stim(1).patternvalue{6})) 
%         repmat([1/(freq1*step^6*2)],1,length(h.Settings.stim(1).patternvalue{7})) 
%         repmat([1/(freq1*step^7*2)],1,length(h.Settings.stim(1).patternvalue{8})) 
%         repmat([1/(freq1*step^8*2)],1,length(h.Settings.stim(1).patternvalue{9})) 
%         %repmat([1/(freq1*step^9*2)],1,length(h.Settings.stim(1).patternvalue{10})) 
%         };
%     
%     % OPTION 6: frequency AND duration
%     train_dur = 2; % seconds
%     freq1 = 10; % first frequency, Hz
%     step = 1.2; % multiplier of freq
%     h.Settings.stim(1).patternvalue = {
%         repmat([400 0],1,freq1)
%         repmat([400 0],1,freq1) 
%         repmat([400 0],1,freq1) 
%         repmat([400 0],1,freq1) 
%         repmat([400 0],1,freq1) 
%         repmat([400 0],1,freq1) 
%         repmat([400 0],1,freq1) 
%         repmat([400 0],1,freq1) 
%         repmat([400 0],1,freq1) 
%         %repmat([400 0],1,round(freq1*step^9)) 
%         };
%     h.Settings.stim(1).dur = {
%         repmat([1/(freq1*2)],1,length(h.Settings.stim(1).patternvalue{1}))
%         repmat([1/(freq1*step*2)],1,length(h.Settings.stim(1).patternvalue{2}))
%         repmat([1/(freq1*step^2*2)],1,length(h.Settings.stim(1).patternvalue{3}))
%         repmat([1/(freq1*step^3*2)],1,length(h.Settings.stim(1).patternvalue{4}))
%         repmat([1/(freq1*step^4*2)],1,length(h.Settings.stim(1).patternvalue{5}))
%         repmat([1/(freq1*step^5*2)],1,length(h.Settings.stim(1).patternvalue{6}))
%         repmat([1/(freq1*step^6*2)],1,length(h.Settings.stim(1).patternvalue{7})) 
%         repmat([1/(freq1*step^7*2)],1,length(h.Settings.stim(1).patternvalue{8})) 
%         repmat([1/(freq1*step^8*2)],1,length(h.Settings.stim(1).patternvalue{9}))
%         %repmat([1/(freq1*step^9*2)],1,length(h.Settings.stim(1).patternvalue{10})) 
%         };
%     % increase the duration of blank space at the end
%     for d = 1:length(h.Settings.stim(1).dur)
%         h.Settings.stim(1).dur{d}(end) = train_dur - sum(h.Settings.stim(1).dur{d}(1:end-1));
%     end
    
     % Freq stimuli for Chris Knaggs
    train_dur = 2; % seconds
    h.Settings.stim(1).patternvalue = {
        repmat([400 0],1,18)
        repmat([400 0],1,18) 
        repmat([400 0],1,19) 
        repmat([400 0],1,19) 
        repmat([400 0],1,21) 
        repmat([400 0],1,21) 
        repmat([400 0],1,22) 
        repmat([400 0],1,22) 
        repmat([400 0],1,23)
        repmat([400 0],1,23) 
        repmat([400 0],1,24) 
        repmat([400 0],1,24) 
        repmat([400 0],1,26) 
        repmat([400 0],1,26) 
        repmat([400 0],1,27) 
        repmat([400 0],1,27)  
        };
    h.Settings.stim(1).dur = {
        repmat([1/(18*2)],1,length(h.Settings.stim(1).patternvalue{1}))
        repmat([1/(18.5*2)],1,length(h.Settings.stim(1).patternvalue{2}))
        repmat([1/(19*2)],1,length(h.Settings.stim(1).patternvalue{3}))
        repmat([1/(19.5*2)],1,length(h.Settings.stim(1).patternvalue{4}))
        repmat([1/(20.5*2)],1,length(h.Settings.stim(1).patternvalue{5}))
        repmat([1/(21*2)],1,length(h.Settings.stim(1).patternvalue{6}))
        repmat([1/(21.5*2)],1,length(h.Settings.stim(1).patternvalue{7})) 
        repmat([1/(22*2)],1,length(h.Settings.stim(1).patternvalue{8})) 
        repmat([1/(23*2)],1,length(h.Settings.stim(1).patternvalue{9}))
        repmat([1/(23.5*2)],1,length(h.Settings.stim(1).patternvalue{10}))
        repmat([1/(24*2)],1,length(h.Settings.stim(1).patternvalue{11}))
        repmat([1/(24.5*2)],1,length(h.Settings.stim(1).patternvalue{12}))
        repmat([1/(25.5*2)],1,length(h.Settings.stim(1).patternvalue{13}))
        repmat([1/(26*2)],1,length(h.Settings.stim(1).patternvalue{14}))
        repmat([1/(26.5*2)],1,length(h.Settings.stim(1).patternvalue{15})) 
        repmat([1/(27*2)],1,length(h.Settings.stim(1).patternvalue{16})) 
        };
    % increase the duration of blank space at the end
    for d = 1:length(h.Settings.stim(1).dur)
        h.Settings.stim(1).dur{d}(end) = train_dur - sum(h.Settings.stim(1).dur{d}(1:end-1));
    end
    
    %% the rest
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
        %h.Settings.stim(1).patternvalue = repmat([400 0],1,ngaps); % one per stimdur in each cell; one cell per oddball value
        h.Settings.stim(1).patternmethod = 'pitch';% Pattern type method: intensity, pitch, duration. 
        h.Settings.stim(1).f0 = 400; % pitch
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
    h.Settings.ntrials = [16]; % if 3 stims per trial should multiples of 3 to allow every combination as programmed in oddballvalue. 
    
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
