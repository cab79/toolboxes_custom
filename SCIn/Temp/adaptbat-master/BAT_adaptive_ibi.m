function[] = BAT_adaptive_ibi(subject)

% Adaptive inter-beat interval version of the Beat Alignment Test, version 2
% Written by Jessica Ross, March 2016
%
% Iversen, JR, Patel, AD. The Beat Alignment Test (BAT): Surveying beat processing
% abilities in the general population. Proceedings of the 10th International Conference
% on Music Perception and Cognition 2008: 465-468.
%
% Discriminate between correct ibi and augmented ibi. The 12 stimuli from
% the BAT are used. Trial 1 and 2 have correct ibi's. Trial 3 has an ibi that
% has been augmented by 10%. For each trial after trial 3, an incorrect
% response will result in the next trial being easier (using a bigger ibi).
% 2 correct responses in a row will result in the next trial being harder
% (using a smaller ibi). A correct response followed by an incorrect
% response will result in the next trial using the same ibi as the previous
% trial, but with a new musical stimulus. Lengthened and shortened ibi's
% are alternated.
%
% INPUT:
%     subject - subject ID string (e.g., 'JR1')
%
% OUTPUT FILE:
%     Column 1 = subject's response (off beat = 0, on beat = 1)
%     Column 2 = augmentation condition (off beat = 0, on beat = 1)
%     Column 3 = augmentation amount (0%, -10%, 10%, etc...)
%     Column 4 = correct/incorrect (incorrect response = 1, correct response = 2)

Screen('Preference', 'VisualDebugLevel', 1);
Screen('Preference', 'SkipSyncTests', 1);
HideCursor;
PsychDefaultSetup(2);

%Set screen number to my one screen (this can detect multiple displays if
%needed)
screenNum = max(Screen('Screens'));

%Define colors
white = WhiteIndex(screenNum);
grey = white / 2;
black = BlackIndex(screenNum);

% Keyboard Information
escapeKey = KbName('ESCAPE');
NO = KbName('q'); % No (off beat) is q during experiment, and 0 in output file
YES = KbName('p'); % Yes (on beat) is p during experiment, and 1 in output file
downKey = KbName('DownArrow');

%instruction_screen____________________________
%Black screen appears first
[window, rect] = PsychImaging('OpenWindow', screenNum, black);
Screen('Flip',window);

%Query frame duration
ifi = Screen('GetFlipInterval',window);

%Query maximum priority level
topPriorityLevel = MaxPriority(window);

%Get the center coordinate of the window
[X,Y] = RectCenter(rect);

%Set text properties
Screen('TextSize',window, 30);
Screen('TextFont',window,'Ariel');
Screen('TextStyle', window, 1);
textString = ['Welcome to the experiment. \n\n\n\n\n Please listen to each' ...
    'musical excerpt and judge whether the \n click track is on beat or off ' ...
    'beat. \n\n Respond at the end of each musical excerpt by \n pressing ' ...
    'q if the tone is off beat \n or p if the tone is on beat. \n\n This ' ...
    'should take about 7 minutes.'];
DrawFormattedText(window, textString, 'center', 'center', white);

Screen('TextSize',window, 15);
textString = ['\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n'...
    '< Press any key to continue >'];
DrawFormattedText(window, textString, 'center', 'center', white);

Screen('Flip', window);
KbStrokeWait;
tic
%_____________________________

responses = zeros(26,4);

% screenNum=0;
% [window, rect] = Screen('OpenWindow', screenNum, 1);
% [X,Y] = RectCenter(rect);
FixCross = [X-1,Y-40,X+1,Y+40;X-40,Y-1,X+40,Y+1];
Screen('FillRect', window, [255,255,255], FixCross');
Screen('Flip', window);
WaitSecs(2)

% Initialize Sounddriver
InitializePsychSound(1);

% Number of channels and Frequency of the sound
nrchannels = 2;
freq = 48000;

% How many times we wish to play the sound
repetitions = 1;

% Start immediately (0 = immediately)
startCue = 0;

% Should we wait for the device to really start (1 = yes)
% INFO: See help PsychPortAudio
waitForDeviceStart = 1;

% Open Psych-Audio port, with the follow arguements
% (1) [] = default sound device
% (2) 1 = sound playback only
% (3) 1 = default level of latency
% (4) Requested frequency in samples per second
% (5) 2 = stereo putput
pahandle = PsychPortAudio('Open', [], 1, 1, freq, nrchannels);

% Set the volume to half for this demo
PsychPortAudio('Volume', pahandle, 0.5);

% Get wav file names and categorize them by augmentation interval
data_loc =  'stimuli_adaptive_short/';
on0 = dir([data_loc,'*_B0_v2.wav']);dirlist_on0 = {on0.name};
offn1 = dir([data_loc,'*_B-1_v2.wav']);dirlist_offn1 = {offn1.name};
offn2 = dir([data_loc,'*_B-2_v2.wav']);dirlist_offn2 = {offn2.name};
offn3 = dir([data_loc,'*_B-3_v2.wav']);dirlist_offn3 = {offn3.name};
offn4 = dir([data_loc,'*_B-4_v2.wav']);dirlist_offn4 = {offn4.name};
offn5 = dir([data_loc,'*_B-5_v2.wav']);dirlist_offn5 = {offn5.name};
offn6 = dir([data_loc,'*_B-6_v2.wav']);dirlist_offn6 = {offn6.name};
offn7 = dir([data_loc,'*_B-7_v2.wav']);dirlist_offn7 = {offn7.name};
offn8 = dir([data_loc,'*_B-8_v2.wav']);dirlist_offn8 = {offn8.name};
offn9 = dir([data_loc,'*_B-9_v2.wav']);dirlist_offn9 = {offn9.name};
offn10 = dir([data_loc,'*_B-10_v2.wav']);dirlist_offn10 = {offn10.name};
offp1 = dir([data_loc,'*_B1_v2.wav']);dirlist_offp1 = {offp1.name};
offp2 = dir([data_loc,'*_B2_v2.wav']);dirlist_offp2 = {offp2.name};
offp3 = dir([data_loc,'*_B3_v2.wav']);dirlist_offp3 = {offp3.name};
offp4 = dir([data_loc,'*_B4_v2.wav']);dirlist_offp4 = {offp4.name};
offp5 = dir([data_loc,'*_B5_v2.wav']);dirlist_offp5 = {offp5.name};
offp6 = dir([data_loc,'*_B6_v2.wav']);dirlist_offp6 = {offp6.name};
offp7 = dir([data_loc,'*_B7_v2.wav']);dirlist_offp7 = {offp7.name};
offp8 = dir([data_loc,'*_B8_v2.wav']);dirlist_offp8 = {offp8.name};
offp9 = dir([data_loc,'*_B9_v2.wav']);dirlist_offp9 = {offp9.name};
offp10 = dir([data_loc,'*_B10_v2.wav']);dirlist_offp10 = {offp10.name};

wavfilename = [];

% no need to hard-code neg & pos values. change to abs val staircase and
% randomly select pos or neg - BH
aug_intervals = [on0 offn1 offn2 offn3 offn4 offn5 offn6 offn7 offn8 offn9 ...
    offn10 offp1 offp2 offp3 offp4 offp5 offp6 offp7 offp8 offp9 offp10];
aug_cols = [1,21,11,20,10,19,9,18,8,1,1,17,17,7,7,16,16,6,6,15,15,1,1,5,5, ...
    14,14,4,4,1,1,13,13,3,3,12,12,2,2,1,1,1,1,1,1,1,1,1,1,1,1];
aug_amount = [0 10 -10 9 -9 8 -8 7 -7 0 0 6 6 -6 -6 5 5 -5 -5 4 4 0 0 -4 -4 ...
    3 3 -3 -3 0 0 2 2 -2 -2 1 1 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0];

trial_order = {};
k = randperm(12);

% This should be 1 loop that performs these different operations
% conditionally depending on trial number - BH

for trial_number = 1; % why is this a loop when it's hard-coded as 1 iteration? - BH
    %Get stimulus and make sure it isn't a repeat
    stim = aug_intervals(k(trial_number),1).name;
    song = stim(1:3);
    trial_order{trial_number} = char(song)
    
    %Start with an on beat stimulus
    wavfilename = [data_loc aug_intervals(k(1),1).name];
    aug_cond = 1;
    aug_trial = 0;
    
    % Read WAV file from filesystem:
    [y, freq] = audioread(wavfilename);
    wavdata = y';
    nrchannels = size(wavdata,1); % Number of rows == number of channels.
    
    % Wait 0.5 s before starting
    WaitSecs(0.5);
    
    % Fill the audio playback buffer with the audio data
    PsychPortAudio('FillBuffer', pahandle, [wavdata; wavdata]);
    
    % Start audio playback
    PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);
    
    % Collect keyboard response: on beat/YES = p, off beat/NO = q
    respToBeMade = true;
    while respToBeMade == true
        
        [keyIsDown,secs, keyCode] = KbCheck;
        if keyCode(YES)
            response = 1;
            respToBeMade = false;
        elseif keyCode(NO)
            response = 0;
            respToBeMade = false;
        elseif keyCode(escapeKey)
            sca;
        end
        
    end
    
    responses(trial_number,1) = response;
    responses(trial_number,2) = aug_cond;
    responses(trial_number,3) = aug_trial;
    
    % 1 is incorrect, 2 is correct
    if response == 1
        score = 2;
    elseif response == 0
        score = 1;
    end
    
    responses(trial_number,4) = score;
    
end
x = 2;
last = responses((trial_number),4);
if last == 1 % 1=incorrect
    aug_col = aug_cols(x);
    aug_trial = aug_amount(x);
elseif last == 2 % 2=correct
    %         previous = responses((trial_number-1),4);
    %         if previous == 2
    aug_col = aug_cols(x);
    aug_trial = aug_amount(x);
else aug_col = aug_cols(x-1);
    aug_trial = aug_amount(x-1);
    %         end
end

for trial_number = 2;
    k = randperm(12);
    
    %Get stimulus and make sure it isn't a repeat
    if trial_order{trial_number-1} == aug_intervals(k(1),aug_col).name(1:3)
        stim = aug_intervals(k(2),aug_col).name;
    else
        stim = aug_intervals(k(1),aug_col).name;
    end
    song = stim(1:3);
    trial_order{trial_number} = char(song)
    
    %Start with an on beat stimulus
    wavfilename = [data_loc aug_intervals(k(1),1).name];
    aug_cond = 1;
    aug_trial = 0;
    
    % Read WAV file from filesystem:
    [y, freq] = audioread(wavfilename);
    wavdata = y';
    nrchannels = size(wavdata,1); % Number of rows == number of channels.
    
    % Wait 0.5 s before starting
    WaitSecs(0.5);
    
    % Fill the audio playback buffer with the audio data
    PsychPortAudio('FillBuffer', pahandle, [wavdata; wavdata]);
    
    % Start audio playback
    PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);
    
    % Collect keyboard response: on beat/YES = p, off beat/NO = q
    respToBeMade = true;
    while respToBeMade == true
        
        [keyIsDown,secs, keyCode] = KbCheck;
        if keyCode(YES)
            response = 1;
            respToBeMade = false;
        elseif keyCode(NO)
            response = 0;
            respToBeMade = false;
        elseif keyCode(escapeKey)
            sca;
        end
        
    end
    
    responses(trial_number,1) = response;
    responses(trial_number,2) = aug_cond;
    responses(trial_number,3) = aug_trial;
    
    % 1 is incorrect, 2 is correct
    if response == 1
        score = 2;
    elseif response == 0
        score = 1;
    end
    
    responses(trial_number,4) = score;
    
end
x = 2;
last = responses((trial_number),4);
if last == 1 % 1=incorrect
    aug_col = aug_cols(x);
    aug_trial = aug_amount(x);
elseif last == 2 % 2=correct
    previous = responses((trial_number-1),4);
    if previous == 2
        aug_col = aug_cols(x);
        aug_trial = aug_amount(x);
    else aug_col = aug_cols(x-1);
        aug_trial = aug_amount(x-1);
    end
end

for trial_number = 3:26;
    k = randperm(12);
    
    %Get stimulus and make sure it isn't a repeat
    if trial_order{trial_number-1} == aug_intervals(k(1),aug_col).name(1:3)
        stim = aug_intervals(k(2),aug_col).name;
    else
        stim = aug_intervals(k(1),aug_col).name;
    end
    song = stim(1:3);
    trial_order{trial_number} = char(song)
    
    %Now start the adaptive part, beginning with 10% augmented ibi (column 21)
    wavfilename = [data_loc aug_intervals(k(1),aug_col).name];
    if aug_col == 1
        aug_cond = 1;
        
    else
        aug_cond = 0;
    end
    
    % Read WAV file from filesystem:
    [y, freq] = audioread(wavfilename);
    wavdata = y';
    nrchannels = size(wavdata,1); % Number of rows == number of channels.
    
    % Wait 0.5 s before starting
    WaitSecs(0.5);
    
    % Fill the audio playback buffer with the audio data
    PsychPortAudio('FillBuffer', pahandle, [wavdata; wavdata]);
    
    % Start audio playback
    PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);
    
    % Collect keyboard response: on beat/YES = p, off beat/NO = q
    respToBeMade = true;
    if aug_col == 1
        while respToBeMade == true
            [keyIsDown,secs, keyCode] = KbCheck;
            if keyCode(YES)
                response = 1;
                respToBeMade = false;
            elseif keyCode(NO)
                response = 0;
                respToBeMade = false;
            elseif keyCode(escapeKey)
                sca;
            end
        end
        
        responses(trial_number,1) = response;
        responses(trial_number,2) = aug_cond;
        responses(trial_number,3) = aug_trial;
        
        % 1 is incorrect, 2 is correct
        if response == 1
            score = 2;
        elseif response == 0
            score = 1;
        end
        
        responses(trial_number,4) = score;
        
        last = responses((trial_number),4);
        if last == 1 % 1=incorrect
            x = x-1;
            aug_col = aug_cols(x);
            aug_trial = aug_amount(x);
        elseif last == 2 % 2=correct
            previous = responses((trial_number-1),4);
            if previous == 2
                x = x + 2;
                aug_col = aug_cols(x);
                aug_trial = aug_amount(x);
            else x = x + 1;
                aug_col = aug_cols(x);
                aug_trial = aug_amount(x);
            end
        end
    else
        while respToBeMade == true
            [keyIsDown,secs, keyCode] = KbCheck;
            if keyCode(YES)
                response = 1;
                respToBeMade = false;
            elseif keyCode(NO)
                response = 0;
                respToBeMade = false;
            elseif keyCode(escapeKey)
                sca;
            end
        end
        
        responses(trial_number,1) = response;
        responses(trial_number,2) = aug_cond;
        responses(trial_number,3) = aug_trial;
        
        % 1 is incorrect, 2 is correct
        if response == 1
            score = 1;
        elseif response == 0
            score = 2;
        end
        
        responses(trial_number,4) = score;
        
        last = responses((trial_number),4);
        if last == 1 % 1=incorrect
            x = x-1;
            aug_col = aug_cols(x);
            aug_trial = aug_amount(x);
        elseif last == 2 % 2=correct
            previous = responses((trial_number-1),4);
            if previous == 2
                x = x + 2;
                aug_col = aug_cols(x);
                aug_trial = aug_amount(x);
            else x = x + 1;
                aug_col = aug_cols(x);
                aug_trial = aug_amount(x);
            end
        end
    end
    
end

sca;

ResultsFolder ='BAT_adaptive_results/';
Outputfile = [ResultsFolder subject '_ibi_results.csv'];
csvwrite(Outputfile,responses);

elapsedTime = toc

%Open summary file for output
Summaryfile = [ResultsFolder subject '_summary.txt'];
fid = fopen(Summaryfile,'a+');
for sf = 1:length(trial_order)
    fprintf(fid,'\n %s',trial_order{sf});
end
fprintf(fid,'\n %f', elapsedTime);
close all;
