function ptb = bat_initialize_ptb(params)

% Itializes PsychToolbox graphics and audio for BAT experiment. Stores
% graphics and audio handles in structure array called ptb.
% 
% INPUT:
%     params: parameter structure array with the following fields:
%          * freq - frequency at which audio should be played
%          * nrchannels - 1 (mono) or 2 (stereo) playback
% OUTPUT:
%     ptb: structure array with PsychToolbox object handles that will be
%     needed in the experiment function
%
% Written 6 OCT 2016 by Brian K. Hurley
%   - Some Psychtoolbox code adapted from BAT_adaptive_ibi.m &
%       BAT_adaptive_phase.m written by Jessica M. Ross

%% Graphics Initialization

Screen('Preference', 'VisualDebugLevel', 1);
Screen('Preference', 'SkipSyncTests', 1);

%HideCursor;
PsychDefaultSetup(2)

%Set screen number to my one screen (this can detect multiple displays if
%needed)
screenNum = max(Screen('Screens'));

%Define colors
ptb.white = WhiteIndex(screenNum);
black = BlackIndex(screenNum);

%% Instruction Screen

%Black screen appears first
[ptb.window, ~] = PsychImaging('OpenWindow', screenNum, black);
Screen('Flip',ptb.window);

% Deliver instruction text
Screen('TextSize',ptb.window, 30);
Screen('TextFont',ptb.window,'Ariel');
Screen('TextStyle', ptb.window, 1);
textString = ['Please listen to each musical excerpt and judge whether ' ...
  'the click track is on-beat or off-beat. \n\n\n Respond at the end of ' ... 
  'each musical excerpt by pressing Q if the tone is off-beat or P if ' ... 
  'the tone is on-beat.'];
DrawFormattedText(ptb.window, textString, 'center', 'center', ptb.white, 60);

Screen('TextSize',ptb.window, 15);
textString = '< Press any key to continue >';
DrawFormattedText(ptb.window, textString, 'center', 1000, ptb.white);

Screen('Flip', ptb.window);
KbStrokeWait;

% subject will see abbreviated on-screen instructions throughout task
Screen('TextSize',ptb.window, 30);
trial_textString = 'When music ends\n\n\nPress Q for off-beat, P for on-beat';
DrawFormattedText(ptb.window, trial_textString, 'center', 'center', ptb.white);
Screen('Flip', ptb.window);
WaitSecs(2)

%% Audio Initialization

% try 
%   PsychPortAudio('Close');
% end

% Initialize sound driver. If driver open from previous session, will be
% blocked for current session and will throw error. In that case, close 
% PsychPortAudio and re-initialize driver
try
  InitializePsychSound(1);
catch
  PsychPortAudio('Close');
  InitializePsychSound(1);
end

% Open Psych-Audio port, with the follow arguements
% (1) [] = default sound device
% (2) 1 = sound playback only
% (3) 1 = default level of latency
% (4) Requested frequency in samples per second
% (5) 2 = stereo putput
ptb.pahandle = PsychPortAudio('Open', [], 1, 1, params.freq, params.nrchannels);

% Keyboard Information
ptb.escapeKey = KbName('ESCAPE');
ptb.NO = KbName('q'); % Q key = NO
ptb.YES = KbName('p'); % P key = YES

% Start immediately (0 = immediately)
ptb.startCue = 0;

% Should we wait for the device to really start (1 = yes)
% INFO: See help PsychPortAudio
ptb.waitForDeviceStart = 1;