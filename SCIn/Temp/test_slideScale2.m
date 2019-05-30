% NEED TO:
% Add in code to start with a scanner trigger (KB button 's')
% Integrate into SCIn, including within trial design

% Preliminary stuff
% Clear Matlab/Octave window:
clc;
PsychDefaultSetup(2);

% check for Opengl compatibility, abort otherwise:
AssertOpenGL;

bgColor   = [128 128 128];

HideCursor;

% Get the screen numbers. This gives us a number for each of the screens
% attached to our computer. For help see: Screen Screens?
%Screen('Preference','SkipSyncTests', 1);
screens = Screen('Screens');

% Draw we select the maximum of these numbers. So in a situation where we
% have two screens attached to our monitor we will draw to the external
% screen. When only one screen is attached to the monitor we will draw to
% this. For help see: help max
screenNumber = max(screens);

% Get information about the screen and set general things
Screen('Preference', 'SuppressAllWarnings',0);
Screen('Preference', 'SkipSyncTests', 1);

rect          = Screen('Rect',0);
screenRatio   = rect(3)/rect(4);
pixelSizes    = Screen('PixelSizes', 0);
startPosition = round([rect(3)/2, rect(4)/2]);

% Creating screen etc.
[myScreen, rect] = Screen('OpenWindow', screenNumber, bgColor);
center           = round([rect(3) rect(4)]/2);


question  = 'How intense is your pain?';
endPoints = {'no pain', 'extreme pain'};

[position, RT, answer,out] = slideScale(...
myScreen, question, rect, endPoints, ...
    'device', 'keyboard', ...
    'scalaposition', 0.5, ...
    'startposition', 'center', ...
    'displayposition', true, ...
    'range', 2, ...
    'incr', 5 ... 
);

Screen('CloseAll') 