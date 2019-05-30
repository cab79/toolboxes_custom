function h=PTBvisual(h)

% Clear the screen
%sca;
%clear all
%close all;
%clearvars;
%% setup
% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

% Seed the random number generator. Here we use the an older way to be
% compatible with older systems. Newer syntax would be rng('shuffle').
% For help see: help rand
%rand('seed', sum(100 * clock));

% Get the screen numbers. This gives us a number for each of the screens
% attached to our computer. For help see: Screen Screens?
if isunix
    Screen('Preference','SkipSyncTests', 0);
else
    Screen('Preference','SkipSyncTests', 1);
end
h.screens = Screen('Screens');

% Draw we select the maximum of these numbers. So in a situation where we
% have two screens attached to our monitor we will draw to the external
% screen. When only one screen is attached to the monitor we will draw to
% this. For help see: help max
h.screenNumber = max(h.screens);

% Define black and white (white will be 1 and black 0). This is because
% luminace values are (in general) defined between 0 and 1.
% For help see: help WhiteIndex and help BlackIndex
h.white = WhiteIndex(h.screenNumber);
h.black = BlackIndex(h.screenNumber);

% Open an on screen window and color it black.
% For help see: Screen OpenWindow?
[h.window, h.windowRect] = PsychImaging('OpenWindow', h.screenNumber, h.black);

% Get the size of the on screen window in pixels.
% For help see: Screen WindowSize?
[h.screenXpixels, h.screenYpixels] = Screen('WindowSize', h.window);

% Get the centre coordinate of the window in pixels
% For help see: help RectCenter
[h.xCenter, h.yCenter] = RectCenter(h.windowRect);

% Enable alpha blending for anti-aliasing
% For help see: Screen BlendFunction?
% Also see: Chapter 6 of the OpenGL programming guide
Screen('BlendFunction', h.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

