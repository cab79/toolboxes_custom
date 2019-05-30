function h=VisualStim(h,opt)
Screen('Preference','SkipSyncTests', 1);

switch opt
    case 'blank'
        % Blank screen: Flip again to sync us to the vertical retrace at the same time as
        % drawing our fixation point
        Screen('DrawDots', h.window, [h.xCenter; h.yCenter], 10, h.black, [], 2);
        h.vbl = Screen('Flip', h.window,0,0,1);
        
    case 'stim'
        
        if isfield(h.Settings.stim(h.sn),'rectColor') && ~isempty(h.Settings.stim(h.sn).rectColor)
            % present a circular shape, e.g. for MoNoPly study
        
            for nr = 1:size(h.Settings.stim(h.sn).rectColor,1)

                if (h.seqtype.thresh && h.Settings.threshold.stim==h.sn && h.Settings.threshold.stimpart==nr)... %threshold
                        || (h.seqtype.adapt && h.Settings.adaptive_general.stim==h.sn && h.Settings.adaptive_general.stimpart==nr)... %adaptive
                        || (isfield(h.Settings,'AL') && any(ismember(h.Settings.AL.stimnums, h.sn)) && h.Settings.AL.stimpart==nr) %associative
                    h.inten_out = h.stim(h.sn).inten;%h.varlevel;
                else
                    h.inten_out = h.Settings.stim(h.sn).rectInt(nr);
                end

                % Set the color of our dot to full red. Color is defined by red green
                % and blue components (RGB). So we have three numbers which
                % define our RGB values. The maximum number for each is 1 and the minimum
                % 0. So, "full red" is [1 0 0]. "Full green" [0 1 0] and "full blue" [0 0
                % 1]. Play around with these numbers and see the result.
                rectColor = h.Settings.stim(h.sn).rectColor(nr,:) * h.inten_out;

                % size in pixels
                rectSize = h.Settings.stim(h.sn).rectSize(nr,:);

                % For Ovals we set a miximum diameter up to which it is perfect for
                maxDiameter = max(rectSize) * 1.01;

                % Center the rectangle on the centre of the screen
                centeredRect = CenterRectOnPointd(rectSize, h.xCenter, h.yCenter);

                % Draw the rect to the screen
                Screen('FillOval', h.window, rectColor, centeredRect, maxDiameter);

                %fix_size = [0 0 20 20];
                %fixation = CenterRectOnPointd(fix_size, h.xCenter, h.yCenter);
                %maxDiameterfixation = max(fix_size) * 1.01;
                %Screen('FillOval', h.window, [1 1 1], fixation, maxDiameterfixation);
            end
        elseif isfield(h.Settings.stim(h.sn),'textColor') && ~isempty(h.Settings.stim(h.sn).textColor)
            % present text, e.g. for PAL study
            
            % which outcome?
            text = h.Settings.stim(h.sn).text{h.Seq.signal(h.sn,h.i)};
            textColor = h.Settings.stim(h.sn).textColor{h.Seq.signal(h.sn,h.i)};
            
            % Draw text in the middle of the screen in Courier in white
            Screen('TextSize', h.window, h.Settings.stim(h.sn).textSize);
            Screen('TextFont', h.window, h.Settings.stim(h.sn).textFont);
            DrawFormattedText(h.window, text, 'center', 'center', textColor);
            
        end

        % Flip to the screen. This command basically draws all of our previous
        % commands onto the screen. See later demos in the animation section on more
        % timing details. And how to demos in this section on how to draw multiple
        % rects at once.
        % For help see: Screen Flip?
        h.vbl = Screen('Flip', h.window,0,0,1);
        %WaitSecs(1)
        %Screen('Close')
        
end
