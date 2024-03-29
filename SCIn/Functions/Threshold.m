function h = Threshold(h)

if ~isfield(h.Settings,'threshold')
    return
end

if isfield(h,'s')
    s=h.s;
else
    s = ThresholdParameters(h); % see function at the end
end

if ~isfield(s,'lasttrialrun')
    s.lasttrialrun = 0;
end
if ~isfield(s,'lastpresstrial')
    s.lastpresstrial = 0;
end


%% RUN


%if strcmp(opt,'responded')
    
%else
%    s.lasttrialrun = h.i;
%end
    
% evaluate the subject's response
resi=0;
if ~isempty(h.out.presstrial)
    % only run once per trial if pressed
    if ~(s.lastpresstrial == h.out.presstrial(end))
        s.lastpresstrial = h.out.presstrial(end);
    else
        % only run once per trial if NOT pressed
        if ~(s.lasttrialrun == h.i)
            s.lasttrialrun = h.i;
        else
            return
        end
    end
    if h.out.presstrial(end) == h.i
        %presstrial=trls_pressed(end);
        pressbutton = h.out.pressbutton(end);
        resi = find(strcmp(pressbutton(end),h.Settings.buttonopt)); % which button was pressed?
        resfun = h.Settings.threshold.signalval(resi); %what is meaning of this response?
    end
else
    % only run once per trial if NOT pressed
    if ~(s.lasttrialrun == h.i)
        s.lasttrialrun = h.i;
    else
        return
    end
end
if isempty(resi)
    return
end
disp('Running Threshold')
if resi && resfun==1
%s.SubjectAccuracy(s.adaptive.trial)= EvaluateAnswer(CorrectAnswer,s.feedback,Question);   %evaluing the subject answer (right or wrong)
    s.change= -1;
elseif resi && resfun==2
    s.change = 0;
else
    s.change= 1;
end

% UPDATE THE ROWOFOUTPUT
s.rowofoutput (1, 1) = s.trial;
s.rowofoutput (1, 2) = s.StimulusLevel;
s.rowofoutput (1, 3) = s.actStimulusLevel;
s.rowofoutput (1, 4) = s.change;
% here I update the count for the up-s.down motion
if s.change == 1
    if s.isstep == 1
        s.StimulusLevel = s.StimulusLevel + s.actualstep;
    else
        s.StimulusLevel = s.StimulusLevel * s.actualstep;
    end
elseif s.change == -1
    s.StimulusLevel = h.Settings.threshold.startinglevel;
end

% additional attenutation from GUI or settings
if strcmp(h.Settings.stim(h.sn).control,'PsychPortAudio') || strcmp(h.Settings.stim(h.sn).control,'audioplayer')
    if isfield(h,'vol_atten')
        try
            inten_atten = str2double(get(h.vol_atten,'string'));
        catch
            inten_atten = str2double(h.vol_atten);
        end
    else
        inten_atten = h.Settings.atten; 
    end
    s.actStimulusLevel = min(h.Settings.threshold.maxinten,inten_atten+s.StimulusLevel);
end
if strcmp(h.Settings.stim(h.sn).control,'labjack') || strcmp(h.Settings.stim(h.sn).control,'LJTick-DAQ') || strcmp(h.Settings.stim(h.sn).control,'T7')
    % additional intensity (e.g. GUI control of base tactile intensity)
    if isfield(h,'baseinten')
        baseinten = h.baseinten(1);
    else
        baseinten = 0; 
    end
    s.actStimulusLevel = min(h.Settings.threshold.maxinten,baseinten+s.StimulusLevel);
end

% UPDATE THE ROWOFOUTPUT
s.trial = s.trial + 1;

%disp(['length_blockthresh = ' num2str(length(s.blockthresholds))]);
%disp(['Next trial level = ' num2str(s.actStimulusLevel)]);

% threshold calculation
s.sensthresholds=0;
s.painthresholds=0;
if ~isempty(s.out.threshold)
    
    % sensory threshold
    presstrialsall = find(s.out.threshold(:,4)==0);
    if any(presstrialsall)
        s.sensthresholds = s.out.threshold(presstrialsall(end)-1,3);
        fprintf('Sensory threshold on this trial equal to %1.3f\n', s.sensthresholds);
    end
    
    % pain threshold
    presstrialsall = find(s.out.threshold(:,4)==-1);
    if any(presstrialsall)
        s.painthresholds = s.out.threshold(presstrialsall(end)-1,3);
        fprintf('Pain threshold on this trial equal to %1.3f\n', s.painthresholds);
    end
end
s.rowofoutput (1, 5) = s.sensthresholds;
s.rowofoutput (1, 6) = s.painthresholds;

% UPDATE THE GLOBAL OUTPUT VARIABLE
s.out.threshold = [s.out.threshold; s.rowofoutput];

presstrialsall = find(s.out.threshold(:,4)==0);
if length(presstrialsall)>2 && presstrialsall(end)<s.trial-2
    try
        fprintf('Sensory Threshold over last 3 trials equal to %1.3f\n', nanmean(s.out.threshold(presstrialsall(end-2:end)+1,5)));
    catch
        fprintf('Sensory Threshold over last 3 trials equal to %1.3f\n', mean(s.out.threshold(presstrialsall(end-2:end)+1,5)));
    end
end
presstrialsall = find(s.out.threshold(:,4)==-1);
if length(presstrialsall)>2 && presstrialsall(end)<s.trial-2
    try
        fprintf('Pain Threshold over last 3 trials equal to %1.3f\n', nanmean(s.out.threshold(presstrialsall(end-2:end)+1,6)));
    catch
        fprintf('Pain Threshold over last 3 trials equal to %1.3f\n', mean(s.out.threshold(presstrialsall(end-2:end)+1,6)));
    end
end

h.s =s;
h.out.threshold = s.out.threshold;
%fprintf('Press return to continue\n');
%pause
%fprintf ('\nBLOCK ENDED\n');
%pause(2)
end

function s = ThresholdParameters(h)

%s.stepsize = 2;
%s.SaveResults =1;
s.numForthresh = 3;
s.isstep = 1;


% if this is the first run, do some setup
if ~isfield(s,'trial')
    s.out.threshold = [];
    s.rowofoutput = zeros(1, 6);
    s.expthresholds = zeros(1, 1);
    s.trial = 1;
    s.blockthresholds = zeros(length(s.numForthresh), 1);
    s.StimulusLevel = h.Settings.threshold.startinglevel;
    s.actStimulusLevel = NaN;
    s.actualstep = h.Settings.threshold.step;
end

end