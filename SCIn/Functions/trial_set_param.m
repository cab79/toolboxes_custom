function h = trial_set_param(h)

switch h.Settings.stim(h.sn).control
    
    case {'LJTick-DAQ','labjack','ptb_visual'}

        if ~isfield(h.Settings.stim(h.sn),'inten_diff')
            h.Settings.stim(h.sn).inten_diff = 0;
        end
        if h.seqtype.thresh
            if ~isfield(h.Settings.threshold,'stim')
                h.Settings.threshold.stim = 1;
            end
        end
        h.dur = h.Settings.stim(h.sn).dur; 
       
        % set initial values from GUI or settings if not already
        % defined
        set_inten_mean = 0;
        if ~isfield(h.stim,'inten_mean') || length(h.stim)<h.sn
            set_inten_mean = 1;
        elseif isempty(h.stim(h.sn).inten_mean)
            set_inten_mean = 1;
        end
        if set_inten_mean
            if ~isfield(h,'inten_mean_gui')
                h.stim(h.sn).inten_mean=0;
            else
                h.stim(h.sn).inten_mean = str2double(h.inten_mean_gui);
            end
            if ~any(h.stim(h.sn).inten_mean) && ~isempty(h.Settings.stim(h.sn).inten)
                h.stim(h.sn).inten_mean = h.Settings.stim(h.sn).inten;
            end
            h.stim(h.sn).inten_mean_start = h.stim(h.sn).inten_mean;
        end
        set_inten_diff = 0;
        if ~isfield(h.stim,'inten_diff') || length(h.stim)<h.sn
            set_inten_diff = 1;
        elseif isempty(h.stim(h.sn).inten_diff)
            set_inten_diff = 1;
        end
        if set_inten_diff
            if ~isfield(h,'inten_diff_gui')
                h.stim(h.sn).inten_diff=0;
            else
                h.stim(h.sn).inten_diff = str2double(h.inten_diff_gui);
            end
            if ~any(h.stim(h.sn).inten_diff) && ~isempty(h.Settings.stim(h.sn).inten_diff)
                h.stim(h.sn).inten_diff = h.Settings.stim(h.sn).inten_diff;
            end
            if ~isfield(h,'inten_diff_max_gui')
                h.stim(h.sn).inten_diff_max=0;
            else
                h.stim(h.sn).inten_diff_max = str2double(h.inten_diff_max_gui);
            end
            if ~any(h.stim(h.sn).inten_diff_max) && ~isempty(h.Settings.stim(h.sn).inten_diff_max)
                h.stim(h.sn).inten_diff_max = h.Settings.stim(h.sn).inten_diff_max;
            end
            h.stim(h.sn).inten_diff_start = h.stim(h.sn).inten_diff_max;
        end

        % modify according to procedure
        if h.seqtype.thresh && h.Settings.threshold.stim==h.sn
            if ~isfield(h,'s')
                h.stim(h.sn).inten = h.stim(h.sn).inten_mean+h.Settings.threshold.startinglevel;
            else
                h.stim(h.sn).inten = h.stim(h.sn).inten_mean+h.s.StimulusLevel;
            end
        % adaptive trial
        elseif h.seqtype.adapt && h.Settings.adaptive_general.stim==h.sn
            detect = find(strcmp(h.atypes,'detect'));
            discrim = find(strcmp(h.atypes,'discrim'));
            % if this trial is adaptive
            if ~isnan(h.Seq.adapttype(h.i))
                % update mean from adaptive procedure.
                % do this even if it's a discrim trial
                if ~isempty(detect) && isfield(h,'s') 
                    if isfield(h.s.a(detect),'StimulusLevel') % if a level has been determined
                        h.stim(h.sn).inten_mean = h.s.a(detect).StimulusLevel;
                        if isempty(h.stim(h.sn).inten_mean)
                            h.stim(h.sn).inten_mean = h.stim(h.sn).inten_mean_start;
                        end
                    end
                % or set the adaptive starting level otherwise
                elseif ~isempty(detect) && ~isfield(h,'s')
                    h.Settings.adaptive(detect).startinglevel = h.stim(h.sn).inten_mean;
                end

                % update diff from adaptive
                if ~isempty(discrim)
                    if isfield(h,'s') && length(h.s.a)>=discrim % if a level has been determined
                        h.stim(h.sn).inten_diff = h.s.a(discrim).StimulusLevel;
                    elseif ~isfield(h,'s') && ~(h.stim(h.sn).inten_diff == 0 || isempty(h.stim(h.sn).inten_diff))
                        h.Settings.adaptive(discrim).startinglevel = h.stim(h.sn).inten_diff;
                    end
                    if h.stim(h.sn).inten_diff == 0 || isempty(h.stim(h.sn).inten_diff) % if not set in GUI or settings
                        h.stim(h.sn).inten_diff = h.Settings.adaptive(discrim).startinglevel;
                    end
                end

                % only do this if it's a discrim trial, or a ratio detect
                % trial only
                if h.Seq.adapttype(h.i)==discrim || (h.Seq.adapttype(h.i)==detect && strcmp(h.Settings.adaptive(detect).subtype,'ratio'))
                    % calculate intensity
                    if h.Seq.signal(h.sn,h.i)==h.sn
                        h.stim(h.sn).inten = h.stim(h.sn).inten_mean + h.Settings.adaptive(discrim).stepdir * h.stim(h.sn).inten_diff / 2;
                    else
                        h.stim(h.sn).inten = h.stim(h.sn).inten_mean - h.Settings.adaptive(discrim).stepdir * h.stim(h.sn).inten_diff / 2;
                    end
                else
                    h.stim(h.sn).inten = h.stim(h.sn).inten_mean;
                end

            % if adaptive is part of the sequence, but not this trial
            elseif isnan(h.Seq.adapttype(h.i)) && h.Settings.adaptive_general.stim==h.sn
                detect_thresh =  find(h.out.adaptive(:,10)==detect);
                discrim_thresh =  find(h.out.adaptive(:,10)==discrim);
                if isempty(detect_thresh)
                    meanval = h.stim(h.sn).inten_mean;
                else
                    meanval = h.out.adaptive(detect_thresh(end),7);
                end
                if isempty(discrim_thresh)
                    diffval = h.stim(h.sn).inten_diff;
                else
                    diffval = h.out.adaptive(discrim_thresh(end),7);
                end
                if h.Seq.signal(h.sn,h.i)==1
                    h.stim(h.sn).inten = meanval - diffval / 2;
                else
                    h.stim(h.sn).inten = meanval + diffval / 2;
                end
            end
        elseif h.seqtype.adapt && ~any(h.Settings.adaptive_general.stim==h.sn)
            h.stim(h.sn).inten = h.stim(h.sn).inten_mean;
            
        % Otherwise, use sequence to determine intensity
        else
            if strcmp(h.Settings.oddballmethod,'intensity')
                if iscell(h.Settings.oddballvalue)
                    if size(h.Settings.oddballvalue,1)==1
                        h.stim(h.sn).inten = h.Settings.oddballvalue{h.Seq.signal(h.sn,h.tr)};
                    else
                        h.stim(h.sn).inten = h.Settings.oddballvalue{h.Seq.signal(h.sn,h.tr),:};
                    end
                else
                    h.stim(h.sn).inten = h.Settings.oddballvalue(h.Seq.signal(h.sn,h.tr),:);
                end
            elseif strcmp(h.Settings.oddballmethod,'index') && size(h.Seq.signal,1)>=h.sn
                if isfield(h.Settings,'AL') && any(ismember(h.Settings.AL.stimnums, h.sn)) && isempty(h.stim(h.sn).inten_diff_max)
                    h.stim(h.sn).inten = h.stim(h.sn).inten_mean - h.Settings.AL.intenstim(h.Seq.signal(h.sn,h.i)) * h.stim(h.sn).inten_diff/2;
                elseif isfield(h.Settings,'AL') && any(ismember(h.Settings.AL.stimnums, h.sn)) && ~isempty(h.stim(h.sn).inten_diff_max)
                    %h.stim(h.sn).inten = h.stim(h.sn).inten_mean - h.Settings.AL.intenstim(h.Seq.signal(h.sn,h.i)) * h.stim(h.sn).inten_diff/2;
                    % calculate intensity
                    if h.Seq.signal(h.sn,h.i)==1
                        h.stim(h.sn).inten = h.stim(h.sn).inten_mean - h.stim(h.sn).inten_diff / 2;
                    elseif h.Seq.signal(h.sn,h.i)==2
                        h.stim(h.sn).inten = h.stim(h.sn).inten_mean + h.stim(h.sn).inten_diff / 2;
                    elseif h.Seq.signal(h.sn,h.i)==3
                        h.stim(h.sn).inten = h.stim(h.sn).inten_mean - h.stim(h.sn).inten_diff_max / 2;
                    elseif h.Seq.signal(h.sn,h.i)==4
                        h.stim(h.sn).inten = h.stim(h.sn).inten_mean + h.stim(h.sn).inten_diff_max / 2;
                    end
                elseif isfield(h.Settings,'AL') && any(ismember(h.Settings.AL.stimnums, h.sn))
                    % calculate intensity
                    if h.Seq.signal(h.sn,h.i)==1
                        h.stim(h.sn).inten = h.stim(h.sn).inten_mean - h.stim(h.sn).inten_diff / 2;
                    else
                        h.stim(h.sn).inten = h.stim(h.sn).inten_mean + h.stim(h.sn).inten_diff / 2;
                    end
                elseif isfield(h.Settings,'PL') && isempty(h.stim(h.sn).inten_diff_max)
                    h.stim(h.sn).inten = h.stim(h.sn).inten_mean - h.Settings.PL.intenstim(h.Seq.signal(h.sn,h.i)) * h.stim(h.sn).inten_diff/2;
                elseif isfield(h.Settings,'PL') && ~isempty(h.stim(h.sn).inten_diff_max)
                    %h.stim(h.sn).inten = h.stim(h.sn).inten_mean - h.Settings.AL.intenstim(h.Seq.signal(h.sn,h.i)) * h.stim(h.sn).inten_diff/2;
                    % calculate intensity
                    if h.Seq.signal(h.sn,h.i)==1
                        h.stim(h.sn).inten = h.stim(h.sn).inten_mean - h.stim(h.sn).inten_diff / 2;
                    elseif h.Seq.signal(h.sn,h.i)==2
                        h.stim(h.sn).inten = h.stim(h.sn).inten_mean + h.stim(h.sn).inten_diff / 2;
                    elseif h.Seq.signal(h.sn,h.i)==3
                        h.stim(h.sn).inten = h.stim(h.sn).inten_diff_max;
                    elseif h.Seq.signal(h.sn,h.i)==4
                        h.stim(h.sn).inten = h.stim(h.sn).inten_mean + h.stim(h.sn).inten_diff / 2;
                    end
                elseif isfield(h.Settings,'PL')
                    % calculate intensity
                    if h.Seq.signal(h.sn,h.i)==1
                        h.stim(h.sn).inten = h.stim(h.sn).inten_mean - h.stim(h.sn).inten_diff / 2;
                    else
                        h.stim(h.sn).inten = h.stim(h.sn).inten_mean + h.stim(h.sn).inten_diff / 2;
                    end
                else
                    h.stim(h.sn).inten = h.Settings.stim(h.sn).inten;
                end
            end
        end

        % set max intensity
        if ~isfield(h.Settings.stim,'maxinten') || length(h.stim)<h.sn
            h.Settings.stim(h.sn).maxinten = inf;
        end
        h.stim(h.sn).inten = min(h.stim(h.sn).inten,h.Settings.stim(h.sn).maxinten); 
        

    case {'PsychPortAudio','audioplayer'}
    
        % create temporary variables for intensity, pitch and duration
        h.stim(h.sn).inten = h.Settings.stim(h.sn).inten;
        h.freq = h.Settings.stim(h.sn).f0;
        h.dur = h.Settings.stim(h.sn).dur; 
        h.varlevel = 0;
        % if calculating all trials, requires totdur:
        if strcmp(h.Settings.design,'continuous') && ~isfield(h,'totdur') % Use calculated duration by default.
            if isfield(h.Settings,'totdur')
                h.totdur = h.Settings.totdur; 
            end
        end
        if ~isfield(h.Settings.stim(h.sn),'inten_type')
            h.Settings.stim(h.sn).inten_type = 'abs';
        end
        if ~isfield(h.stim,'inten_diff') || length(h.stim)<h.sn || isempty(h.stim(h.sn).inten_diff)
            if ~isfield(h,'aud_diff_gui')
                h.stim(h.sn).inten_diff=0;
            else
                h.stim(h.sn).inten_diff = str2double(h.aud_diff_gui);
            end
            if ~any(h.stim(h.sn).inten_diff) && ~isempty(h.Settings.stim(h.sn).inten_diff)
                h.stim(h.sn).inten_diff = h.Settings.stim(h.sn).inten_diff;
            end
            h.stim(h.sn).inten_diff_start = h.stim(h.sn).inten_diff;
        end
        
        % thresholding
        if h.seqtype.thresh && ~h.seqtype.oddball
            if ~isfield(h,'s')
                h.varlevel = h.Settings.threshold.startinglevel;
            else
                h.varlevel = h.s.StimulusLevel;
            end
        end

        % condition method
        if isfield(h.Settings,'conditionmethod')
            if ~isempty(h.Settings.conditionmethod)
                if iscell(h.Settings.conditionmethod)
                    for i = 1:length(h.Settings.conditionmethod)
                        conditionmethod = h.Settings.conditionmethod{i};
                        if strcmp(conditionmethod,'pitch') || strcmp(conditionmethod,'freq')
                            if iscell(h.Settings.conditionvalue)
                                h.freq = h.Settings.conditionvalue{h.Seq.signal(h.sn,h.tr),i};
                            else
                                h.freq = h.Settings.conditionvalue(i,h.Seq.signal(h.sn,h.tr));
                            end
                        end
                        if strcmp(conditionmethod,'intensity')
                            if iscell(h.Settings.conditionvalue)
                                h.stim(h.sn).inten = h.Settings.conditionvalue{h.Seq.signal(h.sn,h.tr),i};
                            else
                                h.stim(h.sn).inten = h.Settings.conditionvalue(i,h.Seq.signal(h.sn,h.tr));
                            end
                        end
                        if strcmp(conditionmethod,'phase')
                            if iscell(h.Settings.conditionvalue)
                                h.alignphase = h.Settings.conditionvalue{h.Seq.signal(h.sn,h.tr),i};
                            else
                                h.alignphase = h.Settings.conditionvalue(i,h.Seq.signal(h.sn,h.tr));
                            end
                        end
                    end
                else
                    error('h.Settings.conditionmethod must be a cell');
                end
            end
        end

        % oddball method
        if h.seqtype.oddball && any(strcmp(h.Settings.oddballmethod,{'channel','intensity','index','pitch','duration','freq'}))
            if ~h.seqtype.adapt && ~h.seqtype.thresh
                if iscell(h.Settings.oddballvalue)
                    if size(h.Settings.oddballvalue,1)==1
                        oddval = h.Settings.oddballvalue{h.Seq.signal(h.sn,h.tr)};
                    else
                        oddval = h.Settings.oddballvalue{h.Seq.signal(h.sn,h.tr),:};
                    end
                else
                    oddval = h.Settings.oddballvalue(h.Seq.signal(h.sn,h.tr),:);
                end
                if strcmp(h.Settings.oddballmethod,'channel')
                    h.chan=oddval;
                elseif strcmp(h.Settings.oddballmethod,'intensity') && strcmp(h.Settings.stim(h.sn).inten_type,'abs') && h.stim(h.sn).inten_diff==0
                    h.stim(h.sn).inten = oddval;
                elseif strcmp(h.Settings.oddballmethod,'index') && strcmp(h.Settings.stim(h.sn).inten_type,'dB') && h.stim(h.sn).inten_diff~=0 && strcmp(h.Settings.stim(h.sn).patternmethod,'intensity') && isempty(h.Settings.stim(h.sn).patternvalue)
                    if isfield(h.Settings,'AL')
                        h.varlevel = h.Settings.AL.intenstim(h.Seq.signal(h.sn,h.i)) * h.stim(h.sn).inten_diff/2;
                    else
                        if h.Seq.signal(h.sn,h.i)==1
                            h.varlevel = h.stim(h.sn).inten_diff/2;
                        else
                            h.varlevel = -h.stim(h.sn).inten_diff/2;
                        end
                    end
                elseif strcmp(h.Settings.oddballmethod,'pitch') || strcmp(h.Settings.oddballmethod,'freq')
                    h.freq = oddval;
                elseif strcmp(h.Settings.oddballmethod,'duration')
                    h.dur = oddval;
                end
            elseif (h.seqtype.adapt && h.Settings.adaptive_general.stim==h.sn) || (h.seqtype.thresh && h.Settings.threshold.stim==h.sn)
                if isfield(h,'s')
                    h.varlevel = h.s.a(h.Seq.adapttype(h.i)).StimulusLevel;
                else
                    if h.seqtype.adapt
                        h.varlevel = h.Settings.adaptive.startinglevel;
                    else
                        h.varlevel = h.Settings.threshold.startinglevel;
                    end  
                end
                if strcmp(h.Settings.oddballmethod,'pitch') || strcmp(h.Settings.oddballmethod,'freq')
                    h.freq = [h.Settings.oddballvalue(1), (h.Settings.oddballvalue(1)+h.varlevel)]; % create new pitch pair
                    h.freq = h.freq(h.Seq.signal(h.sn,h.tr));
                % LEGACY CODE FROM HARRY'S EXPT: NO LONGER USE ABS FOR AUDIO
                %elseif strcmp(h.Settings.oddballmethod,'intensity') && strcmp(h.Settings.stim(h.sn).inten_type,'abs')
                %    h.stim(h.sn).inten = [h.Settings.oddballvalue(1), (h.Settings.oddballvalue(1)+h.varlevel)]; % create new pitch pair
                %    h.stim(h.sn).inten = h.stim(h.sn).inten(h.Seq.signal(h.sn,h.tr));
                elseif strcmp(h.Settings.oddballmethod,'duration') && (strcmp(h.Settings.stim(h.sn).patternmethod,'pitch') || strcmp(h.Settings.stim(h.sn).patternmethod,'freq'))
                    if iscell(h.Settings.oddballvalue)
                        h.dur = h.Settings.oddballvalue{h.Seq.signal(h.sn,h.tr),:};
                    else
                        h.dur = h.Settings.oddballvalue(h.Seq.signal(h.sn,h.tr),:);
                    end
                    h.freq = [h.Settings.stim(h.sn).patternvalue(1), (h.Settings.stim(h.sn).patternvalue(1)+h.varlevel)]; % create new pitch pair
                end
            end
        end

        %apply pitch pattern?
        h.trialtype.freqpattern=0;
        if isfield(h.Settings.stim(h.sn),'patternmethod')
            if strcmp(h.Settings.stim(h.sn).patternmethod,'pitch') || strcmp(h.Settings.stim(h.sn).patternmethod,'freq') % pitch changes
                h.trialtype.freqpattern=1;
                if ~((h.seqtype.adapt || h.seqtype.thresh) && (strcmp(h.Settings.oddballmethod,'pitch') || strcmp(h.Settings.oddballmethod,'freq'))) && ~(~isempty(strcmp(h.Settings.conditionmethod,'pitch')) || ~isempty(strcmp(h.Settings.conditionmethod,'freq'))) % pitch already defined above in this case
                    if strcmp(h.Settings.PL.oddballmethod,'pattern') && strcmp(h.Settings.stim(h.sn).durtype,'FromSequence')
                        %load from sequence file
                        global d
                        file = load(fullfile(d.seq, h.Settings.stim(h.sn).patternvalue));
        %                 fds = fieldnames(file.settings.stim(1));
        %                 for fd = 1:length(fds)
        %                     if strcmp(fds{fd},'durtype')
        %                         continue
        %                     end
        %                     h.Settings.stim(h.sn).(fds{fd}) = file.settings.stim(1).(fds{fd});
        %                 end
                        h.dur = file.settings.stim(1).dur;
                        h.freq = file.settings.stim(1).patternvalue{h.Seq.signal(h.sn,h.i)};
                    else
                        if isnumeric(h.Settings.stim(h.sn).patternvalue)
                            h.freq = h.Settings.stim(h.sn).patternvalue;
                        elseif iscell(h.Settings.stim(h.sn).patternvalue) %&& strcmp(h.Settings.oddballmethod,'pattern')
                            h.freq = h.Settings.stim(h.sn).patternvalue{h.Seq.signal(h.sn,h.i)};
                        end
                        % OLD CODE: NOT NEEDED?
                        %elseif iscell(h.Settings.stim(h.sn).patternvalue)
                        %    nDur = length(h.dur);
                        %    nPit = cellfun(@length,h.Settings.stim(h.sn).patternvalue);
                        %    h.freq = h.Settings.stim(h.sn).patternvalue{nPit==nDur};
                        %end
                    end
                end
            end
            
        end
        % apply response probe?
        if isfield(h.Settings,'RPmethod')
            if strcmp(h.Settings.RPmethod,'pitch') || strcmp(h.Settings.RPmethod,'freq')
                if h.Seq.RP(h.tr)==1
                    h.trialtype.freqpattern=1;
                    h.freq = h.Settings.RPvalue;
                    h.dur = h.Settings.RPdur;
                end
            end
        end
        %apply intensity pattern?
        h.trialtype.intenpattern=0;
        if isfield(h.Settings.stim(h.sn),'patternmethod')
            if strcmp(h.Settings.stim(h.sn).patternmethod,'intensity') && ~isempty(h.Settings.stim(h.sn).patternvalue) % intensity changes
                h.trialtype.intenpattern=1;
                if ~any(strcmp(h.Settings.conditionmethod,'intensity')) % then already defined
                    h.stim(h.sn).inten = h.Settings.stim(h.sn).patternvalue;
                end
            end
        end
        % apply response probe?
        h.resp_probe=0;
        if isfield(h.Settings,'RPmethod')
            if strcmp(h.Settings.RPmethod,'intensity')
                if h.Seq.RP(h.tr)==1
                    h.resp_probe=1;
                    h.stim(h.sn).inten = h.Settings.RPvalue;
                    h.dur = h.Settings.RPdur;
                end
            end
        end
        
end