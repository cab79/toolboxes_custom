function h = CreateDesign(h)

%% Perceptual learning
if isfield(h.Settings,'PL') && isfield(h.Settings,'AL')
    h = CreateSequence(h,{'PL','AL'});
elseif isfield(h.Settings,'PL')
    h = CreateSequence(h,{'PL'});
elseif isfield(h.Settings,'AL')
    h = CreateSequence(h,{'AL'});
end

%% Associative learning: add cue sequence
if isfield(h.Settings,'AL')
    %h = CreateSequence(h,'AL');
    h = CreateAssociative(h);
%elseif length(h.Settings.stim)>1
%    h.Seq.signal = repmat(h.Seq.signal,length(h.Settings.stim),1);
end

%% other
if ~isfield(h,'Seq')
    h = CreateSequence(h);
end

%% create Adaptive type order
if isfield(h.Settings,'adaptive')
    h = CreateAdaptive(h);
end

% create all trials if design is continuous
if isfield(h.Settings,'stimcontrol') && strcmp(h.Settings.design,'continuous') && h.Settings.savesinwave
    if ~isempty(h.Settings.stimcontrol)
        %if ~strcmp(h.Settings.stimcontrol,'labjack')
            opt = 'create';
            h = stimtrain(h,opt);
        %end
    end
end