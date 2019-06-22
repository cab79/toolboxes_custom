function h = TriangularSequence(h)

dbstop if error

disp('Creating sequence...');
h.Seq =struct;
h.Seq.condnum = [];
seq=[];

if ~isfield(h.Settings,'nstim_trial')
    h.Settings.nstim_trial=3;
end

if ~isfield(h.Settings,'leading_stims')
    h.Settings.leading_stims=0;
end

if length(h.Settings.PL.oddballvalue{1})==1
    nodds = length(h.Settings.PL.oddballvalue);
else
    nodds = length(h.Settings.PL.oddballvalue{1});
end

% create list of which stim to present on each trial   
if nodds<=h.Settings.nstim_trial
% two possible stims (patterns)
        for i = 1:size(h.Settings.PL.oddprob,1)
            ntrial_per_stimtype(i,:) = floor(h.Settings.ntrials(i)*h.Settings.PL.oddprob(i,:));
            thisseq=[];
            for ns = 1:size(ntrial_per_stimtype,2)
                thisseq = [thisseq repmat(ns,1,ntrial_per_stimtype(i,ns))];
            end
            h.Seq.condnum=[h.Seq.condnum i*ones(1,length(thisseq))];
            seq = [seq thisseq];
        end
        randind = randperm(length(seq));
        seq = seq(randind); % oddballs or standards? In this case, which stim pattern to start with.
        h.Seq.condnum = h.Seq.condnum(randind);

        % use oddball values to create signal
        for i = 1:size(h.Settings.PL.oddprob,1)
            condind = h.Seq.condnum==i;
            h.Seq.signal(1:h.Settings.nstim_trial,condind)=h.Settings.PL.nonoddballstimvalue{i}*ones(h.Settings.nstim_trial,sum(condind));
            for ns = h.Settings.target_stims
                ind = find(h.Settings.target_stims==ns);
                h.Seq.signal(ns,condind) = h.Settings.PL.oddballvalue{i}(ind);
            end
        end
else % four possible stims (patterns)
        for i = 1:size(h.Settings.PL.oddprob,1)
            ntrial_per_stimtype(i,:) = floor(h.Settings.ntrials(i)*h.Settings.PL.oddprob(i,:));
            thisseq=[];
            for ns = 1:size(ntrial_per_stimtype,2)
                thisseq = [thisseq repmat(ns,1,ntrial_per_stimtype(i,ns))];
            end
            h.Seq.condnum=[h.Seq.condnum i*ones(1,length(thisseq))];
            seq = [seq thisseq];
        end
        randind = randperm(length(seq));
        seq = seq(randind);
        h.Seq.condnum = h.Seq.condnum(randind);
        h.Seq.signal(1:h.Settings.nstim_trial,length(h.Seq.condnum))=NaN;
        
        % for each stim pattern being left out
        for i = 1:size(h.Settings.PL.oddprob,1)
            
            % get all values except value i being left out
            oddind = 1:size(h.Settings.PL.oddprob,1);
            oddind(i) = [];
            oddvals = h.Settings.oddballvalue(oddind);
            oddvals = oddvals(randperm(length(oddvals)));
            
            % get the indices of this 'condition' (i.e. the pattern being left out)
            condind = h.Seq.condnum==i;
            for ns = h.Settings.target_stims
                h.Seq.signal(ns,condind) = oddvals{find(h.Settings.target_stims==ns)};
            end
        end
        
end
h.Seq.blocks=ones(1,length(h.Seq.condnum));

% create patterns by randomising sub-components of the stimulus train
if strcmp(h.Settings.stim(1).durtype,'sequence_rand')
    if ~isempty(h.Settings.stim(1).stimrandind)
        disp('Creating patterns...');
        ind1 = 1:length(h.Settings.stim(1).dur);
        ind1 = ind1(randperm(length(ind1)));
        ind2 = ind1;
        ind3 = 1:length(h.Settings.stim(1).dur);
        ind3 = ind3(randperm(length(ind3)));
        ind4 = ind3;
        x=1;
        while x
            randind1=randperm(length(h.Settings.stim(1).stimrandind));
            randind2=randperm(length(h.Settings.stim(1).stimrandind));
            randind3=randperm(length(h.Settings.stim(1).stimrandind));
            randind4=randperm(length(h.Settings.stim(1).stimrandind));
            if ~any(randind1==randind2) || ~any(randind3==randind4)
                x=0;
            end
        end
        ind1(h.Settings.stim(1).stimrandind) = h.Settings.stim(1).stimrandind(randind1);
        ind2(h.Settings.stim(1).stimrandind) = h.Settings.stim(1).stimrandind(randind2);
        ind3(h.Settings.stim(1).stimrandind) = h.Settings.stim(1).stimrandind(randind3);
        ind4(h.Settings.stim(1).stimrandind) = h.Settings.stim(1).stimrandind(randind4);

        patternvalue = h.Settings.stim(1).patternvalue; 
        h.Settings.stim(1).patternvalue = {};
        h.Settings.stim(1).patternvalue = {patternvalue(ind1),patternvalue(ind2),patternvalue(ind3),patternvalue(ind4)};
        for ns = 2:h.Settings.nstim_trial
            h.Settings.stim(ns).patternvalue = h.Settings.stim(1).patternvalue;
        end
    end
elseif strcmp(h.Settings.stim(1).durtype,'FromSequence')
    %load from sequence file
    global d
    file = load(fullfile(d.seq, h.Settings.stim(1).patternvalue));
    h.Settings.stim(1) = file.settings.stim(1);
    for ns = 2:h.Settings.nstim_trial
        h.Settings.stim(ns) = file.settings.stim(2);
        h.Settings.stim(ns) = file.settings.stim(3);
    end
end

%% create Adaptive type order
if isfield(h.Settings,'adaptive')
    h = CreateAdaptive(h);
end