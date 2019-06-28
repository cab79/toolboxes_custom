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

if ~isfield(h.Settings,'randomise_conds_evenly')
    h.Settings.randomise_conds_evenly=1;
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
        if h.Settings.randomise_conds_evenly
            
            % create temporary cond/signal array
            condsig = [h.Seq.condnum;seq]';
            
            % Create index of all rows
            condsig_i = 1:length(condsig);
            
            % get the number of unique cond/stims; get the number of
            % repeats/rows of each combination
            [uniq,uA,uB] = unique(condsig,'rows');
            
            % Calculate the number needed from each unique combination to be randomised
            for n = 1:length(uniq)
                nuniq(n) = length(find(uB==n));
            end
            
            % divide set of numbers by the minimum number to generate an error if not divisible.
            minuniq = nuniq/min(nuniq);
            if ~all(round(minuniq)==minuniq)
                error('numbers of each unique cond/stim should be divisible')
            end
            
            % Master sequence
            randind=[];
            
            % row indices in cell array, so that indices can be removed
            % with each sub-sequence randomisation
            for n = 1:length(uniq)
                %get index of all matching rows
                rowind{n} = find(uB==n);
            end
            
            %For the number of sub-sequences to be randomised
            for mn = 1:min(nuniq)
                % sub sequence
                sub_seq = [];
                %For each unique combination
                for n = 1:length(uniq)
                    %get index of all matching rows
                    rowrand = rowind{n};
                    %take first n row indices � add to a sequence
                    sub_seq = [sub_seq rowrand(1:minuniq(n))];
                    % remove indices already includes
                    rowind{n}(1:minuniq(n))=[];
                end
                %randomise and add to master sequence
                randind = [randind sub_seq(randperm(length(sub_seq)))];
            end

        else
            randind = randperm(length(seq));
        end
        seq = seq(randind); % oddballs or standards? In this case, which stim pattern to start with.
        h.Seq.condnum = h.Seq.condnum(randind);

        % use oddball values to create signal
        for i = 1:size(h.Settings.PL.oddprob,1)
            condind = h.Seq.condnum==i;
            n_nonodd = h.Settings.nstim_trial-length(h.Settings.target_stims);
            if length(h.Settings.PL.nonoddballstimvalue{i})<n_nonodd
                h.Settings.PL.nonoddballstimvalue{i} = repmat(h.Settings.PL.nonoddballstimvalue{i},1,n_nonodd);
            end
            if n_nonodd
                h.Seq.signal(1:n_nonodd,condind)=repmat(h.Settings.PL.nonoddballstimvalue{i}',1,sum(condind));
            end
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
elseif strcmp(h.Settings.stim(1).durtype,'oddballvalue')
    for s = 1:length(h.Settings.stim)
        patternvalue = h.Settings.stim(s).patternvalue;
        dur = h.Settings.stim(s).dur;
        for i = 1:length(patternvalue)
            h.Settings.stim(s).patternvalue{i*2-1} = patternvalue{i};
            h.Settings.stim(s).patternvalue{i*2} = patternvalue{i};
            h.Settings.stim(s).dur{i*2-1} = dur{i};
            h.Settings.stim(s).dur{i*2} = dur{i};
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