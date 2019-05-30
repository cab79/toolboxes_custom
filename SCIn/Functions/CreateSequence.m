 function h = CreateSequence(h,varargin)
dbstop if error
disp('Creating sequence... wish me luck!');

setconds = h.Settings.conds;

if nargin>1
    % design type
    dtypes = varargin{1};
else 
    dtypes = {''};
end

for dt = 1:length(dtypes)
    if ~isempty(dtypes{1})
        % design type
        dtype = dtypes{dt};

        % the following fields are required in settings for this to run:
        setoddprob = h.Settings.(dtype).oddprob;
        setstandardind = h.Settings.(dtype).standardind; 
        setoddind = h.Settings.(dtype).oddind;
        setsep_odd = h.Settings.(dtype).sep_odd;
        setsep_odd_ind = h.Settings.(dtype).sep_odd_ind;
        setstd_lead = h.Settings.(dtype).std_lead;
        setn_set = h.Settings.(dtype).n_set;
        setn_odd = h.Settings.(dtype).n_odd;
        setn_odd_set = h.Settings.(dtype).n_odd_set;
        setrand_set = h.Settings.(dtype).rand_set;
        oddballmethod = h.Settings.(dtype).oddballmethod;
        oddballvalue = h.Settings.(dtype).oddballvalue;
        oddballtype = h.Settings.(dtype).oddballtype;
        try
            setcondnum = h.Settings.(dtype).condnum;
        catch
            setcondnum = [];
        end
        settol = h.Settings.(dtype).sep_odd_tol;
    else
        dtype='X';
        setoddprob = h.Settings.oddprob;
        setstandardind = h.Settings.standardind; 
        setoddind = h.Settings.oddind;
        setsep_odd = h.Settings.sep_odd;
        setsep_odd_ind = h.Settings.sep_odd_ind;
        setstd_lead = h.Settings.std_lead;
        setn_set = h.Settings.n_set;
        setn_odd = h.Settings.n_odd;
        setn_odd_set = h.Settings.n_odd_set;
        setrand_set = h.Settings.rand_set;
        oddballmethod = h.Settings.oddballmethod;
        oddballvalue = h.Settings.oddballvalue;
        oddballtype = h.Settings.oddballtype;
        try
            setcondnum = h.Settings.condnum;
        catch
            setcondnum = [];
        end
        settol = h.Settings.sep_odd_tol;
    end

    %% define orthgonal conditions (CURRENTLY HAS NO IMPACT ON REST OF THE FUNCTION)
    % e.g. ensuring that condition levels and oddball levels are orthogonal, if
    % needed
    % conds = [];
    % ncond = 0;
    % nonprobcond = [];
    % if ~isempty(setconds)
    %     allcondind='';
    %     for i = 1:length(setconds)
    %         eval(['condval{i} = set' setconds{i} ';']);
    %         if size(condval{i},1)<2
    %             continue
    %         end
    %         condind{i} = 1:size(condval{i},1);
    %         allcondind = [allcondind ' condind{' num2str(i) '},'];
    %     end
    %     eval(['allcond = allcomb(' allcondind(1:end-1) ');']);
    %     allcond = 1:i;
    %     probcond = find(strcmp(setconds,'oddprob'));
    %     nonprobcond = allcond; nonprobcond(probcond) = [];
    % end

    %% work out parameters of each sequence set (length, number of oddballs, probability of oddballs, etc.)

    %if ~isfield(h.Settings,'probdepend')
    %    h.Settings.probdepend = ones(size(setoddprob));
    %end

    % ensure probs are integers for gcd calculation
    probs = setoddprob*1000000;

    % number of CP conditions
    nCP = size(probs,1);

    % find greatest common divisor to find least number of repetitions of each
    % oddball
    gd=[];
    for i = 1:nCP
        gd(i,1) = probs(i,1);
        for ii=2:numel(probs(i,:))
            gd(i,1) = gcd(gd(i,1),probs(i,ii));
        end
    end
    if nCP==0
        mult=1;
    else
        mult = probs./repmat(gd,1,size(probs,2)); % multiplier is the min number of repetitions of each option (column) for each row (CP)
    end

    % adjust for minimum number of oddballs per set, with a min of 1
    minmult = ceil(setn_odd_set./max(min(mult'),1))';
    mult = mult.* repmat(minmult,1,size(mult,2));

    % increase set size
    %mult = repmat(setn_set',1,size(mult,2)).*mult;

    % calculate total duration of one set
    %tot = sum(mult(:)); % total number of the minimum set of stim types
    if strcmp(oddballmethod,'duration')
        stimdur = oddballvalue;
    else
        stimdur = [h.Settings.stim(:).dur];
    end
    if iscell(stimdur)
        dursum = cellfun(@sum,stimdur);
    else
        dursum = sum(stimdur,2);
    end
    totdur = sum(max(dursum,mean(h.Settings.trialdur)) .* mult(mult>0));% total duration of one set of all stim types

    % calculate number of sets that can provide at least h.Settings.totdur of stimulation AND n_odd oddballs per CP
    num_sets{dt}=0;
    if isfield(h.Settings,'totdur')
        num_sets{dt} = ceil(h.Settings.totdur/totdur);
    end
    if isempty(setn_set) && ~isempty(setn_odd)
        try
            num_sets{dt} = setn_odd./min(mult');
        catch
            num_sets{dt} = setn_odd./min(mult);
        end
    elseif ~isempty(setn_set)
        num_sets{dt} = setn_set;
    end

    h.Seq.totdur = sum(max(dursum,h.Settings.trialdur) * sum(num_sets{dt}*mult));% total duration of one set of all stim types



    %if ~isfield(h.Settings,'rand_within_set')
    %    h.Settings.rand_within_set = 1;
    %end
    
    % check
    if dt>1
        if ~isequal(num_sets{1},num_sets{2})
            error('num sets not equal')
        end
    end


    %% create sequence sets
    setx{dt} = []; cpx{dt}=[]; setnum=0;
    if nCP>0
        for cp = 1:nCP
            setind{cp} = [];
            for ii = 1:length(mult(cp,:))
                setind{cp} = [setind{cp} ii*ones(1,mult(cp,ii))];
            end

            % create a different randomised list (block) for each repeat of the set
            randind{cp} = [];
            for s = 1:num_sets{dt}(cp)

                % find sequence in which oddball trials are apart by at least nX standards
                nX = setsep_odd(cp);
                nXi = setsep_odd_ind;

                % remove first nR standards - not to be randomised, but added to the
                % start of each set later
                nR = setstd_lead(cp);
                setindnX = setind{cp}(nR+1:end);

                sequence_found = false;
                while ~sequence_found

                    candidate = setindnX(randperm(length(setindnX)));

                    no_conseq=1;min_stan=1;

                    for ii = 1:length(nXi)
                        % find indices of oddball types NOT to consider and
                        % make them 1 (as if standards)
                        sub_cand = candidate;
                        oi = nXi{ii};
                        oi_ind = ~ismember(candidate,oi);
                        sub_cand(oi_ind) = 1;


                        % How long is each sequence of standards?
                        w = [false sub_cand==setstandardind false]; %// "close" w with zeros, and transform to logical
                        starts = find(w(2:end) & ~w(1:end-1)); %// find starts of runs of non-zeros
                        ends = find(~w(2:end) & w(1:end-1))-1; %// find ends of runs of non-zeros
                        result = cell2mat(arrayfun(@(s,e) length(sub_cand(s:e)), starts, ends, 'uniformout', false)); %// build result
                        if sum(result>=nX)<length(result)*settol(1)
                            min_stan=0;
                        end

                        % must also be no consecutive oddballs
                        if nX>0
                            cand_odd = sub_cand>1;
                            diffcand = [diff(cand_odd) 0];
                            if sum(diffcand(cand_odd) ~= 0)<length(diffcand(cand_odd))*settol(2) %// check if no repeated values
                                no_conseq=0;
                            end
                        end
                    end

                    if min_stan && no_conseq 
                        sequence_found = true;
                    end
                end

                disp(['SETUP SEQUENCE: design type ' num2str(dt) ', CP ' num2str(cp) '/' num2str(nCP) ', Set ' num2str(s) '/' num2str(num_sets{dt}(cp)) ' complete']);

                randind{cp}{s} = [setind{cp}(1:nR) candidate];
                %if length(h.Settings.rand_within_set)>1
                %    if h.Settings.rand_within_set(cp)
                %end
                setnum = setnum+1;
                setx{dt} = [setx{dt} setnum*ones(1,length(candidate))];
                cpx{dt} = [cpx{dt} cp*ones(1,length(candidate))]; 
            end

            % for roving oddball, each oddball stimulus is a persistent change in the
            % stimulus characteristic (until the next oddball when it changes back).
            % Hence, the stimulus types need updating here to be consistent with this:

            for s = 1:num_sets{dt}(cp)
                if strcmp(oddballtype,'roving')
                    for i = 1:length(randind{cp}{s})
                        if i==1
                            h.(dtype).condnum{cp}{s}(i) = 0; % first stimulus in a run should be indicated with 0
                        elseif i==2 % initialise values
                            if randind{cp}{s}(i)==setstandardind % standard
                                h.(dtype).condnum{cp}{s}(i) = setstandardind;
                            elseif randind{cp}{s}(i)==setoddind % oddball
                                h.(dtype).condnum{cp}{s}(i) = setoddind;
                            end
                        elseif i>2 % values now depend on the direction of change
                            if randind{cp}{s}(i)==setstandardind % standard
                                if h.(dtype).condnum{cp}{s}(i-1)==1
                                    h.(dtype).condnum{cp}{s}(i) = 1;
                                elseif h.(dtype).condnum{cp}{s}(i-1)==2
                                    h.(dtype).condnum{cp}{s}(i) = 3;
                                elseif h.(dtype).condnum{cp}{s}(i-1)==3
                                    h.(dtype).condnum{cp}{s}(i) = 3;
                                elseif h.(dtype).condnum{cp}{s}(i-1)==4
                                    h.(dtype).condnum{cp}{s}(i) = 1;
                                end
                            elseif randind{cp}{s}(i)==setoddind % oddball
                                if h.(dtype).condnum{cp}{s}(i-1)==1
                                    h.(dtype).condnum{cp}{s}(i) = 2;
                                elseif h.(dtype).condnum{cp}{s}(i-1)==2
                                    h.(dtype).condnum{cp}{s}(i) = 4;
                                elseif h.(dtype).condnum{cp}{s}(i-1)==3
                                    h.(dtype).condnum{cp}{s}(i) = 4;
                                elseif h.(dtype).condnum{cp}{s}(i-1)==4
                                    h.(dtype).condnum{cp}{s}(i) = 2;
                                end
                            end
                        end
                    end

                    h.(dtype).stimtype{cp}{s} = nan(1,length(h.(dtype).condnum{cp}{s}));
                    h.(dtype).stimtype{cp}{s}(ismember(h.(dtype).condnum{cp}{s},[0])) = oddballvalue{cp}(1);
                    h.(dtype).stimtype{cp}{s}(ismember(h.(dtype).condnum{cp}{s},[1 4])) = oddballvalue{cp}(1);
                    h.(dtype).stimtype{cp}{s}(ismember(h.(dtype).condnum{cp}{s},[2 3])) = oddballvalue{cp}(2);
                else
                    stan_odd_val = [1 2];
                    if all(ismember(unique(randind{cp}{s}),stan_odd_val))
                        %update values according to oddball method
                        if strcmp(oddballmethod,'index')
                            if ~isempty(oddballvalue)
                                stan_odd_val = oddballvalue{cp};
                            end
                        end
                        h.(dtype).stimtype{cp}{s}(randind{cp}{s}==1)=stan_odd_val(1);
                        h.(dtype).stimtype{cp}{s}(randind{cp}{s}==2)=stan_odd_val(2);
                        h.(dtype).condnum{cp}{s}=h.(dtype).stimtype{cp}{s};
                    else
                        h.(dtype).stimtype{cp}{s}=randind{cp}{s};
                        h.(dtype).condnum{cp}{s}=randind{cp}{s};
                    end
                end
            end

            if strcmp(oddballtype,'roving')
                condnum_allset = cat(2,h.(dtype).condnum{cp}{:});
                % identify number of standards and oddballs
                cn=unique(condnum_allset);
                stan_ind = find(mod(cn,2)); % odd numbered
                odd_ind = find(~mod(cn,2)); % even numbered
                j=0;
                for i = cn
                    j=j+1;
                    ni(j) = length(find(condnum_allset==i));
                end
                % finds the numbers of two types of oddballs: those occurring
                % in change from stim1 to stim2, and those vice versa.
                h.Seq.(dtype).nodd{cp} = ni(odd_ind);
                h.Seq.(dtype).nstan{cp} = ni(stan_ind);

                % store indices
                h.Seq.(dtype).stan_ind{cp} = cn(stan_ind); % odd numbered
                h.Seq.(dtype).odd_ind{cp} = cn(odd_ind); % even numbered

                % find a new sequence if nstan are too different; otherwise complete
                %if max(h.Seq.nstan{cp})-min(h.Seq.nstan{cp}) > 1
                %    h = CreateSequence(h);
                %end
            end

            % make condition numbers distinct for different CP conditions
            for s = 1:num_sets{dt}(cp)
                if strcmp(oddballtype,'classical') && ~isempty(setcondnum)
                    ucon = unique(h.(dtype).condnum{cp}{s});
                    ucon = ucon(ucon>0);
                    for i = 1:length(ucon)
                        h.(dtype).condnum{cp}{s}(h.(dtype).condnum{cp}{s}==ucon(i)) = setcondnum(cp,ucon(i));
                    end
                elseif strcmp(oddballtype,'roving') && cp>1
                    for s = 1:num_sets{dt}(cp)
                        % find max value of previous CP condition
                        maxval = max([h.(dtype).condnum{cp-1}{:}]);
                        % find min value of current CP condition
                        vals=[h.(dtype).condnum{cp}{:}];
                        minval = min(vals(vals>0));
                        % is current value too low?
                        if minval<=maxval
                            % gat value to add on
                            valadd = maxval-minval+1;
                            % add this value to the h.(dtype).condnum
                            if strcmp(oddballtype,'classical')
                                h.(dtype).condnum{cp}{s} = h.(dtype).condnum{cp}{s}+valadd; 
                            elseif strcmp(oddballtype,'roving')
                                h.(dtype).condnum{cp}{s}(2:end) = h.(dtype).condnum{cp}{s}(2:end)+valadd; % keep first value as 0
                            end
                        end
                    end
                elseif isempty(setcondnum)
                    h.(dtype).condnum{cp}{s} = h.(dtype).stimtype{cp}{s};
                end
            end
        end
    end
end

% check
if dt>1
    if ~isequal(cpx{1},cpx{2})
        error('cp indices not equal')
    end
    if ~isequal(setx{1},setx{2})
        error('set indices not equal')
    end
end


%% duplicate sequences for non-probability conditions
% if ~isempty(nonprobcond)
%     
% end


%% create final sequences/blocks
%if ~isfield(h.Seq,'signal')
    
%randomise sets
if any(setrand_set)
    if length(setrand_set)>1
        % get set indices to randomise
        setrand = find(ismember(cpx{1},find(setrand_set)));
    else
        setrand = 1:length(setx{1});
    end
    % get set indices not beign randomised
    setnorand = 1:length(setx{1});
    setnorand(setrand)=[];

    % create set index
    us = unique(setx{1}(setrand));
    us_pos = ismember(unique(setx{1}),us);
    usnr = unique(setx{1}(setnorand));
    usnr_pos = ismember(unique(setx{1}),usnr);
    randus = us(randperm(length(us)));
    randnorand(us_pos) = randus;
    randnorand(usnr_pos) = usnr;
    setx_ind = [];
    for s = 1:length(randnorand)
        setx_ind = [setx_ind setx{1}(setx{1}==randnorand(s))];
    end
else
    setx_ind=setx{1};
end

%sets_all = unique(setx_ind,'stable');
    
for dt = 1:length(dtypes)
    if ~isempty(dtypes{1})
        dtype = dtypes{dt};
    else
        dtype='X';
    end
    
    h.Seq.(dtype).signal=nan(1,length(setx_ind));
    h.Seq.(dtype).condnum=nan(1,length(setx_ind));
    h.Seq.(dtype).blocks=nan(1,length(setx_ind));
    h.Seq.(dtype).cp_cond=nan(1,length(setx_ind));
    %if nCP>0
        cps=0;
        for cp = 1:nCP
            for s = 1:num_sets{dt}(cp)
                cps=cps+1;
                h.Seq.(dtype).signal(setx_ind==cps) = h.(dtype).stimtype{cp}{s}; % type of signal for each trial: intensity, pitch, duration or channel
                %h.Seq.(dtype).pattern = ; % type of temporal pattern of the signal within each trial
                h.Seq.(dtype).condnum(setx_ind==cps) = h.(dtype).condnum{cp}{s};
                h.Seq.(dtype).cp_cond(setx_ind==cps)=cp;
                
                % if we are working on a cp condition whose sets are
                % randomised, balance the conds between blocks:
                randcp = ismember(cp,find(setrand_set));
                nrandcp = length(find(setrand_set));
                if any(randcp)
                    if strcmp(h.Settings.blockopt,'cond')
                        h.Seq.(dtype).blocks(setx_ind==cps) = cp*ones(1,length(h.(dtype).stimtype{cp}{s}));
                    elseif strcmp(h.Settings.blockopt,'divide')
                        bv = s+(nCP-nrandcp); % block value
                        h.Seq.(dtype).blocks(setx_ind==cps) = bv*ones(1,length(h.(dtype).stimtype{cp}{s}));
                    else
                        h.Seq.(dtype).blocks(setx_ind==cps) = ones(1,length(h.(dtype).stimtype{cp}{s}));
                    end
                % otherwise assign block value according to CP value
                % not ideal - current only works if non-rand CP==1
                else
                    nrand = length(setrand_set);
                    if strcmp(h.Settings.blockopt,'cond')
                        h.Seq.(dtype).blocks(setx_ind==cps) = cp*ones(1,length(h.(dtype).stimtype{cp}{s}));
                    elseif strcmp(h.Settings.blockopt,'divide')
                        if cp==1
                            bv = 1; % block value
                            h.Seq.(dtype).blocks(setx_ind==cps) = bv*ones(1,length(h.(dtype).stimtype{cp}{s}));
                        elseif cp==nCP
                            bv = 1000; % block value
                            h.Seq.(dtype).blocks(setx_ind==cps) = bv*ones(1,length(h.(dtype).stimtype{cp}{s}));
                        end
                    else
                        h.Seq.(dtype).blocks(setx_ind==cps) = ones(1,length(h.(dtype).stimtype{cp}{s}));
                    end
                end
            end
        end
        % correct max block number
        blockuni=unique(h.Seq.(dtype).blocks);
        if length(blockuni)>1
            h.Seq.(dtype).blocks(h.Seq.(dtype).blocks==max(h.Seq.(dtype).blocks)) = blockuni(end-1)+1;
        end
        
        if strcmp(h.Settings.blockopt,'cond')
            % make block nums monotonic without changing their order
            blockuni=unique(h.Seq.(dtype).blocks);
            blockunistable=unique(h.Seq.(dtype).blocks,'stable');
            for b = 1:length(blockuni)
                block_ind{b} = find(h.Seq.(dtype).blocks==blockunistable(b));
            end
            for b = 1:length(blockuni)
                h.Seq.(dtype).blocks(block_ind{b}) = blockuni(b);
            end
        end
        
        % sort (if needed)
        [~,blockind] = sort(h.Seq.(dtype).blocks);
        h.Seq.(dtype).blocks = h.Seq.(dtype).blocks(blockind);
        h.Seq.(dtype).signal = h.Seq.(dtype).signal(blockind);
        h.Seq.(dtype).condnum = h.Seq.(dtype).condnum(blockind);
        h.Seq.(dtype).cp_cond=h.Seq.(dtype).cp_cond(blockind);
    %else
    %    h.Seq.(dtype).signal = ones(1,num_sets{dt}); % type of signal for each trial: intensity, pitch, duration or channel
    %    %h.Seq.(dtype).pattern = ; % type of temporal pattern of the signal within each trial
    %    h.Seq.(dtype).condnum = ones(1,num_sets{dt});
    %    %h.Seq.(dtype).changedist = design(3,:);
    %    if strcmp(h.Settings.blockopt,'cond')
    %        h.Seq.(dtype).blocks = ones(1,num_sets{dt});
    %    end
    %end
    
    figure
    plot(h.Seq.(dtype).condnum)
    figure
    plot(h.Seq.(dtype).signal,'g')
        
    
end

if strcmp(dtype,'X') || dt==1
    h.Seq.blocks=h.Seq.(dtype).blocks;
    h.Seq.signal=h.Seq.(dtype).signal;
    h.Seq.condnum=h.Seq.(dtype).condnum;
    h.Seq.cp_cond=h.Seq.(dtype).cp_cond;
end


% response probe sequence
if isfield(h.Settings,'RPconds')
    RPconds = h.Settings.RPconds;
    RPi = find(ismember(h.Seq.condnum,RPconds)); % indices of all possible stims to use
    nSelect = round(h.Settings.RPprob*length(RPi));
    RPrand = RPi(randperm(length(RPi)));
    RPselect = RPrand(1:nSelect);
    h.Seq.RP = zeros(1,length(h.Seq.condnum));
    h.Seq.RP(RPselect) = 1;
end