function h = CreateAssociative(h)

if isfield(h.Seq,'PL')
    use_PL=1;
else
    use_PL=0;
end

% cue types
stimtypes = [1 2];

% initialise h.Seq
h.Seq.signal = nan(length(stimtypes),size(h.Seq.AL.signal,2));
h.Seq.blocks = h.Seq.AL.blocks;
h.Seq.condnum = h.Seq.AL.condnum;
h.Seq.cp_cond = h.Seq.AL.cp_cond;

% get block indices
blockuni=unique(h.Seq.AL.blocks);
for b = 1:length(blockuni)
    block_ind{b} = find(h.Seq.AL.blocks==blockuni(b));
end

for b = 1:length(blockuni)
    
    % assign values to second stim based on condnum
    conduni=unique(h.Seq.condnum(block_ind{b}));
    
    for c = 1:length(conduni)
        % for each block, get index of trials for each condnum
        cond_ind{b}{c} = find(h.Seq.blocks==blockuni(b) & h.Seq.condnum==conduni(c));
        
        if use_PL
            h.Seq.signal(2,cond_ind{b}{c}) = h.Seq.PL.signal(1,cond_ind{b}{c});
        else
            % assign random values for 2nd stim
            %\\\\make this conditional on not already havign been defined
            randval = repmat(stimtypes,1,ceil(length(cond_ind{b}{c})/length(stimtypes)));
            randval = randval(randperm(length(randval)));
            h.Seq.signal(2,cond_ind{b}{c}) = randval(1:length(cond_ind{b}{c}));
        end
        
        % then decide which pairing to use.
        % Depending on the probability during this block, one pairing may
        % be more likely than another
        pair_idx = h.Settings.AL.pair(conduni(c));
        
        % for the number of unique values in each cell of
        % h.Settings.AL.pairing (i.e. the number of stimulus outcomes
        % per pairing/cue), create a randomised sequence according to a
        % certain probability.
        % need to be balanced across different cue types.
        
        %initialise index
        cuestim_idx = nan(1,length(cond_ind{b}{c}));
        %for each stim type, get index of each type 
        stimtype = unique(h.Seq.signal(2,cond_ind{b}{c}));
        for st = 1:length(stimtype)
            stim_idx{st} = find(h.Seq.signal(2,cond_ind{b}{c})==stimtype(st));
            % get number of stims per stimtype and cuetype
            stim_num{st} = round(h.Settings.AL.probstim{conduni(c)} * length(stim_idx{st}));
            cue_idx{st} = [];
            for i = 1:length(stim_num{st})
                cue_idx{st} = [cue_idx{st} i*ones(1,stim_num{st}(i))];
            end
            cue_idx{st} = cue_idx{st}(randperm(length(cue_idx{st})));
            cuestim_idx(stim_idx{st}) = cue_idx{st};
        end
        
        % finally, decide actual cue according to which stim it is (i.e.
        % second row of signal)
        for t = 1:length(cond_ind{b}{c})
            h.Seq.signal(1,cond_ind{b}{c}(t)) = h.Settings.AL.pairing{pair_idx,h.Seq.signal(2,cond_ind{b}{c}(t))}(cuestim_idx(t));
        end
        
    end
    
    
    % diagnostics for debugging
    if 1
        for p = 1:4 % for each value of signal (in either row)
            for st = 1:2
                num_stimtype{b}(p,st) = length(find(h.Seq.signal(st,block_ind{b})==p));
            end
        end
    end
    
end

% update condnum so that there are unique values for high and low stims
conduni=unique(h.Seq.condnum);
sigrow=2;% stims
siguni=unique(h.Seq.signal(sigrow,:)); 
newcond=0;
newcondnum=nan(1,length(h.Seq.condnum));
newcondnumset = h.AL.condnum;
for c = 1:length(conduni)
    cond_ind = find(h.Seq.condnum==conduni(c));
    for s = 1:length(siguni)
        newcond=newcond+1;
        sig_ind = find(h.Seq.signal(sigrow,:)==siguni(s));
        comb_ind = ismember(cond_ind,sig_ind);
        newcondnum(cond_ind(comb_ind)) = newcond;
    end
end
h.Seq.condnum=newcondnum;
h.condnum = newcondnumset;

% update condnum if there is also PL
if use_PL
    [ur, i1, i2] = unique(h.Settings.PL.condnum,'rows','stable');
    cp_ind = nan(1,length(h.Seq.PL.cp_cond));
    for cp = 1:length(i1)
        cp_val = find(i2==cp);
        cp_ind(ismember(h.Seq.PL.cp_cond,cp_val)) = cp;
    end
    maxcon = max(h.Seq.condnum);
    h.Seq.condnum = h.Seq.condnum + (cp_ind-1)*maxcon;
    
    conduni=unique(h.Seq.condnum);
    newcondnum=h.Seq.condnum;
    for cu=1:length(conduni)
        newcondnum(h.Seq.condnum==conduni(cu))=cu;
    end
    h.Seq.condnum=newcondnum;
end

figure
plot(h.Seq.condnum)
hold on; 
low = ismember(h.Seq.signal(2,:),[1 3]);
condplot = nan(1,length(h.Seq.condnum));
condplot(low) = h.Seq.condnum(low);
s1 = scatter(1:length(h.Seq.condnum),condplot,'k');
hold on; 
high = ismember(h.Seq.signal(2,:),[2 4]);
condplot = nan(1,length(h.Seq.condnum));
condplot(high) = h.Seq.condnum(high);
s2 = scatter(1:length(h.Seq.condnum),condplot,'r');
title('design')
legend([s1,s2],'low intensity stim','high intensity stim')
    

