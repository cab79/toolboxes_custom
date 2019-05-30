function h = CreateAdaptive(h)

%% create Adaptive type order

if length(h.Settings.adaptive)>1 && isfield(h.Settings,'adaptive_general')

    % Create order of runs
    if strcmp(h.Settings.adaptive_general.seqtype,'alt')
        order = repmat(h.Settings.adaptive_general.adapttypes,1,h.Settings.adaptive(1).nRuns);
    elseif strcmp(h.Settings.adaptive_general.seqtype,'rand')
        order=[];
        nRuns=[h.Settings.adaptive.nRuns];
        bs=h.Settings.adaptive_general.seqrandblocksize;
        num_blocks = sum(nRuns)/bs;
        for nb=1:num_blocks
            bi = bs*(nb-1)+1:bs*nb;
            runs = round(([nRuns(1)/sum(nRuns),nRuns(2)/sum(nRuns)] * bs));
            order_bi=[];
            for i = 1:length(h.Settings.adaptive_general.adapttypes)
                order_bi = [order_bi h.Settings.adaptive_general.adapttypes(i)*ones(1,runs(i))];
            end
            order(bi) = order_bi(randperm(length(order_bi)));
        end
    elseif strcmp(h.Settings.adaptive_general.seqtype,'cond')
        order=[];
        for cp = 1:nCP
            idx = ismember(h.Seq.condnum,h.condnum{cp}{:});
            order(idx)=h.Settings.adaptive_general.seqtypecond(cp);
        end
    end

    % Create trial order
    uorder = unique(order)';
    for i = 1:length(uorder)
        ntype(i) = h.Settings.adaptive(i).trialsperrun;
    end

    % create sequence
    aseq=[];
    for i = 1:length(order)
        aseq = [aseq order(i)*ones(1,ntype(order(i)))];
    end

    % get condition numbers and their indices
    condind=[];
    if isfield(h.Settings.adaptive_general.selectcond,'cp') % use cp condition
        condind = ismember(h.Seq.cp_cond,h.Settings.adaptive_general.selectcond.cp);
%         for cpi = 1:length(h.Settings.adaptive_general.selectcond.cp)
%             condval = [condval unique(h.condnum{h.Settings.adaptive_general.selectcond.cp(cpi)}{1})];
%         end
    end
    if isempty(condind) || ~any(condind)
%         condind = find(ismember(h.Seq.condnum,condval));
%     else
        condind = 1:length(h.Seq.condnum);
    end

    % check Seq.signal is long enough
    if length(h.Seq.signal)<length(condind)
        error('sequence is not long enough for this number of adaptive trials')
    end
    if length(aseq)<length(condind)
        error('number of adaptive trials is not long enough for this sequence')
    end

    % add to Seq
    h.Seq.adapttype = nan(1,length(h.Seq.condnum));
    h.Seq.adapttype(condind) = aseq(1:length(condind));
else
    h.Seq.adapttype = ones(1,length(h.Seq.condnum));
end
hold on; plot(h.Seq.adapttype,'r')


