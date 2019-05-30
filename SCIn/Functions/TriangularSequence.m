function h = TriangularSequence(h)

dbstop if error

disp('Creating sequence...');
h.Seq =struct;
h.Seq.condnum = [];
seq=[];

% create list of which stim to present on each trial   
switch length(h.Settings.PL.oddballvalue)
    case 2
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

        for i = 1:size(h.Settings.PL.oddprob,1)
            condind = h.Seq.condnum==i;
            h.Seq.signal(1:3,condind)=h.Settings.PL.oddballvalue{i}(1)*ones(3,sum(condind));
            h.Seq.signal(1,seq==1 & condind) = h.Settings.PL.oddballvalue{i}(2);
            h.Seq.signal(2,seq==2 & condind) = h.Settings.PL.oddballvalue{i}(2);
            h.Seq.signal(3,seq==3 & condind) = h.Settings.PL.oddballvalue{i}(2);
        end
    case 4
        for i = 1:size(h.Settings.oddprob,1)
            ntrial_per_stimtype(i,:) = floor(h.Settings.ntrials(i)*h.Settings.oddprob(i,:));
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
        h.Seq.signal(1:3,length(h.Seq.condnum))=NaN;
        
        % for each stim pattern being left out
        for i = 1:size(h.Settings.oddprob,1)
            
            % get all values except value i being left out
            oddind = 1:size(h.Settings.oddprob,1);
            oddind(i) = [];
            oddvals = h.Settings.oddballvalue(oddind);
            oddvals = oddvals(randperm(length(oddvals)));
            
            % get the indices of this 'condition' (i.e. the pattern being left out)
            condind = h.Seq.condnum==i;
            h.Seq.signal(1,condind) = oddvals{1};
            h.Seq.signal(2,condind) = oddvals{2};
            h.Seq.signal(3,condind) = oddvals{3};
        end
        
end
h.Seq.blocks=ones(1,length(h.Seq.condnum));

% create patterns
if strcmp(h.Settings.stim(1).durtype,'sequence_rand')
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
    h.Settings.stim(2).patternvalue = h.Settings.stim(1).patternvalue;
    h.Settings.stim(3).patternvalue = h.Settings.stim(1).patternvalue;
elseif strcmp(h.Settings.stim(1).durtype,'FromSequence')
    %load from sequence file
    global d
    file = load(fullfile(d.seq, h.Settings.stim(1).patternvalue));
    h.Settings.stim(1) = file.settings.stim(1);
    h.Settings.stim(2) = file.settings.stim(2);
    h.Settings.stim(3) = file.settings.stim(3);
end