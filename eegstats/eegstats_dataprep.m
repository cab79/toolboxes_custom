function D = eegstats_dataprep(S)
% inputs: S is a settings structure (see example script)
% outputs: D is the prepared EEG data and predictors as a struct

%% preliminaries
D = struct;
if strcmp(S.prep.parallel.option,'condor')  
    delete('input*.mat');
    delete('output*.mat');
    delete('*.out');
    delete('*.log');
    delete('*.err');
end
% change to the input directory
eval(sprintf('%s', ['cd(''' S.prep.path.inputs.eeg ''')']));

%% find data using getfilelist.m function
GFL=S;
GFL.func = 'prep';
GFL.path.file = S.prep.path.inputs.eeg;
GFL.path.datfile = S.prep.path.inputs.subjects;
if isfield(S.prep.select,'prefixes') && ~isempty(S.prep.select.prefixes); GFL.(GFL.func).load.prefix = S.prep.select.prefixes; end
if isfield(S.prep.select,'suffixes') && ~isempty(S.prep.select.suffixes); GFL.(GFL.func).load.suffix = S.prep.select.suffixes; end
GFL = getfilelist(GFL);
if isfield(GFL.(GFL.func),'designtab')
    designmat=GFL.(GFL.func).designtab;
else
    designmat=cell2table(GFL.(GFL.func).designmat(2:end,:));
    designmat.Properties.VariableNames = GFL.(GFL.func).designmat(1,:);
end
filelist = designmat.file;
subjects = GFL.(GFL.func).select.subjects;
if isempty(filelist)
    error('No files found to import!\n');
end

% run though all EEG files in a loop
for d = 1:length(subjects)
    
    % create output D
    D(d).prep.subname = subjects(d);
        
    % FIND THE FILES FOR THIS SUBJECT
    subfiles = filelist(find(not(cellfun('isempty', strfind(filelist,D(d).prep.subname{:})))));
    
    if isempty(subfiles)
        continue
    end
    
    disp(['subject ' num2str(d) '/' num2str(length(subjects))])  

    % CYCLE THROUGH EACH FILE FOR THIS SUBJECT
    % add data to D structure which has a common format for any data type
    for f = 1:length(subfiles)
        filename = subfiles{f};
        loadsuffix = S.prep.select.suffixes{1};
        fprintf('\nImporting %s.\n\n', filename);

        switch S.prep.fname.ext{:}
            case 'mat'
                temp = load(filename);
            case 'set'
                temp = pop_loadset(filename);
        end
        fename = fieldnames(temp(f));
        
        for fn = 1:length(fename)
            D(d).prep.(fename{fn})(f).dat = temp.(fename{fn});
        end

        % get other filetypes with same name
        for fn = 1:length(S.prep.select.suffixes)-1
            loadsuffix = S.prep.select.suffixes{fn+1};
            filename = strrep(filename,S.prep.select.suffixes{1},loadsuffix);

            temp = load(filename);
            fename = fieldnames(temp(f));
            for fn = 1:length(fename)
                if isstruct(temp.(fename{fn}))
                    D(d).prep.(fename{fn})(f) = temp.(fename{fn});
                else
                    D(d).prep.(fename{fn})(f).dat = temp.(fename{fn});
                end
            end
        end
    end

    % get data information
    switch S.prep.fname.ext{:}
        case 'mat'
            if isfield(S.prep.select,'freq')
                timecourse = squeeze(D(d).prep.fdata.dat{1, 1}.powspctrm(:,:,S.prep.select.freq,:));
                n_chans = length(D(d).prep.fdata.dat{1, 1}.label);
                n_samples = length(D(d).prep.fdata.dat{1, 1}.time{1});
                n_trials = length(D(d).prep.fdata.dat{1, 1}.time);
                timecourse = reshape(permute(timecourse,[2 3 1]),n_chans,[]);
                eventType = D(d).prep.fdata.dat{1, 1}.conds; tnums = D(d).prep.fdata.dat{1, 1}.tnums; fnums = D(d).prep.fdata.dat{1, 1}.fnums; bnums = D(d).prep.fdata.dat{1, 1}.bnums;
            elseif isfield(S.prep.fname,'struct')
                eeg = D(d).prep.(fename{1})(1).dat.(S.prep.fname.struct{2}); 
                n_chans = size(eeg,1);
                n_samples = size(eeg,2);
                n_trials = size(eeg,3);
                eventType = ones(1,n_trials);
                tnums=1:length(eventType); % assumes all trials present
                D(d).prep.dim = [n_chans,n_samples,n_trials];
                timecourse = reshape(eeg,n_chans,[]); % make 2D: chans x (time*trials)
                % save other info
                other = fieldnames(D(d).prep.(fename{1})(1).dat);
                other = setdiff(other,S.prep.fname.struct{2});
                for fo = 1:length(other)
                    D(d).prep.other.(other{:}) = D(d).prep.(fename{1})(1).dat.(other{:});
                end
            else % this is for processing group ICA data files
                try
                    topography = D(d).prep.topography.dat;
                end
                timecourse = D(d).prep.timecourse.dat;
                eventType = D(d).prep.eventType.dat;
                n_chans = D(d).prep.n_chans.dat;
                n_samples = D(d).prep.n_samples.dat;
                n_trials = D(d).prep.n_trials.dat;
                tnums=1:length(eventType); % assumes all trials present
            end
        case 'set'
            % For data recorded from an EGI system with STIM markers for
            % timing events
            if any(find(strcmp('STIM',temp.epoch(1).eventtype)))
                [eventType,tnums, fnums, bnums] = get_markers(temp,'EGI');
            else % from other systems, e.g. Brainamps
                [eventType,tnums, fnums, bnums] = get_markers(temp,'BV');
            end
            
            % This is a study-specific patch for CORE - can remove
            % eventually.
            % Make correction for CORE006, part2 data
                % block 2 tnums incorrectly restart from 0, when should start
                % from 452
            if strcmp(D(d).prep.subname,'CORE006') && length(tnums)<1000 % i.e. it is part2 data
                tnums(bnums==2)=tnums(bnums==2)+451;
            end
            
            n_chans = D(d).prep.nbchan.dat;
            n_samples = D(d).prep.pnts.dat;
            n_trials = D(d).prep.trials.dat;
            D(d).prep.dim = [n_chans,n_samples,n_trials];
            timecourse = reshape(D(d).prep.data.dat,n_chans,[]); % make 2D: chans x (time*trials)
    end
    D(d).prep.tnums = tnums;
    
    if exist('topography','var')
        comps = 1:size(timecourse,1);
        data= topography(:,comps) * timecourse(comps,:);  
    else
        data= timecourse;  
    end
    if isfield(D(d).prep,'times')
        total_samples = D(d).prep.times.dat;
    else
        total_samples=S.prep.original_samples;
    end
    if ~isempty(S.prep.select.samples)
        select_samples=S.prep.select.samples;
    else
        select_samples=total_samples;
    end
    
    
    %% Predictors: factors from EEG markers - ASSUMES THEY ARE NUMERIC IN EEG DATA FILE but this can be updated if needed
    % create predictor variables
    dtab=table; % initialise data table for this subject
    
    if ~isempty(S.prep.pred.factor.markers) 
        factor_markers=S.prep.pred.factor.markers;
        
        if ~iscell(factor_markers{1})
            factor_markers = {factor_markers}; % for backwards compatability with some old scripts
        end
        for fac = 1:length(factor_markers)
        
            % for each level of the factor
            for fac_lev = 1:length(factor_markers{fac})
                condidx{fac_lev}=[];
                % for each marker, get indices over EEG trials
                for marker = 1:length(factor_markers{fac}{fac_lev})
                    condidx{fac_lev} = [condidx{fac_lev} find(ismember(eventType,factor_markers{fac}{fac_lev}(marker)))];
                end
                % pool into conditions defined by factor_markers, indexing
                % from zero
                predtemp(condidx{fac_lev},1)=fac_lev-1;
            end
            if any(S.prep.pred.factor.ordinal==fac)
                pred = categorical(predtemp,unique(predtemp),'ordinal',1);
            else
                pred = categorical(predtemp,unique(predtemp));
            end
            if isfield(S.prep.pred.factor,'label') && ~isempty(S.prep.pred.factor.label)
                pred_label=S.prep.pred.factor.label{fac};
            else
                pred_label=['fac' num2str(fac)];
            end
            dtab.(pred_label) = pred;
            clear predtemp predlabel
        end
    end
    
    %% select factor levels
    orig_tnums = tnums;
    if ~isempty(S.prep.pred.factor.select)
        for f = 1:size(S.prep.pred.factor.select,1)
            pred_label = S.prep.pred.factor.select{f,1};
            disp(['selecting conditions from factor ' pred_label])
            fac_levels = categorical(S.prep.pred.factor.select{f,2});
            dat = categorical(dtab.(pred_label));
            idx = ismember(dat,fac_levels);
            dtab = dtab(idx,:);
            tnums = tnums(idx);
            fnums = fnums(idx);
            bnums = bnums(idx);
            eventType = eventType(idx);
            if length(unique(dtab.(pred_label)))==1
                dtab.(pred_label) =[];
            elseif length(unique(dtab.(pred_label)))==2
                dtab.(pred_label) = categorical(dtab.(pred_label)); % remove ordinality
                dtab.(pred_label) = removecats(dtab.(pred_label));
            end
        end
    end
    D(d).prep.tnums = tnums;
    
    
    %% EEG data operations / transformations / calculations
    % this is done after defining predictors in case subtractions are
    % required

    % smooth
    if S.prep.calc.eeg.smooth_samples
        disp('smoothing...')
        for i = 1:size(data,1)
            data(i,:) = smooth(data(i,:),S.prep.calc.eeg.smooth_samples,'moving');
        end
    end

    % downsample
    if S.prep.calc.eeg.dsample
        orig_samples = total_samples;
        total_samples = downsample(total_samples',S.prep.calc.eeg.dsample)';
        select_samples = downsample(select_samples',S.prep.calc.eeg.dsample)';
        idx = repmat(ismember(orig_samples,total_samples),1,n_trials);
        data = data(:,idx);
    end
    
    % make 3D
    data=reshape(data,size(data,1),[],n_trials);
    
    % select data
    data = data(:,:,ismember(orig_tnums,tnums));
    n_trials = size(data,3);

    % flip channels right to left
    if ~isempty(S.prep.calc.eeg.flipchan)
        load(S.prep.path.inputs.chanlocs);
        flipidx = ismember(eventType,S.prep.calc.eeg.flipchan);
        data(:,:,flipidx) = flipchan(data(:,:,flipidx),chanlocs,S.prep.path.inputs.GSNlocs);
    end
    
    % subtractions
%     if ~isempty(S.prep.pred.factor.subtract)
%         ...
%     end
    
    % Data transformation
    if ~isempty(S.prep.calc.eeg.transform)
        x=data;
        if strcmp(S.prep.calc.eeg.transform,'arcsinh')
            data=log(x+sqrt(x.^2+1));
        end
    end
    
    % baseline correct
    if ~isempty(S.prep.calc.eeg.baseline_correct)
        bp = dsearchn(total_samples',S.prep.calc.eeg.baseline_correct');
        data = bsxfun(@minus, data, mean(data(:,bp(1):bp(2),:),2));
    end

    % select samples
    if exist('select_samples','var') && ~isempty(select_samples)
        select = dsearchn(total_samples',[select_samples(1),select_samples(end)]');
        data = data(:,select(1):select(2),:);
    end
    
    % trim
%     S.prep.calc.eeg.ndec=8; % trim data to a number of decimal places
        
    %% reshape data into samples and combine over subjects
    D(d).prep.dim = size(data);
    D(d).prep.data=reshape(data,[],n_trials)';
    clear data % to free memory

    %% separate EEG/predictor data into training and testing fractions, only for testing decoders
    trainfrac=S.prep.traintest.frac;
    switch S.prep.traintest.type
        case 'random'
            % Split into training (encoding) and testing (decoding) fractions
            % Produces S.trainidx and S.testidx. These are indices
            % of conData-ordered trials that will be used for training and
            % testing.
            trainidx = {[]};
            testidx = {[]};
            if S.prep.traintest.balance_conds
                ucond = unique(eventType);
                for u = 1:length(ucond)
                    cond_idx = find(eventType==ucond(u)); % for each condition, it's index within conData
                    % DO THIS FOR A NUMBER OF RUNS - MULTIPLE RANDOMISED INDICES
                    for run = 1:S.prep.traintest.num_runs
                        rng(run); % seed random number generator for consistency
                        trainidx{run} = [trainidx{run} randsample(cond_idx,round(trainfrac*length(cond_idx)))]; 
                    end
                end
            else
                % DO THIS FOR A NUMBER OF RUNS - MULTIPLE RANDOMISED INDICES
                for run = 1:S.prep.traintest.num_runs
                    rng(run); % seed random number generator for consistency
                    trainidx{run} = randsample(length(eventType),round(trainfrac*length(eventType)));
                end
            end
            % create test index
            train_matrix = zeros(length(tnums),length(trainidx));
            test_matrix=zeros(length(tnums),length(testidx));
            for run = 1:length(trainidx)
                trainidx{run} = sort(trainidx{run});
                if trainfrac<1
                    testidx{run} = 1:length(eventType);
                    testidx{run}(trainidx{con}) = [];
                else
                    % otherwise, duplicate - same trials for each
                    testidx{run} = trainidx{run};
                end
                % convert to design matrix
                train_matrix(trainidx{run},run) = 1;
                test_matrix(testidx{run},run) = 1;
            end

    end
    
    % store info in table format
    D(d).prep.dtab=dtab;
    D(d).prep.dtab.ID=repmat(D(d).prep.subname,length(tnums),1);
    D(d).prep.dtab.group=repmat(designmat.groups(d),length(tnums),1);
    D(d).prep.dtab.eventTypes=categorical(eventType');
    if exist('fnums','var'); D(d).prep.dtab.fnums=categorical(fnums'); end
    if exist('bnums','var'); D(d).prep.dtab.bnums=categorical(bnums'); end
    if exist('tnums','var'); D(d).prep.dtab.tnums=tnums'; end
    D(d).prep.dtab.train=categorical(train_matrix);
    D(d).prep.dtab.test=categorical(test_matrix);
    %D(d).prep.tnums = tnums;
    
end

% clean up D
E=struct;
for d = 1:length(D)
    E(d).prep.subname = D(d).prep.subname;
    E(d).prep.data = D(d).prep.data;
    E(d).prep.dtab = D(d).prep.dtab;
    E(d).prep.dim = D(d).prep.dim;
    if isfield(D(1).prep,'other')
        E(d).prep.other = D(d).prep.other;
    end
    %E(d).prep.tnums = D(d).prep.tnums;
end
D=E;
clear E

if S.prep.output.save
    disp('saving data to disk...')
    save(fullfile(S.prep.path.outputs,[S.prep.sname '.mat']),'D','S','-v7.3')
    disp('...done')
end


function [conds, tnums, fnums, bnums] = get_markers(EEG,type)

conds = nan(1,length(EEG.epoch));
tnums = nan(1,length(EEG.epoch));
fnums = nan(1,length(EEG.epoch));
bnums = nan(1,length(EEG.epoch));

switch type
    case 'EGI'
        for ep = 1:length(EEG.epoch)
            stimevidx = find(strcmp('STIM',EEG.epoch(ep).eventtype));
            if ep<length(EEG.epoch); stimevidx1 = find(strcmp('STIM',EEG.epoch(ep+1).eventtype));end
            if ~isempty(stimevidx)
                stimcodes = EEG.epoch(ep).eventcodes{stimevidx(end)};
                if ~any(strcmp('CNUM',stimcodes(:,1)))
                    error('change CNUM to FNUM to analyse conditions');
                end
                conds(1,ep) = stimcodes{strcmp('CNUM',stimcodes(:,1)),2};
                tnums(1,ep) = stimcodes{strcmp('TNUM',stimcodes(:,1)),2};
                fnums(1,ep) = stimcodes{strcmp('FNUM',stimcodes(:,1)),2};
                bnums(1,ep) = stimcodes{strcmp('BNUM',stimcodes(:,1)),2};

            else
                if length(stimevidx1)==2
                    stimcodes = EEG.epoch(ep+1).eventcodes{stimevidx1(1)};
                    conds(1,ep) = stimcodes{strcmp('CNUM',stimcodes(:,1)),2};
                    tnums(1,ep) = stimcodes{strcmp('TNUM',stimcodes(:,1)),2};
                    fnums(1,ep) = stimcodes{strcmp('FNUM',stimcodes(:,1)),2};
                    bnums(1,ep) = stimcodes{strcmp('BNUM',stimcodes(:,1)),2};
                else
                    error(['too many / too few STIMs on trial ' num2str(ep+1)])
                end
            end
        end
    case 'BV'
        evtypes = unique(horzcat(EEG.epoch(:).eventtype));
        for ep = 1:length(EEG.epoch)
            stimevidx = find([EEG.epoch(ep).eventlatency{:}]==0);
            if ~isempty(stimevidx)
                str=EEG.epoch(ep).eventtype{stimevidx};
                conds(1,ep) = str2double(regexp(str,'.* (\d+)','tokens','once'));
                tnums(1,ep) = EEG.epoch(ep).eventurevent{stimevidx};
            end
        end
end





