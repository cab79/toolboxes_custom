function [S,D] = eegstats_dataprep(S)
% inputs: S is a settings structure (see example script)
% outputs: D is the prepared EEG data and predictors as a long-format table
% (all subjects) or series of tables (one per subject)

%% preliminaries
D = struct;
if strcmp(S.parallel.option,'condor')  
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
GFL.load.suffix = S.prep.select.suffixes(1);
GFL.func = 'prep';
GFL.path.file = S.prep.path.inputs.eeg;
GFL.path.datfile = S.prep.path.inputs.subjects;
GFL = getfilelist(GFL);
designmat=cell2table(GFL.designmat(2:end,:));
designmat.Properties.VariableNames = GFL.designmat(1,:);
filelist = GFL.filelist;
subjects = GFL.select.subjects;
if isempty(filelist)
    error('No files found to import!\n');
end

% run though all EEG files in a loop
for d = 1:length(subjects)
    
    % create output D
    D.prep(d).subname = subjects{d};
        
    % FIND THE FILES FOR THIS SUBJECT
    subfiles = filelist(find(not(cellfun('isempty', strfind(filelist,D.prep(d).subname)))));
    
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

        switch S.prep.fname.eeg.ext{:}
            case 'mat'
                temp = load(filename);
            case 'set'
                temp = pop_loadset(filename);
        end
        fename = fieldnames(temp(f));
        
        for fn = 1:length(fename)
            D.prep(d).(fename{fn})(f).dat = temp.(fename{fn});
        end

        % get other filetypes with same name
        for fn = 1:length(S.prep.select.suffixes)-1
            loadsuffix = S.prep.select.suffixes{fn+1};
            filename = strrep(filename,S.prep.select.suffixes{1},loadsuffix);

            temp = load(filename);
            fename = fieldnames(temp(f));
            for fn = 1:length(fename)
                if isstruct(temp.(fename{fn}))
                    D.prep(d).(fename{fn})(f) = temp.(fename{fn});
                else
                    D.prep(d).(fename{fn})(f).dat = temp.(fename{fn});
                end
            end
        end
    end

    % get data information
    switch S.prep.fname.eeg.ext{:}
        case 'mat'
            if isfield(S.prep.select,'freq')
                timecourse = squeeze(D.prep(d).fdata.dat{1, 1}.powspctrm(:,:,S.prep.select.freq,:));
                n_chans = length(D.prep(d).fdata.dat{1, 1}.label);
                n_samples = length(D.prep(d).fdata.dat{1, 1}.time{1});
                n_trials = length(D.prep(d).fdata.dat{1, 1}.time);
                timecourse = reshape(permute(timecourse,[2 3 1]),n_chans,[]);
                eventType = D.prep(d).fdata.dat{1, 1}.conds; tnums = D.prep(d).fdata.dat{1, 1}.tnums; fnums = D.prep(d).fdata.dat{1, 1}.fnums; bnums = D.prep(d).fdata.dat{1, 1}.bnums;
            else % this is for processing group ICA data files
                try
                    topography = D.prep(d).topography.dat;
                end
                timecourse = D.prep(d).timecourse.dat;
                eventType = D.prep(d).eventType.dat;
                n_chans = D.prep(d).n_chans.dat;
                n_samples = D.prep(d).n_samples.dat;
                n_trials = D.prep(d).n_trials.dat;
                tnums=1:length(eventType); % assumes all trials present
            end
        case 'set'
            % For data recorded from an EGI system with STIM markers for
            % timing events
            if any(find(strcmp('STIM',temp.epoch(1).eventtype)))
                [eventType,tnums, ~, bnums] = get_markers(temp);
            else % from other systems, e.g. Brainamps
                eventType = S.prep.select.markers{:};
                tnums=1:length(eventType); % assumes all trials present
            end
            
            % This is a study-specific patch for CORE - can remove
            % eventually.
            % Make correction for CORE006, part2 data
                % block 2 tnums incorrectly restart from 0, when should start
                % from 452
            if strcmp(D.prep(d).subname,'CORE006') && length(tnums)<1000 % i.e. it is part2 data
                tnums(bnums==2)=tnums(bnums==2)+451;
            end
            
            n_chans = D.prep(d).nbchan.dat;
            n_samples = D.prep(d).pnts.dat;
            n_trials = D.prep(d).trials.dat;
            timecourse = reshape(D.prep(d).data.dat,n_chans,[]); % make 2D: chans x (time*trials)
    end
    
    if exist(topography,'var')
        comps = 1:size(timecourse,1);
        data= topography(:,comps) * timecourse(comps,:);  
    else
        data= timecourse;  
    end
    total_samples=S.prep.original_samples;
    
    
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
                pred_label=[pred_label S.prep.pred.factor.label{fac}];
            else
                pred_label=[pred_label ['fac' num2str(fac)]];
            end
            dtab.(pred_label) = pred;
            clear predtemp predlabel
        end
    end
    
    %% Covariates
    S.prep.path.covariate(1).type = 'HGF'; % type of data input file
    S.prep.path.covariate(1).input = 'Q:\Projects\CORE\behaviour\hgf\fitted\D_fit_r1_it9_n30r.mat'; % path and file
    S.prep.path.covariate(1).level = 'trial'; % subject, trial
    S.prep.path.covariate(1).supp = 'Q:\Projects\CORE\behaviour\hgf\fitted\CORE_fittedparameters_percmodel2_respmodel4_fractrain0_20190211T074650.mat'; % supplementary file if needed
    S.prep.path.covariate(1).def = { % definition of covariates within the file
        {'PL','epsi',[0],[]}; % precision-weighted prediction errors
    %     {'RT'}
    %     {'choice'}
        };
    CVs=S.prep.path.covariate;
    for ci = 1:length(CVs)
        switch covariate(ci).type
            case 'HGF'
                pred_out = HGF_covariates(...
                    covariate(1),...
                    D.prep(d).subname,...
                    tnums,...
                    S.prep.path.inputs.subjects);
        end
        % add suffixes if labels are replicated
        suffixes={'a','b','c','d','e'};
        % add pred_out to dtab
        ...
    end
    
    
     if any(strcmp(S.pred_type,'RT'))
        HGF=load(fullfile(S.path.hgf,S.file.hgf{1})); D_fit=HGF.D_fit; % prevents S loading.
        ind = find(strcmp({D_fit(:).subname},D.prep(d).subname));
        D.prep(d).HGF.dat = D_fit(ind).HGF;
        predtemp = D.prep(d).HGF.dat.y(:,2);
        
        % remove trials not in EEG and put in EEG trial order
        if ~any(diff(tnums)>1)
            error('tnums is probably wrong')
        end
        predtemp=predtemp(tnums,:);
        
        pred = [pred predtemp];
        pred_label=[pred_label 'RT'];
        clear predtemp  predlabel
        
        
    end
    if any(strcmp(S.pred_type,'HGF'))
        % load and format HGF data
        % for each traj group... reformat into predictor matrix
        suffixes={'a','b','c','d','e'};
        for nC = S.use_hgf
            HGF=load(fullfile(S.path.hgf,S.file.hgf{nC})); D_fit=HGF.D_fit; % prevents S loading.
            if ~isfield(D_fit,'dt')
                dtfile = load(fullfile(S.path.hgf,S.file.get_dt));

                pdata = readtable(S.path.datfile);
                subs = pdata.Subject(find(pdata.Include));
                subind = ismember({dtfile.D_fit(:).subname},subs);
    %             D_fit = D_fit(subind);
    %             D_fit.subname=dtfile.D_fit(subind).subname;
    %             D_sub.dt=dtfile.D_fit(subind).dt;
                [D_fit(:).dt]=deal(dtfile.D_fit(subind).dt);
                [D_fit(:).subname]=deal(dtfile.D_fit(subind).subname);
            end
            ind = find(strcmp({D_fit.subname},D.prep(d).subname));
            if isempty(ind)
                disp('USING ONE HGF FOR ALL SUBJECTS')
                ind=1;
            end
            D.prep(d).HGF = D_fit.HGF;
            if length(S.use_hgf)>1
                [Stemp] = HGF_traj2mat(S,D,suffixes{nC});
            else
                [Stemp] = HGF_traj2mat(S,D);
            end
            predtemp=Stemp.pred;
            predgroup=Stemp.pred_group;
            predlabel=Stemp.pred_label;
            predtt=repmat(pred_type_traintest(strcmp(S.pred_type,'HGF')),1,length(predlabel));

            % remove trials not in EEG and put in order of EEG data trials
            predtemp=predtemp(tnums,:);

            pred = [pred predtemp];
            pred_group=[pred_group predgroup+length(pred_group)];
            pred_label=[pred_label predlabel];
            pred_format = [pred_format repmat({'cont'},1,length(predlabel))];
            pred_traintest=[pred_traintest predtt];
            clear predtemp predgroup predlabel predtt
        end
    end  

    % remove all-NaN predictors
    ip = all(isnan(pred),1);
    pred(:,ip) = [];
    pred_group(ip)=[];
    pred_label(ip)=[];
     
    %% continuous covariate predictor transforms
    S.prep.calc.pred.transform = 'notrans'; % options: arcsinh, rank or notrans
    S.prep.calc.pred.zscore = 1;
    
    
    % Pred data transformation
    if ~isempty(S.pred_transform)
        contpred=find(strcmp(pred_format,'cont'));
        for pr = contpred
            if strcmp(S.pred_transform,'arcsinh') % need to modify this to apply to only selected predictors
                x=pred(:,pr);
                pred(:,pr)=log(x+sqrt(x.^2+1));
            
            elseif strcmp(S.pred_transform,'rank') % need to modify this to apply to only selected predictors
                x=pred(:,pr);
                [~,pred(:,pr)]=sort(x);
            end
        end
    end
    


    %% PCA on conds/covariates that are multicollinear. PCA results in common variance component and max unique variance components for each
    S.prep.calc.pred.PCA_cov = {'ODD','epsi'}; % will include all covariates that include this text for PCA
    S.prep.calc.pred.PCA_minvar = 0.9; % minimum fraction of variance explained by sparse PCA relative to standard PCA. Standard PCA values returned if not meeting this criterion.
    S.prep.calc.pred.PCA_models = [1:4]; % specifies how many models to run (columns) and how many PCs to estimate for each

    % re-do z-scoring
    S.prep.calc.pred.zscore
    
    
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
        data = downsample(data',S.prep.calc.eeg.dsample)';
        total_samples = downsample(total_samples',S.dsample)';
        select_samples = downsample(S.prep.select.samples',S.dsample)';
    end
    
    % make 3D
    data=reshape(data,size(data,1),[],n_trials);

    % flip channels right to left
    if ~isempty(S.prep.calc.eeg.flipchan)
        load(S.prep.path.inputs.chanlocs);
        flipidx=[];
        for row = 1:length(S.flipchan)
            flipidx = [flipidx find(ismember(eventType,S.cond_idx{S.flipchan(row)}))];
        end
        data(:,:,flipidx) = flipchan(data(:,:,flipidx),chanlocs,S);
    end

    % select samples
    select = dsearchn(total_samples',[select_samples(1),select_samples(end)]');
    data = data(:,select(1):select(2),:);
    
    % subtractions
    if ~isempty(S.prep.pred.factor.subtract)
        ...
    end
    
    % Data transformation
    if ~isempty(S.prep.calc.eeg.transform)
        x=data;
        if strcmp(S.prep.calc.eeg.transform,'arcsinh')
            data=log(x+sqrt(x.^2+1));
        end
    end
    
    % transform to z-score over trials (unit variance) so that statistical coefficients are normalised
    if S.prep.calc.eeg.zscore % 
        data = zscore(data,[],3);
    end
    
    % finally, trim
%     S.prep.calc.eeg.ndec=8; % trim data to a number of decimal places
        


%% separate EEG/predictor data into training and testing fractions, only for testing decoders
S.prep.traintest.type = 'random'; % method to split data into training and testing. options: 'random' selection of trials (S.trainfrac must be >0); 'cond': select by condition
S.prep.traintest.frac = 1; % fraction of trials to use as training data (taken at random). Remainder used for testing. E.g. training for encoding vs. testing for decoding. Or, MPVA training vs. MVPA testing.
S.prep.traintest.balance_conds =1; % random training trials balanced across conditions
S.prep.traintest.num_runs = 1; % run analysis this many times, each time with a different random choice of training/test data


    % create sets of trials split into conditions, according to S.cond_idx
    % this enables subtraction of the mean of one condition from trials
    % of another condition (next step)
    for set = 1:size(S.cond_idx,1)
        setidx{set} = find(ismember(eventType,[S.cond_idx{set}])); % tnum indices organised in sets
        setData{set} = data(:,:,setidx{set}); % data ordered in sets
    end

    % perform subtraction of mean of one condition from trials of
    % another. Useful for isolating mismatch responses, i.e. mismatch
    % trials minus the mean of standard trials.
    if ~isempty(S.row_subtract)
        for sb = 1:length(S.row_subtract{1})
            meandat = mean(setData{S.row_subtract{2}(sb)},3);
            subData{sb} = bsxfun(@minus,setData{S.row_subtract{1}(sb)},meandat);
            subidx{sb}=setidx{S.row_subtract{1}(sb)};
        end
        setData=subData;
        setidx = subidx;
    end

    % pool sets into "contrasts" (con) according to S.contrast_rows
    % this reduces the number of sets down to only those of interest,
    % e.g. by removing the distinction between left and right-sided
    % stimuli by maintaining distinction between mismatch and standard
    % trials (if standards not subtracted already).
    conData={};
    idx={};
    if ~isempty(S.contrast_rows)
        for con = 1:length(S.contrast_rows)
            conData{con} = cat(3,setData{S.contrast_rows{con}}); % data still ordered according to original sub-sets
            idx{con} = [setidx{S.contrast_rows{con}}]; % indices of tnum for each contrast, so we know what trial order the data is in
        end
    else
        conData = {cat(3,setData{:})};
        idx = {[setidx{:}]};
    end


    % duplicate data over a number of runs
    if S.num_runs>1 && length(conData)==1
        for con = 1:S.num_runs
            conData{con} = conData{1};
            idx{con} = idx{1};
        end
    end

    switch S.traintest
        case 'random'
            for con = 1:length(conData)
                % Split into training (encoding) and testing (decoding) fractions
                % Produces S.trainidx and S.testidx. These are indices
                % of conData-ordered trials that will be used for training and
                % testing.
                if S.balance_conds
                    set_events = eventType(idx{con}); % eventypes organised into sets within this contrast (current conData trial order)
                    ucond = unique(set_events);
                    S.trainidx{con} = [];
                    for u = 1:length(ucond)
                        cond_idx = find(set_events==ucond(u)); % for each condition, it's index within conData
                        rng(con); % seed random number generator for consistency
                        S.trainidx{con} = [S.trainidx{con} randsample(cond_idx,round(S.trainfrac*length(cond_idx)))]; 
                    end
                    S.trainidx{con} = sort(S.trainidx{con});
                else
                    rng(1); % seed random number generator for consistency
                    S.trainidx{con} = sort(randsample(length(idx{con}),round(S.trainfrac*length(idx{con}))));
                end
                if S.trainfrac<1
                    S.testidx{con} = 1:length(idx{con});
                    S.testidx{con}(S.trainidx{con}) = [];
                else
                    S.testidx{con} = S.trainidx{con};
                end
            end
        case 'cond'
            % split into training (encoding) and testing (decoding)
            % trials by trial type (condition/contrast)
            if length(idx)~=2
                error('wrong number of contrasts: needs 2')
            end
            conData{1} = cat(3,conData{1},conData{2});
            conData(2)=[]; 
            %pred = cat(1,pred(idx{1},:),pred(idx{2},:));
            S.trainidx{1}= 1:length(idx{1});
            S.testidx{1} = length(idx{1})+1 : length(idx{1})+length(idx{2});
            idx{1} = cat(2,idx{1},idx{2});
            idx(2)=[];

    end

    % each pred need train/test elements
    % 

    % UPDATE SO THAT WE DO NOT SPECIFY TRAINING ANFD TESTING PREDICTORS
    % HERE. THEY ARE SELECTED BY ENCODING AND DECODING FUNCTIONS, WHICH
    % UTILISIED TRAINING AND TESTING FRACTIONS
    % create con-specific predictors and update trainidx to remove NaN
    % predictor value indices
    if c==1
        pred=pred;
    end
    pred_train={};
    pred_test={};
    for con = 1:length(conData)
        % idx is an index of tnum derived from eventType
        pred_train{con}=pred(idx{con},strcmp(pred_traintest,'train'));
        pred_test{con}=pred(idx{con},strcmp(pred_traintest,'test'));
        if isempty(pred_test{con})
            pred_test{con}=pred_train{con};
        end
        % remove NaNs
        pred_nonan_train{con}=find(~any(isnan(pred_train{con}),2));
        pred_nonan_test{con}=find(~any(isnan(pred_test{con}),2));
        %pred{con}=pred{con}(D.prep(d).pred_nonan{con},:);
        S.trainidx{con} = S.trainidx{con}(ismember(S.trainidx{con},pred_nonan_train{con}));
        S.testidx{con} = S.testidx{con}(ismember(S.testidx{con},pred_nonan_test{con}));% if no need to remove nans from testidx, this can be commented out

        % store info so we know which trials and how many were analysed
        stats.subname{d,1}=D.prep(d).subname;
        stats.trialinfo{con}.tnums{d,c}=tnums;
        stats.trialinfo{con}.idx{d,c}=idx{con};
        stats.trialinfo{con}.trainidx{d,c}=S.trainidx{con};
        stats.trialinfo{con}.testidx{d,c}=S.testidx{con};
        stats.trialinfo{con}.ntrials_traintest(d,c,:)= [length(S.trainidx{con}),length(S.testidx{con})];
        stats.trialinfo{con}.pred_train{d,c}= pred_train{con};
        stats.trialinfo{con}.pred_test{d,c}= pred_test{con};

    end

    for con = 1:length(conData)
        % reduce to Global Field Power
        if size(conData{con},1)>1
            gfpData{con} = squeeze(std(conData{con},{},1));
        else
            gfpData{con} = squeeze(conData{con});
        end
        stdgfpData{con} = std(gfpData{con},[],1);
    end
    
    
%% output format
S.prep.output.format = 'group'; % options: 'subjects' separately or 'group' data combined
S.prep.output.save = 0; % save data to disk or not

D.prep(d).sample(s).dtab=dtab;

end

function pred = HGF_covariates(cv,subname,tnums,datfile)
% remove trials not in EEG and later use tnums to put in EEG trial order
if ~any(diff(tnums)>1)
    error('tnums is probably wrong')
end
% get data
HGF=load(cv.input);
D_fit=HGF.D_fit; 
ind = find(strcmp({D_fit(:).subname},subname));
dat = D_fit(ind).HGF;
for ci = 1:length(cv.def)
    if strcmp(cv.def{ci}{1},'RT')
        predtemp = dat.y(:,2);
        pred_label='RT';
    elseif strcmp(cv.def{ci}{1},'traj')
        % for each traj group... reformat into predictor matrix

    end
    pred.(pred_label)=predtemp(tnums,:);

end

    HGF=load(fullfile(S.path.hgf,S.file.hgf{nC})); D_fit=HGF.D_fit; % prevents S loading.
    if ~isfield(D_fit,'dt')
        dtfile = load(fullfile(S.path.hgf,S.file.get_dt));

        pdata = readtable(S.path.datfile);
        subs = pdata.Subject(find(pdata.Include));
        subind = ismember({dtfile.D_fit(:).subname},subs);
%             D_fit = D_fit(subind);
%             D_fit.subname=dtfile.D_fit(subind).subname;
%             D_sub.dt=dtfile.D_fit(subind).dt;
        [D_fit(:).dt]=deal(dtfile.D_fit(subind).dt);
        [D_fit(:).subname]=deal(dtfile.D_fit(subind).subname);
    end
    ind = find(strcmp({D_fit.subname},D.prep(d).subname));
    if isempty(ind)
        disp('USING ONE HGF FOR ALL SUBJECTS')
        ind=1;
    end
    D.prep(d).HGF = D_fit.HGF;
    if length(S.use_hgf)>1
        [Stemp] = HGF_traj2mat(S,D,suffixes{nC});
    else
        [Stemp] = HGF_traj2mat(S,D);
    end
    predtemp=Stemp.pred;
    predgroup=Stemp.pred_group;
    predlabel=Stemp.pred_label;
    predtt=repmat(pred_type_traintest(strcmp(S.pred_type,'HGF')),1,length(predlabel));

    % remove trials not in EEG and put in order of EEG data trials
    predtemp=predtemp(tnums,:);

    pred = [pred predtemp];
    pred_group=[pred_group predgroup+length(pred_group)];
    pred_label=[pred_label predlabel];
    pred_format = [pred_format repmat({'cont'},1,length(predlabel))];
    pred_traintest=[pred_traintest predtt];
    clear predtemp predgroup predlabel predtt
