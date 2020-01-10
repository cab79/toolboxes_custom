function D = eegstats_dataprep(S)
% inputs: S is a settings structure (see example script)
% outputs: D is the prepared EEG data and predictors as a long-format table
% (all subjects) or series of tables (one per subject)

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
GFL.(GFL.func).load.suffix = S.prep.select.suffixes(1);
GFL = getfilelist(GFL);
designmat=cell2table(GFL.(GFL.func).designmat(2:end,:));
designmat.Properties.VariableNames = GFL.(GFL.func).designmat(1,:);
filelist = GFL.(GFL.func).filelist;
subjects = GFL.(GFL.func).select.subjects;
if isempty(filelist)
    error('No files found to import!\n');
end

% run though all EEG files in a loop
for d = 1:length(subjects)
    
    % create output D
    D(d).prep.subname = subjects{d};
        
    % FIND THE FILES FOR THIS SUBJECT
    subfiles = filelist(find(not(cellfun('isempty', strfind(filelist,D(d).prep.subname)))));
    
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
                [eventType,tnums, ~, bnums] = get_markers(temp,'EGI');
            else % from other systems, e.g. Brainamps
                [eventType,tnums, ~, bnums] = get_markers(temp,'BV');
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
    
    if exist('topography','var')
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
                pred_label=S.prep.pred.factor.label{fac};
            else
                pred_label=['fac' num2str(fac)];
            end
            dtab.(pred_label) = pred;
            clear predtemp predlabel
        end
    end
    
    %% Covariates
    CVs=S.prep.pred.covariate;
    for ci = 1:length(CVs)
        switch CVs(ci).type
            case 'HGF'
                pred_out = HGF_covariates(...
                    CVs(1),...
                    D(d).prep.subname,...
                    tnums,...
                    S.prep.path.inputs.subjects);
        end
        repl=find(ismember(pred_out.Properties.VariableNames,dtab.Properties.VariableNames));
        if ~isempty(repl)
            % add suffixes if labels are replicated
            for rp = 1:length(repl)
                pred_out.([pred_out.Properties.VariableNames(repl(rp)) 'a']) = pred_out.(pred_out.Properties.VariableNames(repl(rp)));
                pred_out.(pred_out.Properties.VariableNames(repl(rp)))=[];
            end
        end
        % add pred_out to dtab
        dtab=[dtab pred_out];
        % currently unused:
        CVs(ci).level; % subject, trial
    end

    % remove all-NaN predictors
    vt = vartype('numeric');
    numericvar_names=dtab(:,vt).Properties.VariableNames;
    ip = all(isnan(table2array(dtab(:,vt))),1);
    dtab(:,numericvar_names(ip)) = [];
     
    %% continuous covariate predictor transforms
    % Pred data transformation
    if ~isempty(S.prep.calc.pred.transform)
        for pr = 1:length(numericvar_names)
            if strcmp(S.prep.calc.pred.transform,'arcsinh') % need to modify this to apply to only selected predictors
                x=table2array(dtab(:,numericvar_names));
                dtab(:,numericvar_names)=array2table(log(x+sqrt(x.^2+1)));
            
            elseif strcmp(S.prep.calc.pred.transform,'rank') % need to modify this to apply to only selected predictors
                x=dtab(:,numericvar_names);
                [~,dtab(:,numericvar_names)]=array2table(sort(x));
            end
        end
    end
    
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
        total_samples = downsample(total_samples',S.prep.calc.eeg.dsample)';
        select_samples = downsample(S.prep.select.samples',S.prep.calc.eeg.dsample)';
    end
    
    % make 3D
    data=reshape(data,size(data,1),[],n_trials);

    % flip channels right to left
    if ~isempty(S.prep.calc.eeg.flipchan)
        load(S.prep.path.inputs.chanlocs);
        flipidx = ismember(eventType,S.prep.calc.eeg.flipchan);
        data(:,:,flipidx) = flipchan(data(:,:,flipidx),chanlocs,S.prep.path.inputs.GSNlocs);
    end

    % select samples
    select = dsearchn(total_samples',[select_samples(1),select_samples(end)]');
    data = data(:,select(1):select(2),:);
    
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
    
    
    % trim
%     S.prep.calc.eeg.ndec=8; % trim data to a number of decimal places
        
    %% reshape data into samples and combine over subjects
    D(d).prep.dim = size(data);
    grpdata{d}=reshape(data,[],n_trials)';
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
    D(d).prep.dtab.trial=categorical(tnums');
    D(d).prep.dtab.eventTypes=categorical(eventType');
    D(d).prep.dtab.train=categorical(train_matrix);
    D(d).prep.dtab.test=categorical(test_matrix);
    
end

%% outputs: grouped data
switch S.prep.output.format 
    case 'group'
        % create a single vertical concatenated design matrix
        dim=D(1).prep.dim;
        prep.dtab=D(1).prep.dtab;
        for d = 2:length(D)
            prep.dtab = vertcat(prep.dtab,D(d).prep.dtab);
        end
        dim(3)=height(prep.dtab);
        D=struct;
        D.prep=prep;
        D.prep.dim=dim;
%         D(2:end)=[];
        
        
        % vertically concatenate EEG data over subjects
        temp=grpdata;
        grpdata={};
        grpdata{1}= vertcat(temp{:}); 
    case 'subject'
        % clean up D
        E=struct;
        for d = 1:length(D)
            E(d).prep.dtab = D(d).prep.dtab;
            E(d).prep.dim = D(d).prep.dim;
        end
        D=E;
        clear E
end

%% final transformations: must take place on grouped data (over subjects) if it is grouped
for d = 1:length(D)
    disp(['data transformations: ' num2str(d) '/' num2str(length(D))])
    K=1;
    if ~isempty(S.prep.calc.pred.PCA_cov)
            
        % find columns with relevant predictors
        col_idx=contains(D(d).prep.dtab.Properties.VariableNames,S.prep.calc.pred.PCA_cov);
        col_names = D(d).prep.dtab.Properties.VariableNames(col_idx);
        dtab_PCA=D(d).prep.dtab;
        dtab_PCA(:,col_idx)=[];

        % extract predictors
        X=[];
        for nc = 1:length(col_names)
            X(:,nc)=D(d).prep.dtab.(col_names{nc});
        end
        X=center(X);
        X=normalize(X);

        % PCA parameters
        delta = inf;
        maxiter = 3000;
        convergenceCriterion = 1e-9;
        verbose = false;

        % identify range of L1 norm thresholds to test
        [~,~,A] = svd(X, 'econ');
        AXX=abs((A'*X')*X);
        AXX=sort(AXX(:));
        range = [max(AXX)*2 fliplr(prctile(AXX,[1:100]))];

        % test all and find optimal solution for sparse PCA
        out=struct; 
        nC = intersect(S.prep.calc.pred.PCA_models,1:length(col_names));
        for K = nC
            for rn = 1:length(range)
                stop=range(rn);
                [B,SV,Bsvd,SVsvd] = spca(X, [], K, delta, stop, maxiter, convergenceCriterion, verbose);
                scores = X * B;
                % test collinearity
                rs=corr(scores,'type','Pearson');
                rs=abs(triu(rs,1));
                coli=find(rs>0.8);
                % ensure correct number of PCs is obtained with min variance
                % criterion met
                if any(B,'all') && isempty(coli) && ~any(isnan(rs),'all') && sum(SV)/sum(SVsvd)>=S.prep.calc.pred.PCA_minvar
                    out(K).rn=rn;
                    out(K).B=B;
                    out(K).SV = SV;
                    out(K).scores = scores;
                    break
                end
            end

            % if criterion is not met, use standard PCA
            try out(K).B
            catch
                out(K).rn=0;
                out(K).B=Bsvd;
                out(K).SV = SVsvd;
                out(K).scores = X * Bsvd;
            end

            % add PCA components to data table
            for k = 1:K
                disp(['adding component M' num2str(K) 'PC' num2str(k)])
                dtab_PCA.(['M' num2str(K) 'PC' num2str(k)]) = out(K).scores(:,k);
            end
        end
        D(d).prep.dtab=dtab_PCA;
    end
    
    vt = vartype('numeric');
    numericvar_names=D(d).prep.dtab(:,vt).Properties.VariableNames;
    if S.prep.calc.pred.test_collinearity && K==1
        
        % for pr = 1:length(numericvar_names)
%             if strcmp(S.prep.calc.pred.transform,'arcsinh') % need to modify this to apply to only selected predictors
%                 x=table2array(dtab(:,numericvar_names));
%                 dtab(:,numericvar_names)=array2table(log(x+sqrt(x.^2+1)));
        
        rs=corr(table2array(D(d).prep.dtab(:,numericvar_names)),'type','Pearson');
        rs=abs(triu(rs,1));
        save(fullfile(S.prep.path.outputs,'predictor_correlations.mat'),'rs','col_names');
        coli=find(rs>S.prep.calc.pred.test_collinearity);
        [row,col] = ind2sub(size(rs),coli);
        for ci = 1:length(coli)
            var1=numericvar_names{row(ci)};
            var2=numericvar_names{col(ci)};
            disp(['collinear predictors: ' var1 ', ' var2])
        end
        if ~isempty(coli)
            error('collinear predictors - see above')
        end
    end
    
    % z-scoring on continuous predictors
    if S.prep.calc.pred.zscore
        for pr = 1:length(numericvar_names)
            [zdata, D(d).prep.pred_means(pr), D(d).prep.pred_stds(pr)] = zscore(table2array(D(d).prep.dtab(:,numericvar_names{pr})));
            D(d).prep.dtab(:,numericvar_names{pr}) = array2table(zdata);
        end
    end
    
    % check design matrix is full rank
%     rk = rank(D(d).prep.dtab);
%     disp(['Design matrix rank, subject ' num2str(d) ', = ' num2str(rk)]);
%     if rk==size(D(d).prep.dtab,2)
%         disp('Design matrix is full rank')
%     end
    
    % FOR EVERY SAMPLE OVER ALL SUBJECTS:
    for s = 1:size(grpdata{d},2)
        % duplicate dtab over data samples
%         D(d).prep.Y(s).dtab = D(d).prep.dtab;
        D(d).prep.Y(s).dtab=table;
        
        % add EEG after z-scoring over trials (unit variance) so that statistical coefficients are normalised
        if S.prep.calc.eeg.zscore
            [D(d).prep.Y(s).dtab.data, D(d).prep.Y(s).data_mean, D(d).prep.Y(s).data_std] = zscore(double(grpdata{d}(:,s)));
        else
            D(d).prep.Y(s).dtab.data = double(grpdata{d}(:,s));
        end
    end
    
end

% save data to disk
if S.prep.output.save
    disp('saving data to disk...')
    save(fullfile(S.prep.path.outputs,['eegstats_dtab' S.prep.sname '.mat']),'D','-v7.3')
    disp('...done')
end

function pred = HGF_covariates(cv,subname,tnums,datfile)
% remove trials not in EEG and later use tnums to put in EEG trial order
if ~any(diff(tnums)>1)
    error('tnums is probably wrong')
end
pred=table;
% get data
HGF=load(cv.input);
D_fit=HGF.D_fit; 
if ~isfield(D_fit,'dt')
    dtfile = load(cv.supp);
    pdata = readtable(datfile);
    subs = pdata.Subject(find(pdata.Include));
    subind = ismember({dtfile.D_fit(:).subname},subs);
    [D_fit(:).dt]=deal(dtfile.D_fit(subind).dt);
    [D_fit(:).subname]=deal(dtfile.D_fit(subind).subname);
end
ind = find(strcmp({D_fit(:).subname},subname));
if isempty(ind)
    disp('USING ONE HGF FOR ALL SUBJECTS')
    ind=1;
end
dat = D_fit(ind).HGF;
for ci = 1:length(cv.def)
    if strcmp(cv.def{ci}{1},'RT')
        predtemp = dat.y(:,2);
        pred_label='RT';
    elseif strcmp(cv.def{ci}{1},'traj')
        % for each traj group... reformat into predictor matrix
        D.HGF = dat;
        S.traj = {cv.def{ci}(2:end)}; 
        [Stemp] = HGF_traj2mat(S,D);
        predtemp=Stemp.pred;
        pred_label=Stemp.pred_label;
    end

    % remove trials not in EEG and put in order of EEG data trials
    for pr = 1:length(pred_label)
        pred.(pred_label{pr})=predtemp(tnums,pr);
    end
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
