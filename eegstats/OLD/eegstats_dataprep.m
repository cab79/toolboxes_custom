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
    D(d).prep.tnums = tnums;
    
    if exist('topography','var')
        comps = 1:size(timecourse,1);
        data= topography(:,comps) * timecourse(comps,:);  
    else
        data= timecourse;  
    end
    total_samples=S.prep.original_samples;
    select_samples=S.prep.select.samples;
    
    
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
    
    %% Covariates
    if isfield(S.prep.pred,'covariate')
        CVs=S.prep.pred.covariate;
        for ci = 1:length(CVs)
            switch CVs(ci).type
                case 'HGF'
                    pred_out = HGF_covariates(...
                        CVs(ci),...
                        D(d).prep.subname{:},...
                        tnums,...
                        S.prep.path.inputs.subjects);
                case 'trials'
                    pred_out = trial_covariate(...
                        CVs(ci),...
                        D(d).prep.subname{:},...
                        tnums);
                case 'eegfile'
                    pred_out = eegfile_covariate(...
                        CVs(ci),...
                        D(d).prep.(fename{1})(1).dat,...
                        tnums);
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
    
    % select data
    data = data(:,:,ismember(orig_tnums,tnums));
    n_trials = size(data,3);

    % flip channels right to left
    if ~isempty(S.prep.calc.eeg.flipchan)
        load(S.prep.path.inputs.chanlocs);
        flipidx = ismember(eventType,S.prep.calc.eeg.flipchan);
        data(:,:,flipidx) = flipchan(data(:,:,flipidx),chanlocs,S.prep.path.inputs.GSNlocs);
    end

    % select samples
    if exist('select_samples','var')
        select = dsearchn(total_samples',[select_samples(1),select_samples(end)]');
        data = data(:,select(1):select(2),:);
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
        prep.tnums{1} = D(1).prep.tnums;
        for d = 2:length(D)
            prep.dtab = vertcat(prep.dtab,D(d).prep.dtab);
            prep.tnums{d} = D(d).prep.tnums;
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
%         dtab_PCA(:,col_idx)=[];

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
            Btable = table;
            Btable.cov = dtab_PCA(:,col_idx).Properties.VariableNames';
            for k = 1:K
                disp(['adding component M' num2str(K) 'PC' num2str(k)])
                pcname=['M' num2str(K) 'PC' num2str(k)];
                dtab_PCA.(pcname) = out(K).scores(:,k);
                Btable.(pcname) = out(K).B(:,k);
            end
            D(d).prep.pca{K}=Btable;
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

    
    if S.prep.calc.eeg.cca.on
        
        [U,~,iU]=unique(D(d).prep.dtab.ID,'stable');
        origgrpdata=grpdata;
        newdim = D(d).prep.dim;
        
        for pt = 1:length(S.prep.calc.eeg.cca.type)
            pcatype = S.prep.calc.eeg.cca.type{pt};
            separate = S.prep.calc.eeg.cca.type_sep(pt);
%             if exist('NUM_FAC','var')
%                 NUM_FAC = min([NUM_FAC; S.prep.calc.eeg.cca.numfac{pt}]);
%             else
                NUM_FAC = S.prep.calc.eeg.cca.numfac{pt};
%             end
            newgrpdata{pt}={};
            currdim = newdim;
            
            subdata={};
            dat={};
            switch S.prep.calc.eeg.cca.data_in
                case 'conds'
                    % Average over trials within each condition, per subject.
                    % Currently assumes there are no missing conditions - needs
                    % updating to allow missing conditions.
                    [C,~,iC]=unique(double(D(d).prep.dtab.eventTypes));
                    nreps = C(end);
                    for u=1:length(U)
                        for c=1:length(C)
                            dat{u}(c,:,:) = reshape(mean(grpdata{d}(iU==u & iC==c,:),1),1,currdim(1),currdim(2));
                        end
                        repidx{u} = C;
                    end
                case 'trials'
                    trialidx = unique([D.prep.tnums{:}]); % repetitions
                    nreps = trialidx(end);
                    for u=1:length(U)
                        dat{u} = reshape(grpdata{d}(iU==u,:),[],currdim(1),currdim(2));
                        repidx{u} = D.prep.tnums{u};
                    end
            end
        
%         % baseline correct
%         if any(S.prep.calc.eeg.cca.base)
%             for u=1:length(U)
%                 if exist('select_samples','var')
%                     select = dsearchn(select_samples',[S.prep.calc.eeg.cca.base(1),S.prep.calc.eeg.cca.base(end)]');
%                 else
%                     select = dsearchn(total_samples',[S.prep.calc.eeg.cca.base(1),S.prep.calc.eeg.cca.base(end)]');
%                 end
%                 subdata{u} = bsxfun(@minus,subdata{u},mean(subdata{u}(:,select(1):select(2)),2));
%             end
%         end
            
            if strcmp(pcatype,'spatial')
                obs_dim = 2; % observation dimension
                var_dim = 1; % variable dimension
                % split observations?
                if separate
                    for obs = 1:currdim(2)
                        for u=1:length(U)
                            subdata{obs}{u} = reshape(permute(dat{u}(:,:,obs),[2 1 3]),currdim(1),[]); 
                            subdata{obs}{u} = permute(subdata{obs}{u},[2 1]);
                        end
                    end
%                     newdim(obs_dim)=1;
                else
                    % reshape
                    for u=1:length(U)
                        subdata{1}{u} = reshape(permute(dat{u},[2 1 3]),currdim(1),[]); % chan x trials*time
                        subdata{1}{u} = permute(subdata{1}{u},[2 1]); % trials*time x chan 
                    end
                end
            elseif strcmp(pcatype,'temporal')
                obs_dim = 1; % observation dimension
                var_dim = 2; % variable dimension
                % split observations?
                if separate
                    for obs = 1:currdim(1)
                        for u=1:length(U)
                            subdata{obs}{u} = reshape(dat{u}(:,obs,:),[],currdim(2));
                        end
                    end
%                     newdim(obs_dim)=1;
                else
                    for u=1:length(U)
                        subdata{1}{u} = reshape(dat{u},[],currdim(2)); % trials*chan x time
                    end
                end
            elseif strcmp(pcatype,'both')
                obs_dim = 0; % observation dimension
                var_dim = 0; % variable dimension
                for u=1:length(U)
                    subdata{1}{u} = reshape(dat{u},[],currdim(1)*currdim(2)); % trials x chan*time
                end
            end
          
            % repetitions (usually time or chan, multiplied by conds or trials)
            for u=1:length(U)
                rep{u} = repidx{u};
                if obs_dim && length(subdata)==1
                    for di = 2:newdim(obs_dim)
                        rep{u} = [rep{u}, repidx{u} + (di-1)*nreps];
                    end
                end
            end
            
            %%% PCA/CCA function %%%
            O=struct;
            for o = 1:length(subdata)
                if ~obs_dim || length(subdata)>1
                    obsdim=1;
                else
                    obsdim=newdim(obs_dim);
                end
                [O(o).W,O(o).mdata,O(o).COEFF,O(o).mu,O(o).sigma,NUM_FAC] = perform_CCA(subdata{o},NUM_FAC,S,rep,nreps*obsdim);
                sz=size(O(o).mdata);
                O(o).mdata_avg = squeeze(mean(reshape(O(o).mdata,sz(1),nreps,obsdim,sz(3)),2));
            end
            %newdim(var_dim) = sz(1); %UPDATE

%             
            % get PCA and CCA scores for single trials
            for u=1:length(U)
                if strcmp(pcatype,'spatial')
                    if length(O)==1
                        temp = permute(reshape(grpdata{d}(iU==u,:),[],currdim(1),currdim(2)),[1 3 2]);% trials x time x chan
                        temp = reshape(temp,[],size(temp,3)); 
                    else
                        temp = reshape(grpdata{d}(iU==u,:),[],currdim(1),currdim(2));% trials x chan x comp 
                    end
                elseif strcmp(pcatype,'temporal')
                    if length(O)==1
                        temp = reshape(grpdata{d}(iU==u,:),[],currdim(1),currdim(2));% trials x chan x time
                        temp = reshape(temp,[],size(temp,3)); 
                    else
                        temp = permute(reshape(grpdata{d}(iU==u,:),[],currdim(1),currdim(2)),[1 3 2]);% trials x time x comp
                    end
                elseif strcmp(pcatype,'both')
                    if length(O)==1
                        temp = reshape(grpdata{d}(iU==u,:),[],currdim(1),currdim(2));% trials x chan x time
                        temp = reshape(temp,[],size(temp,2)*size(temp,3)); 
                    else
                        temp = permute(reshape(grpdata{d}(iU==u,:),[],currdim(1),currdim(2)),[1 3 2]);% trials x time x chan
                    end
                end
                
                CCAs{u} = [];
                oind{u}=[];
                for o = 1:length(O)
                    
                    if S.prep.calc.eeg.cca.centre
                        tmu = repmat(O(o).mu{u},size(temp,1),1);
                        O(o).PCAs{u} = (temp(:,:,o)-tmu)*O(o).COEFF{u};
                    else
                        O(o).PCAs{u} = temp(:,:,o)*O(o).COEFF{u}; 
                    end
                    
                    if S.prep.calc.eeg.cca.standardise
                        tsigma = repmat(O(o).sigma{u},size(temp(:,:,o),1),1);
                        O(o).PCAs{u} = O(o).PCAs{u}./tsigma;   
                    end
                    
                    CCAs{u} = [CCAs{u}, O(o).PCAs{u}*O(o).W(:,:,u)];
                    oind{u} = [oind{u}, o*ones(1,size(O(o).PCAs{u},2))];
                    
                end
                if var_dim
%                     newdim(var_dim) = length(O)*NUM_FAC(2); % or size(CCAs,2)
                    newdim(var_dim) = NUM_FAC(2);
                end

                if length(O)>1
                    temp = permute(reshape(CCAs{u},length(repidx{u}),NUM_FAC(2),[]),[1 3 2]); % trials x time/chan x cc
                else
                    temp = reshape(CCAs{u},length(repidx{u}),[],NUM_FAC(2)); % trials x time/chan x cc
                end

                tempsz = [size(temp,2),size(temp,3)];

                if strcmp(pcatype,'spatial')
                    temp = reshape(permute(temp,[1 3 2]),[],tempsz(1)*tempsz(2));
                elseif strcmp(pcatype,'temporal')
                    temp = reshape(temp,[],tempsz(1)*tempsz(2));
                elseif strcmp(pcatype,'both')
                    temp = squeeze(temp);
                end

                newgrpdata{pt}{d}(iU==u,:) = temp;
            end
            grpdata=newgrpdata{pt};
            if ~obs_dim
                newdim([1 2]) = tempsz;
            end
            
            % store outputs
            D(d).prep.CCA(pt).pcatype = pcatype;
            D(d).prep.CCA(pt).O = O;
            D(d).prep.CCA(pt).CCA = CCAs;
            D(d).prep.CCA(pt).oind = oind;
            D(d).prep.CCA(pt).currdim = currdim;
            D(d).prep.CCA(pt).options = S.prep.calc.eeg.cca;
            
            % PLOT
            cidx=floor(linspace(1,NUM_FAC(2),min(NUM_FAC(2),16)));
            for o = 1:length(O)
                figure('name',[pcatype ', obs ' num2str(o)]); clf;
                cci=0; chandatacorr=[];
                for cc = 1:NUM_FAC(2)
                    if strcmp(pcatype,'spatial')
                        
                        % topo data
                        O(o).chandata=[];
                        for u=1:length(U)
                            O(o).chandata(cc,:,u) = O(o).COEFF{u}*O(o).W(:,cc,u);
                        end
                        
                        % PLOTS
                        if ismember(cc,cidx)
                            % waveform plot
                            cci=cci+1;
                            mCCA = mean(O(o).mdata_avg(cc,:,:),3);
                            if length(mCCA) == length(S.prep.select.samples)
                                subplot(4,4*2,cci)
                                hold on
                                for u=1:length(U)
                                    % plot(CCA{u}(:,cc));
                                    plot(O(o).mdata_avg(cc,:,u));
                                end
                                % mean
                                % ccall = cellfun(@(X) X(:,cc), CCA,'UniformOutput',0);
                                % mCCA = mean(horzcat(ccall{:}),2);
                                plot(mCCA,'k','LineWidth',2);
                                hold off
                            end
                        
                            % topoplot
                            cci=cci+1;
                            subplot(4,4*2,cci)
                            mchandata = squeeze(mean(O(o).chandata(cc,:,:),3))';
                            topo = topotime_3D(mchandata,S);
                            pcolor(topo), shading interp; axis off
                            title(['comp ' num2str(cc)])
                            
%                             % topoplot from correlation
%                             cci=cci+1;
%                             subplot(4,4*2,cci)
%                             temp = [];
%                             for u=1:length(U)
%                                 s_subdata = reshape(newgrpdata{1}{d}(iU==u,:),[],40,19);
%                                 temp(1,:,u) = corr(O(o).mdata(cc,repidx{u},u)',s_subdata(:,:,o)); % 3x60x20 x trials(57)xComp(1083=57x19)
%                             end
%                             chandatacorr(cc,:) = mean(temp,3); % mean over subjects
%                             topo = topotime_3D(chandatacorr(cc,:)',S);
%                             pcolor(topo), shading interp; axis off
%                             title(['comp corr ' num2str(cc)])
                        end

                    elseif strcmp(pcatype,'temporal')
                        
                        % waveform data
                        for u=1:length(U)
                            O(o).timedata(cc,:,u) = O(o).COEFF{u}(:,1:NUM_FAC(1))*O(o).W(1:NUM_FAC(1),cc,u);  
                        end
                        
                        % PLOTS
                        if ismember(cc,cidx)
                            % waveform plot
                            cci=cci+1;
                            subplot(4,4*2,cci)
                            hold on
                            for u=1:length(U)
                                plot(O(o).timedata(cc,:,u));
                            end
                            mtimedata = squeeze(mean(O(o).timedata(cc,:,:),2));
                            plot(mtimedata,'k','LineWidth',2);
                            hold off
                            
                            % topoplot
                            cci=cci+1;
                            subplot(4,4*2,cci)
                            mCCA = nanmean(O(o).mdata_avg(cc,:,:),3)';
                            if length(mCCA) == length(S.img.chansel)
                                topo = topotime_3D(mCCA,S);
                                pcolor(topo), shading interp; axis off
                            end
                            title(['comp ' num2str(cc)])
                        end

                    end

                end
                D(d).prep.CCA(pt).O(o).chandatacorr = chandatacorr;
            end
               
            
        end
        
        % final chan-time data
        for o = 1:length(O)
            for cc = 1:NUM_FAC(2)
                temp = [];
                for u=1:length(U)
                    temp(1,:,u) = corr(O(o).mdata(cc,repidx{u},u)',origgrpdata{d}(iU==u,:),'type','Spearman'); % 3x60x20 x trials(57)xComp(1083=57x19)
                end
                temp = nanmean(temp,3);
                st_subdata = reshape(temp,D(d).prep.dim(1),D(d).prep.dim(2));
                chantimecorr(cc,:,:,o) = st_subdata; % mean over subjects
            end
        end
        D(d).prep.CCA(pt).chantimecorr = chantimecorr;
        
        % plot - spatiotemporal maxima
        %close all
        for cc = 1:NUM_FAC(2)
            dat=squeeze(chantimecorr(cc,:,:));
            topo = topotime_3D(dat,S);
            thresh = [nanmean(topo(:))-2*nanstd(topo(:)),nanmean(topo(:))+2*nanstd(topo(:))];
            
            dat=topo;
            dat(dat>thresh(1) & dat<thresh(2)) = nan;
            dat(isnan(dat))=0;
            % connected components
            ccon = bwconncomp(dat,26);
            % remove small clusters
            ClusExtent = cellfun(@length,ccon.PixelIdxList);
            ccon.PixelIdxList(ClusExtent<max(ClusExtent)/10)=[];
            nc = length(ccon.PixelIdxList);
            
            % plot
            figure('name',['component ' num2str(cc)]); 
            ax(1) = subplot(nc+1,2,2); hold on
            plot(squeeze(max(topo,[],[1 2])),'b')
            plot([1,size(topo,3)],[thresh(2),thresh(2)],'b--')
            plot(squeeze(min(topo,[],[1 2])),'r')
            plot([1,size(topo,3)],[thresh(1),thresh(1)],'r--')
            
            % plot maps
            pi=2;
            for cci = 1:nc
                % subscripts
                [~,i] = max(abs(topo(ccon.PixelIdxList{cci})));
                [x,y,z]=ind2sub(size(topo),ccon.PixelIdxList{cci}(i));
                
                % line
                plot(ax(1),[z,z],[0,topo(x,y,z)],'k')
                
                % topo
                pi=pi+1;
                subplot(nc+1,2,pi); hold on
                pcolor(topo(:,:,z)), shading interp; axis off
                scatter(y,x,'k','filled')
                hold off
                
                % waveform
                pi=pi+1;
                subplot(nc+1,2,pi); hold on
                plot(squeeze(topo(x,y,:)),'k');
                if topo(x,y,z)>0
                    col='b';
                else
                    col='r';
                end
                plot([z,z],[0,topo(x,y,z)],col)
                hold off
                ylabel(['z=' num2str(z)])
            end
        end
        
        
        % plot - spatial and temporal variance maxima
        % based on the idea that CCAs are multivariate and maximise
        % variance, except this is variance over topotime rather than over
        % trials.
        %close all
        for cc = 1:NUM_FAC(2)
            dat=squeeze(chantimecorr(cc,:,:));
            topo = topotime_3D(dat,S);
            sz=size(topo);
            
            stmask=[];
            sd=2;
            while isempty(find(stmask))
                sd = sd-0.1;
                % temporal regions of high variance over space
                temp = reshape(topo,[],sz(3));
                t_gfp = nanstd(temp,[],1);
                tthresh = mean(t_gfp)+sd*std(t_gfp);

                % spatial regions of high variance over time
                s_gfp = nanstd(topo,[],3);
                s_gfp(isnan(topo(:,:,1))) = nan;
                sthresh = nanmean(s_gfp(:))+sd*nanstd(s_gfp(:));

                % spatiotemporal mask
                tmask = permute(repmat(t_gfp>=tthresh,sz(1),1,sz(2)),[1 3 2]);
                smask = repmat(s_gfp>=sthresh,1,1,sz(3));
                stmask = smask.*tmask;
            end
            
            dat=stmask;
            dat(dat>thresh(1) & dat<thresh(2)) = nan;
            dat(isnan(dat))=0;
            % connected components
            ccon = bwconncomp(dat,26);
            % remove small clusters
            ClusExtent = cellfun(@length,ccon.PixelIdxList);
            ccon.PixelIdxList(ClusExtent<max(ClusExtent)/10)=[];
            nc = length(ccon.PixelIdxList);
            
            % plot
            figure('name',['component ' num2str(cc)]); 
            pi=0;
            % s_gfp
            pi=pi+1;
            subplot(nc+1,2,pi); hold on
            pcolor(s_gfp), shading interp; axis off
            hold off
            title('temporal std')
            % t_gfp
            pi=pi+1;
            subplot(nc+1,2,pi); hold on
            plot(t_gfp,'k');
            plot([1,length(t_gfp)],[tthresh,tthresh],'b--')
            hold off
            title('spatial std')
            for cci = 1:nc
                % subscripts
                [x,y,z]=ind2sub(sz,ccon.PixelIdxList{cci});
                max_zi = find(t_gfp==max(t_gfp(z)));
                max_xyi = find(s_gfp==max(s_gfp(x,y),[],'all'));
                
                % topo
                pi=pi+1;
                subplot(nc+1,2,pi); hold on
                pcolor(topo(:,:,max_zi)), shading interp; axis off
                scatter(y,x,'k')
                hold off
                title('component topography')
                
                % waveform
                pi=pi+1;
                subplot(nc+1,2,pi); hold on
                [xw,yw] = ind2sub(sz([1 2]),max_xyi);
                plot(squeeze(topo(xw,yw,:)),'k');
                plot([max_zi,max_zi],[0,topo(xw,yw,max_zi)],'b')
                hold off
                ylabel(['z=' num2str(max_zi)])
                title('component waveform')
            end
        end
        
        D(d).prep.dim = newdim;


%         % check data looks ok
%         condmean=[];
%         for u=1:length(U)
%             for c=1:length(allrep)
%                 condmean(:,:,c,u) = reshape(mean(grpdata{d}(iU==u & iC==c,:),1),NUM_FAC(2),dim(2));
%             end
%         end
%         figure(); clf;
%         cci=0;
%         cidx=floor(linspace(1,NUM_FAC(2),min(NUM_FAC(2),16)));
%         for cc = 1:NUM_FAC(2)
%             if ismember(cc,cidx)
%                 cci=cci+1;
%                 subplot(4,4,cci)
%                 plot(squeeze(condmean(cc,:,1,:)));
%             end
%         end

        
%         % OLD
%         newgrpdata{d}=nan(size(grpdata{d},1),size(subdata,1));
%         for u=1:length(U)
%             zz=subdata(:,:,u)*AA{u};
%             newgrpdata{d}(iU==u,cidx) = zz;
%         end
%         grpdata=newgrpdata;
        
    end
    
    % FOR EVERY SAMPLE OVER ALL SUBJECTS:
    for s = 1:size(grpdata{d},2)
        % duplicate dtab over data samples
%         D(d).prep.Y(s).dtab = D(d).prep.dtab;
        D(d).prep.Y(s).dtab=table;
        
        % add EEG after z-scoring over trials (unit variance) so that statistical coefficients are normalised
        if S.prep.calc.eeg.zscore
            D(d).prep.Y(s).data_mean = nanmean(double(grpdata{d}(:,s)));
            D(d).prep.Y(s).data_std = nanstd(double(grpdata{d}(:,s)));
            D(d).prep.Y(s).dtab.data = (double(grpdata{d}(:,s)) - D(d).prep.Y(s).data_mean) / D(d).prep.Y(s).data_std;
%             [D(d).prep.Y(s).dtab.data, D(d).prep.Y(s).data_mean, D(d).prep.Y(s).data_std] = zscore(double(grpdata{d}(:,s)));
        else
            D(d).prep.Y(s).dtab.data = double(grpdata{d}(:,s));
        end
    end
    
end

% save data to disk
if S.prep.output.save
    disp('saving data to disk...')
    save(fullfile(S.prep.path.outputs,['eegstats_dtab' S.prep.sname '.mat']),'D','S','-v7.3')
    disp('...done')
end

function pred = trial_covariate(cv,subname,tnums)
% currently works with design info files from CORE study (SCIn outputs)
if ~any(diff(tnums)>1)
    error('tnums is probably wrong')
end
pred=table;
% get data
dtfile = load(cv.input);
subind = ismember({dtfile.D_fit(:).subname},subname);
conds = dtfile.D_fit(subind).dt.design(2,:);
predtemp=[];
for ci = 1:length(cv.def)
    condidx = find(ismember(conds,cv.def{ci}));
    % get consecutive trials
    difftrial = [0 diff(condidx)==1];
    consec=1;
    for df=2:length(difftrial)
        if difftrial(df)==1
            consec(df) = consec(df-1)+1;
        else
            consec(df) = 1;
        end
    end
    predtemp(condidx,1)=consec;
end
% remove trials not in EEG and put in order of EEG data trials
pred.trial=predtemp(tnums);

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


function pred = eegfile_covariate(cv,dat,tnums)
% imports from substructure within eeg data file. Designed for Andrej's
% painloc experiment.
pred=table;
for ci = 1:size(cv.def,1)
    ncov = size(dat.(cv.def{ci}{1}),2);
    for nc = 1:ncov
        predtemp = dat.(cv.def{ci}{1})(:,nc);
        pred.(cv.def{ci,2}{nc})=predtemp(tnums);
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

function [W,mdata,COEFF,mu,sigma,NUM_FAC] = perform_CCA(dataAll,NUM_FAC,S,rep,nobs)
% feed in index of 2nd dimension (trials or time)
               
TEMP_FAC=repmat(size(dataAll{1},2),1,2);
if isempty(NUM_FAC)
    NUM_FAC([1,2]) = TEMP_FAC;
else
    NUM_FAC([1,2]) = min([NUM_FAC;TEMP_FAC]);
end

maxmat = floor((S.prep.calc.eeg.cca.maxGig*1e9)^(1/2));
if size(dataAll,2)>maxmat
    error(['downsample the data by ' num2str(size(dataAll,2)/maxmat) ' to enable CCA'])
end

% centre data
cdata={};
for i = 1:length(dataAll)
    tmpscore = squeeze(dataAll{i});
    mu{i} = mean(tmpscore);
    cdata{i} = tmpscore - repmat(mu{i},size(tmpscore,1),1);
end


% apply MCCA
switch S.prep.calc.eeg.cca.method

    case 'FA'
        % find nfac for each subject
        for i = 1:length(cdata)
            FactorResults = erp_pca(cdata{i},NUM_FAC(1));
            randResults = erp_pca(randn(size(cdata{i})),NUM_FAC(1));
            randResults.scree = randResults.scree*(sum(FactorResults.scree)/sum(randResults.scree)); %scale to same size as real data
            nfac_temp = find(FactorResults.scree<randResults.scree);
            nfac(i) = min(NUM_FAC(1),nfac_temp(1)-1);
        end
        % take median nfac
        NUM_FAC(1) = floor(median(nfac));
        for i = 1:length(cdata)
            FactorResults = erp_pca(cdata{i},NUM_FAC(1));
            COEFF{i} = FactorResults.FacCof;
            W(:,:,i) = COEFF{i}';
            score{i} = FactorResults.FacScr;
        end
    case 'PCA'
       for i = 1:length(cdata)
            [COEFF{i}, score{i},~,~,explained] = pca(cdata{i},'NumComponents',NUM_FAC(1),'Centered',false,'Algorithm','eig');
            [COEFFrand{i}, scorerand{i},~,~,explainedrand] = pca(randn(size(cdata{i})),'NumComponents',NUM_FAC(1));
            nfac_temp = find(explained<explainedrand);
            nfac(i) = min(NUM_FAC(1),nfac_temp(1)-1);
       end
        % take median nfac
        NUM_FAC(1) = floor(median(nfac));
        for i = 1:length(cdata)
            COEFF{i} = COEFF{i}(:,1:NUM_FAC(1));
            W(:,:,i) = COEFF{i}';
            score{i} = score{i}(:,1:NUM_FAC(1));
        end
    case 'eigs'
       for i = 1:length(cdata)
            [COEFF{i},values]=eigs(cdata{i},NUM_FAC(1));
            explained = diag(values)/sum(diag(values));
            [COEFFrand{i},valuesrand]=eigs(randn(size(cdata{i})),NUM_FAC(1));
            explainedrand = diag(valuesrand)/sum(diag(valuesrand));
            nfac_temp = find(explained<explainedrand);
            nfac(i) = min(NUM_FAC(1),nfac_temp(1)-1);
       end
        % take median nfac
        NUM_FAC(1) = floor(median(nfac));
        for i = 1:length(cdata)
            COEFF{i} = COEFF{i}(:,1:NUM_FAC(1));
            W(:,:,i) = COEFF{i}';
            score{i} = cdata{i}*COEFF{i}(:,1:NUM_FAC(1));
        end
end
% obtain subject-specific PCAs
PCdata = nan(NUM_FAC(1),nobs);
sigma = {};
for i = 1:length(cdata)
    if S.prep.calc.eeg.cca.standardise
       sigma{i} = sqrt(var(score{i})); 
       score{i} = score{i}./repmat(sqrt(var(score{i})),size(score{i},1),1); % standardise so that each PC from each subject contributes equally to CCA solution
    end
    PCdata(:,rep{i},i) = score{i}';
end

% apply M-CCA
% https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.23689
reg=1; r= 0.2;
NUM_FAC(2) = min(NUM_FAC(2),NUM_FAC(1));
[W,mdata] = mccas(PCdata,NUM_FAC(2),reg,r,W,S);

function [W, mdata, mWeights] = mccas(data,K,reg,r,weights,Sx)
%% regularized-multiset CCA
%   Date           Programmers               Description of change
%   ====        =================            =====================
%  09/10/2016     Qiong Zhang                 Original code
%% Citation
%  Zhang,Q., Borst,J., Kass, R.E., & Anderson, J.A. (2016) Between-Subject
%  Alignment of MEG Datasets at the Neural Representational Space. 

%% INPUT
% data - input training data to CCA (samples * sensors * subjects)
% K - number of CCAs to retain
% reg - add regularization (1) or not (0)
% r - degree of regularization added
% weights - PCA weight maps for each subject

%% OUTPUT
% W - CCA weights that transform PCAs to CCAs for each subject (PCAs * CCAs * subjects)
% mdata - transformed averaged data (CCAs * sapmles * subjects)
% mWeights - projection weights from sensor to CCAs (CCAs * sensors * subjects)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim = size(data,1);
num = size(data,2);
sub = size(data,3);
dim2 = size(weights,1);
num2 = size(weights,2);
sub2 = size(weights,3);
if(K>dim)
    K = dim;
end
temp=[];
for i=1:sub
    temp=[temp;data(:,:,i)];
end
R=cov(temp.','partialrows'); % cab: added partialrows to allow missing data (nan)
S = zeros(size(R));
for i = 1:dim*sub
    tmp = ceil(i/dim);
    S((tmp-1)*dim+1:tmp*dim,i) = R((tmp-1)*dim+1:tmp*dim,i); 
end

% add regularization 
if(reg==1)    
    if(K>dim2)
        K = dim2;
    end
    temp=[];
    for i=1:sub2
        temp=[temp;weights(:,:,i)];
    end
    R2=cov(temp.');
    S2 = zeros(size(R2));
    for i = 1:dim2*sub2;
        tmp = ceil(i/dim2);
        S2((tmp-1)*dim2+1:tmp*dim2,i) = R2((tmp-1)*dim2+1:tmp*dim2,i); 
    end
    R = R + r*R2;
    S = S + r*S2;
end


% obtain CCA solutions 
%for k = K:-1:1
    [tempW,values]=eigs((R-S),S,K);
%    var_explained = diag(values)/sum(diag(values)); % ARE THEY SORTED?
%    cum_var = cumsum(var_explained);
%    facs = find(var_explained>Sx.prep.calc.eeg.cca.minvar & cum_var>Sx.prep.calc.eeg.cca.frac_explained);
%    if ~isempty(facs)
%        break
%    end
%end

W = zeros(dim,K,sub);
for i=1:sub
    if Sx.prep.calc.eeg.cca.normalise_weights
        W(:,:,i)=tempW((i-1)*dim+1:i*dim,:)./norm(tempW((i-1)*dim+1:i*dim,:));
    else
        W(:,:,i)=tempW((i-1)*dim+1:i*dim,:);
    end
end

% projected data
mdata = zeros(K,num,sub);
for i = 1:sub
    mdata(:,:,i) = W(:,:,i)'*data(:,:,i);
    for j = 1:K
        if Sx.prep.calc.eeg.cca.normalise_weights
            mdata(j,:,i) = mdata(j,:,i)/norm(mdata(j,:,i));
        else
            mdata(j,:,i) = mdata(j,:,i);
        end
    end
end


% projection weights
mWeights = 0;
if(reg==1)
    mWeights = zeros(K,num2,sub2);
    for i = 1:sub2
        mWeights(:,:,i) = W(:,:,i)'*weights(:,:,i);
    end
end

%% OLD
%             case 'PCA'

%                 cca_dat=subdata(:,:); % concatenate channelwise
% 
%                 maxmat = floor((S.prep.calc.eeg.cca.maxGig*1e9)^(1/2));
%                 if size(cca_dat,2)>maxmat
%                     error(['downsample the data by ' num2str(size(cca_dat,2)/maxmat) ' to enable CCA'])
%                 end

%                 NUM_FAC = S.prep.calc.eeg.cca.numfac(1);
%                 cca_dat = zscore(cca_dat')';
%                 CCC=cca_dat'*cca_dat;
%                 [~,score]=nt_mcca(CCC,nChan,[],NUM_FAC);
%                 % keep variance outliers only
%                 keep = find(score>S.prep.calc.eeg.cca.minvar*100);
%                 Nkeep=max(keep);
%                 [A,score]=nt_mcca(CCC,nChan,[],Nkeep);
%                 z=cca_dat*A;
%                 
%                 % CCA: separating into sets
%                 AA=[];
%                 N=nChan;
%                 nblocks=size(CCC,1)/N;
%                 for iBlock=1:nblocks
%                     AA{iBlock}=A(N*(iBlock-1)+(1:N),:);
%                 end
%                 
%                 figure(); clf
%                 hold on
%                 plot(score, '.-');
%                 plot(score(keep), 'r.-');
%                 title('CCA: variance per SC');
%                 ylabel('variance'); xlabel('SC');
%                 
%                 figure(); clf;
%                 cci=0;
%                 cidx=floor(linspace(1,NUM_FAC,min(NUM_FAC,16)));
%                 for cc = 1:NUM_FAC
%                     if ismember(cc,cidx)
%                         cci=cci+1;
%                         subplot(4,4*2,cci)
% %                         yyaxis left
%                         hold on
%                         for u=1:length(U)
%                             subdat = zscore(subdata(:,:,u));
%                             zz(:,:,u)=subdat*AA{u};
%                             plot(zz(:,cc,u));
%                         end
%                         title(['comp ' num2str(cc)])
%                         % mean
% %                         plot(mean(zz(:,cc,:),3),'k','LineWidth',2);
% %                         yyaxis right
%                         plot(z(:,cc),'k','LineWidth',2);
%                         hold off
%                         % topo
%                         cci=cci+1;
%                         subplot(4,4*2,cci)
%                         mchandata = mean(reshape(A(:,cc),[],length(U)),2);
%                         topo = topotime_3D(mchandata,S);
%                         pcolor(topo), shading interp; axis off
%                         title(['comp ' num2str(cc)])
%                     end
%                 end
%                 
%                 % multiply single-trial data by CCA coeff
%                 for u=1:length(U)
%                     temp = permute(reshape(grpdata{d}(iU==u,:),[],dim(1),dim(2)),[1 3 2]);
%                     temp = reshape(temp,[],size(temp,3)); % trials*time x chan
%                     temp = zscore(temp);
%                     temp = temp*AA{u}; % change chans to CCs
%                     temp = reshape(temp,[],dim(2),NUM_FAC); % trials x time x cc
%                     
%                     % baseline correct
%                     if any(S.prep.calc.eeg.cca.base)
%                         if exist('select_samples','var')
%                             select = dsearchn(select_samples',[S.prep.calc.eeg.cca.base(1),S.prep.calc.eeg.cca.base(end)]');
%                         else
%                             select = dsearchn(total_samples',[S.prep.calc.eeg.cca.base(1),S.prep.calc.eeg.cca.base(end)]');
%                         end
%                         for nf = 1:size(temp,3)
%                             temp(:,:,nf) = bsxfun(@minus,temp(:,:,nf),mean(temp(:,select(1):select(2),nf),2));
%                         end
%                     end
%                     
%                     temp = reshape(permute(temp,[1 3 2]),[],NUM_FAC*dim(2));
%                     newgrpdata{d}(iU==u,:) = temp;
%                 end
%                 grpdata=newgrpdata;
%                 D(d).prep.dim(1) = NUM_FAC;
%                 
%                 % check data looks ok
%                 condmean=[];
%                 for u=1:length(U)
%                     for c=1:length(C)
%                         condmean(:,:,c,u) = reshape(mean(grpdata{d}(iU==u & iC==c,:),1),NUM_FAC,dim(2));
%                     end
%                 end
%                 figure(); clf;
%                 for fc = 1:size(condmean,1)
%                     subplot(size(condmean,1),1,fc)
%                     plot(squeeze(condmean(fc,:,1,:)))
%                 end
                