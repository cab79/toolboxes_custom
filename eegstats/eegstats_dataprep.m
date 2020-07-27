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
if isfield(S.prep.select,'prefixes') && ~isempty(S.prep.select.prefixes); GFL.(GFL.func).load.prefix = S.prep.select.prefixes; end
if isfield(S.prep.select,'suffixes') && ~isempty(S.prep.select.suffixes); GFL.(GFL.func).load.suffix = S.prep.select.suffixes; end
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
        prep.grpdata{1} = D(1).prep.data;
        for d = 2:length(D)
            prep.dtab = vertcat(prep.dtab,D(d).prep.dtab);
            prep.tnums{d} = D(d).prep.tnums;
            prep.grpdata{d} = D(d).prep.data;
        end
        dim(3)=height(prep.dtab);
        D=struct;
        D.prep=prep;
        D.prep.dim=dim;
%         D(2:end)=[];
        
        
        % vertically concatenate EEG data over subjects
        temp=prep.grpdata;
        D.prep = rmfield(D.prep,'grpdata');
        D.prep.grpdata{1} = vertcat(temp{:}); 
        clear grpdata temp
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
    
    % z-scoring on continuous predictors, but not if doing PLS on EEG data
    if S.prep.calc.pred.zscore
        %if S.prep.calc.eeg.pca.on == 1 && any(strcmp(S.prep.calc.eeg.pca.PCAmethod,'PLS'))
        %    disp('no z-scoring to properly rescale Y variables for PLS')
        %else
            vt = vartype('numeric');
            numericvar_names=dtab(:,vt).Properties.VariableNames;
            for pr = 1:length(numericvar_names)
                [zdata, D(d).prep.pred_means(pr), D(d).prep.pred_stds(pr)] = zscore(table2array(D(d).prep.dtab(:,numericvar_names{pr})));
                D(d).prep.dtab(:,numericvar_names{pr}) = array2table(zdata);
            end
        %end
    elseif S.prep.calc.eeg.pca.on == 1
        error('must z-score data prior to PCA')
    end
    
    if S.prep.calc.pred.PCA_on
            
        D(d).prep.pred_PCA = predictor_PCA(D(d).prep.dtab,S);

        % add PCA components to data table
        if S.prep.calc.pred.PCA_model_add2dtab
            ym=S.prep.calc.pred.PCA_model_add2dtab;
            for k = 1:size(D(d).prep.pred_PCA.Yscore{ym},2)
                pcname=['PC' num2str(k)];
                D(d).prep.dtab.(pcname) = D(d).prep.pred_PCA.Yscore{ym}(:,k);
            end
        end
    end
    
    if S.prep.calc.pred.test_collinearity
        
        % variables
        var_names = D(d).prep.dtab.Properties.VariableNames;
        var_names = setdiff(var_names,{'ID','group','test','train','eventTypes'},'stable');
        vt = vartype('numeric');
        numericvar_names=D(d).prep.dtab(:,vt).Properties.VariableNames;
        numericvar_ind = find(ismember(var_names,numericvar_names));
        nonnumericvar_ind = find(~ismember(var_names,numericvar_names));
        all_ind = [numericvar_ind,nonnumericvar_ind]; % must be in this order
        
        % combinations
        comb = nchoosek(all_ind,2);   
        comb = comb(ismember(comb(:,1),numericvar_ind),:); 
        corrmat = diag(zeros(1,length(var_names)));
        for ci = 1:length(comb)
            formula = [var_names{comb(ci,1)} '~' var_names{comb(ci,2)} '+(1|ID)'];
            lme = fitlme(D(d).prep.dtab,formula);
            corrmat(comb(ci,1),comb(ci,2)) = sign(double(lme.Coefficients(2,2)))*lme.Rsquared.Ordinary;
        end
        corrmat = triu(corrmat + corrmat',1);
        figure;imagesc(corrmat);colorbar;colormap('jet'); caxis([-1 1])
        plot_names = var_names;
        if ~isempty(S.prep.calc.pred.varname_parts_plot)
            for pn = 1:length(plot_names)
                nparts = strsplit(plot_names{pn},'_');
                if max(S.prep.calc.pred.varname_parts_plot) <= length(nparts)
                    nparts = nparts(S.prep.calc.pred.varname_parts_plot);
                    plot_names{pn} = strjoin(nparts,'_');
                end
            end
            S.prep.calc.pred.varname_parts_plot;
        end
        xticks(1:length(var_names)); xticklabels(plot_names); xtickangle(90);
        yticks(1:length(var_names)); yticklabels(plot_names);
        title('collinearity of predictors')
        
        save(fullfile(S.prep.path.outputs,'predictor_correlations.mat'),'corrmat','var_names');
        coli=find(corrmat>S.prep.calc.pred.test_collinearity);
        [row,col] = ind2sub(size(corrmat),coli);
        for ci = 1:length(coli)
            var1=var_names{row(ci)};
            var2=var_names{col(ci)};
            disp(['collinear predictors: ' var1 ', ' var2])
        end
        if S.prep.calc.pred.allow_collinearity==false && ~isempty(coli)
            error('collinear predictors - see above')
        end
        
        if 0 % old version: continuous predictors only
            vt = vartype('numeric');
            numericvar_names=D(d).prep.dtab(:,vt).Properties.VariableNames;

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
    end
    
    
    % check design matrix is full rank
%     rk = rank(D(d).prep.dtab);
%     disp(['Design matrix rank, subject ' num2str(d) ', = ' num2str(rk)]);
%     if rk==size(D(d).prep.dtab,2)
%         disp('Design matrix is full rank')
%     end

    
    if S.prep.calc.eeg.pca.on
        
        %% use existing data or calculate instead?
        if isfield(S.prep.calc.eeg.pca,'use_existing_data') && ~isempty(S.prep.calc.eeg.pca.use_existing_data)
            Dy = load(S.prep.calc.eeg.pca.use_existing_data,'D');
            D.prep.PCA = Dy.D.prep.PCA;
            D.prep.Y = Dy.D.prep.Y;
        else
        
            %% inputs for PLS/CCA
            Sf.PCAmethod = S.prep.calc.eeg.pca.PCAmethod;
            Sf.type = S.prep.calc.eeg.pca.type;
            Sf.type_obs_sep = S.prep.calc.eeg.pca.type_sep; % analyse "observation" separately. E.g. if temporal PCA, apply separately to each spatial channel?
            Sf.type_numfac = S.prep.calc.eeg.pca.numfac; % max number of factors
            Sf.data_in = S.prep.calc.eeg.pca.data_in; % type of data to analyse: trials or averaged conditions
            Sf.centre_output = S.prep.calc.eeg.pca.centre;
            Sf.standardise_output = S.prep.calc.eeg.pca.standardise;
            Sf.pca_FA_type = S.prep.calc.eeg.pca.FA_type;
            Sf.maxGig = S.prep.calc.eeg.pca.maxGig;
            Sf.normalise_cca_weights = S.prep.calc.eeg.cca.normalise_weights;
            Sf.img=S.img;
            Sf.Y_use4output = S.prep.calc.eeg.pls.Y_use4output;
            Sf.select_ncomp = S.prep.calc.eeg.pca.select_ncomp;
            Sf.cca_reg_weights = S.prep.calc.eeg.cca.reg_weights;
            Sf.cca_test_nPCAcomp = S.prep.calc.eeg.cca.test_nPCAcomp;
            Sf.cca_method = S.prep.calc.eeg.cca.method;
            Sf.cca_reduce_by_scree = S.prep.calc.eeg.cca.reduce_by_scree;
            Sf.cca_select_nPCA_per_group = S.prep.calc.eeg.cca.select_nPCA_per_group;
            Sf.cca_FA_type = S.prep.calc.eeg.cca.FA_type;

            % PLS
            Y={};
            if any(strcmp(Sf.PCAmethod,'PLS'))
                if S.prep.calc.eeg.pca.on
                    Y = D(d).prep.pred_PCA.Yscore;
                    col_names = D(d).prep.pred_PCA.col_names;
                else
                    error('run PCA on predictors prior to PLS')
                end
            end

            % RUN function
            if any(strcmp(Sf.PCAmethod,'PLS')) && S.prep.calc.eeg.pls.Y_select % select a model and output that single model
                [D(d)]=eegstats_components_analysis(D(d),Sf,{Y,col_names});
            elseif any(strcmp(Sf.PCAmethod,'PLS')) % output all models
                D_orig=D;
                for ym = 1:length(Y)
                    [Dtemp]=eegstats_components_analysis(D_orig(d),Sf,{Y(ym),col_names(ym)});
                    D(d).prep(ym) = D(d).prep(1);
                    D(d).prep(ym).grpdata=Dtemp.prep.grpdata;
                    D(d).prep(ym).PCA=Dtemp.prep.PCA;
                    D(d).prep(ym).dim=Dtemp.prep.dim;
                    close all
                    mse_temp=mean(cat(3,D(d).prep(ym).PCA.O.PLS.MSE{:}),3); 
                    msecv_temp=mean(cat(3,D(d).prep(ym).PCA.O.PLS.MSE_CV{:}),3); 
                    pvar_temp=mean(cat(3,D(d).prep(ym).PCA.O.PLS.explained{:}),3); 
                    mse_plot(ym,:) = mse_temp(:,end);
                    msecv_plot(ym,:) = msecv_temp(:,end);
                    pvar_plot(ym,:) = sum(pvar_temp,2);
                end
                figure
                subplot(2,3,1); bar(mse_plot(:,1)); title('X MSE per model')
                subplot(2,3,2); bar(msecv_plot(:,1)); title('X CV MSE per model')
                subplot(2,3,3); bar(pvar_plot(:,1)); title('X % var explained per model')
                subplot(2,3,4); bar(mse_plot(:,2)); title('Y MSE per model')
                subplot(2,3,5); bar(msecv_plot(:,2)); title('Y CV MSE per model')
                subplot(2,3,6); bar(pvar_plot(:,2)); title('Y % var explained per model')
                clear Dtemp D_orig
            else
                [D(d)]=eegstats_components_analysis(D(d),Sf);
            end
        end
    else
        for ip = 1:length(D(d).prep)
            D(d).prep(ip).grpdata = D(d).prep(ip).grpdata{1};
        end
    end
    % FOR EVERY SAMPLE OVER ALL SUBJECTS:
    if ~exist('Dy','var') %|| ~isfield(Dy.prep,'Y') % if Y not added in earlier from existing data file
        for ip = 1:length(D(d).prep)
            for s = 1:size(D(d).prep(ip).grpdata,2)
                % duplicate dtab over data samples
        %         D(d).prep(ip).Y(s).dtab = D(d).prep(ip).dtab;
                D(d).prep(ip).Y(s).dtab=table;

                % add EEG after z-scoring over trials (unit variance) so that statistical coefficients are normalised
                if S.prep.calc.eeg.zscore
                    D(d).prep(ip).Y(s).data_mean = nanmean(double(D(d).prep(ip).grpdata(:,s)));
                    D(d).prep(ip).Y(s).data_std = nanstd(double(D(d).prep(ip).grpdata(:,s)));
                    D(d).prep(ip).Y(s).dtab.data = (double(D(d).prep(ip).grpdata(:,s)) - D(d).prep(ip).Y(s).data_mean) / D(d).prep(ip).Y(s).data_std;
        %             [D(d).prep(ip).Y(s).dtab.data, D(d).prep(ip).Y(s).data_mean, D(d).prep(ip).Y(s).data_std] = zscore(double(D.prep.grpdata(:,s)));
                else
                    D(d).prep(ip).Y(s).dtab.data = double(D(d).prep(ip).grpdata(:,s));
                end
            end
        end
    end
    
end

% save data to disk
if S.prep.output.save
    disp('saving data to disk...')
    if isfield(D.prep,'grpdata')
        D.prep = rmfield(D.prep,'grpdata'); % save memory
    end
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
if ~isempty(cv.def)
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
else
    pred.trial=tnums';
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
% remove almost identical data
dat = corr(table2array(pred));
rm_dat = find(dat(1,2:end)>0.99999);
if ~isempty(rm_dat)
    pred(:,rm_dat+1) = [];
end
% remove data with almost no variance
dat = std(table2array(pred));
rm_dat = dat<0.00001;
if any(rm_dat)
    pred(:,rm_dat) = [];
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



function out = predictor_PCA(dtab,S)
% dtab: all predictors in table format
% out: PCA coefficients etc.
% S: which predictors to include, which PCA method, etc.

% select predictors Y
for ym = 1:length(S.prep.calc.pred.PCA_cov) % for each Y model

    % find columns with relevant predictors
    col_idx=contains(dtab.Properties.VariableNames,S.prep.calc.pred.PCA_cov{ym});
    out.col_names{ym} = dtab.Properties.VariableNames(col_idx);

    for nc = 1:length(out.col_names{ym})
        Y{ym}(:,nc)=dtab.(out.col_names{ym}{nc});
    end
end
out.Y = Y;

% select PCA method
switch S.prep.calc.pred.PCA_type

    case 'PCA'
       for ym = 1:length(Y) % for each Y model
           
            % inital num components
            ncomp=size(Y{ym},2); 
            
            % PCA
            [Ycoeff{ym}, Yscore{ym},~,~,explained] = pca(Y{ym},'NumComponents',ncomp,'Centered',false,'Algorithm','eig');
            
            % select num components
            if S.prep.calc.pred.PCA_minvar
                nfac_temp = find(cumsum(explained) > S.prep.calc.pred.PCA_minvar*100);
                if isempty(nfac_temp)
                    nfac_temp = ncomp;
                end
                nfac = nfac_exp(1);
            else
                [~,~,~,~,explainedrand] = pca(randn(size(Y{ym})),'NumComponents',ncomp,'Centered',false,'Algorithm','eig');
                nfac_temp = find(explained<explainedrand);
                nfac = nfac_temp(1)-1;
            end
            out.Ycoeff{ym} = Ycoeff{ym}(:,1:nfac);
            out.Yscore{ym} = Yscore{ym}(:,1:nfac);
       end
    
    case 'sPCA'
        for ym = 1:length(Y) % for each Y model
           
            % inital num components
            ncomp=size(Y{ym},2); 
            
            % sPCA settings
            delta = inf;
            maxiter = 1000;
            convergenceCriterion = 1e-9;
            verbose = false;

            % sPCA
            [Ycoeff{ym}, values] = spca(Y{ym}, [], ncomp, delta, -ncomp, maxiter, convergenceCriterion, verbose);
            explained = diag(values)/sum(diag(values));
            
            % select num components
            if S.prep.calc.pred.PCA_minvar
                nfac_temp = find(cumsum(explained) > S.prep.calc.pred.PCA_minvar*100);
                if isempty(nfac_temp)
                    nfac_temp = ncomp;
                end
                nfac = nfac_exp(1);
            else
                [~, values] = spca(randn(size(Y{ym})), [], ncomp, delta, -ncomp, maxiter, convergenceCriterion, verbose);
                explainedrand = diag(values)/sum(diag(values));
                nfac_temp = find(explained<explainedrand);
                nfac = nfac_temp(1)-1;
            end
            out.Ycoeff{ym} = spca(Y{ym}, [], nfac, delta, -nfac, maxiter, convergenceCriterion, verbose);
            out.Yscore{ym} = Y{ym}*out.Ycoeff{ym};
            
        end
        

    case 'FA'
        type = S.prep.calc.pred.FA_type;
        
        for ym = 1:length(Y) % for each Y model
           
            % inital num components
            ncomp=size(Y{ym},2); 
            
            % sPCA
            FactorResults = erp_pca(Y{ym}, ncomp, type);
            explained = FactorResults.scree;
            
            % select num components
            if S.prep.calc.pred.PCA_minvar
                nfac_temp = find(explained > S.prep.calc.pred.PCA_minvar*100);
                if isempty(nfac_temp)
                    nfac_temp = ncomp;
                end
                nfac = nfac_exp(1);
            else
                randResults = erp_pca(randn(size(Y{ym})), ncomp, type);
                explainedrand = randResults.scree;
                explainedrand = explainedrand*(sum(explained)/sum(explainedrand)); %scale to same size as real data
                nfac_temp = find(explained<explainedrand);
                if isempty(nfac_temp)
                    nfac=1;
                else
                    nfac = nfac_temp(1)-1;
                end
            end
            FactorResults = erp_pca(Y{ym},nfac,type);
            out.Ycoeff{ym} = FactorResults.FacCof;
            out.Yscore{ym} = FactorResults.FacScr;
            
        end
        

%         %% OLD sPCA code
%         dtab_PCA=D(d).prep.dtab;
% 
%         % PCA parameters
%         delta = inf;
%         maxiter = 3000;
%         convergenceCriterion = 1e-9;
%         verbose = false;
% 
%         % identify range of L1 norm thresholds to test
%         [~,~,A] = svd(X, 'econ');
%         AXX=abs((A'*X')*X);
%         AXX=sort(AXX(:));
%         range = [max(AXX)*2 fliplr(prctile(AXX,[1:100]))];
% 
%         % test all and find optimal solution for sparse PCA
%         out=struct; 
%         nC = intersect(S.prep.calc.pred.PCA_models,1:length(col_names));
%         for K = nC
%             for rn = 1:length(range)
%                 stop=range(rn);
%                 [B,SV,Bsvd,SVsvd] = spca(X, [], K, delta, stop, maxiter, convergenceCriterion, verbose);
%                 scores = X * B;
%                 % test collinearity
%                 rs=corr(scores,'type','Pearson');
%                 rs=abs(triu(rs,1));
%                 coli=find(rs>0.8);
%                 % ensure correct number of PCs is obtained with min variance
%                 % criterion met
%                 if any(B,'all') && isempty(coli) && ~any(isnan(rs),'all') && sum(SV)/sum(SVsvd)>=S.prep.calc.pred.PCA_minvar
%                     out(K).rn=rn;
%                     out(K).B=B;
%                     out(K).SV = SV;
%                     out(K).scores = scores;
%                     break
%                 end
%             end
% 
%             % if criterion is not met, use standard PCA
%             try out(K).B
%             catch
%                 out(K).rn=0;
%                 out(K).B=Bsvd;
%                 out(K).SV = SVsvd;
%                 out(K).scores = X * Bsvd;
%             end
% 
%             % add PCA components to data table
%             Btable = table;
%             Btable.cov = dtab_PCA(:,col_idx).Properties.VariableNames';
%             for k = 1:K
%                 disp(['adding component M' num2str(K) 'PC' num2str(k)])
%                 pcname=['M' num2str(K) 'PC' num2str(k)];
%                 dtab_PCA.(pcname) = out(K).scores(:,k);
%                 Btable.(pcname) = out(K).B(:,k);
%             end
%             D(d).prep.pca{K}=Btable;
%         end
%         D(d).prep.dtab=dtab_PCA;
end

