function D = eegstats_dataprep_covariates(S,varargin)
% inputs: S is a settings structure (see example script)
% outputs: D is the prepared EEG data and predictors as a struct

%% find prepared data
if isempty(varargin)
   % try loading the data if not supplied in varargin
   try
       disp('loading data...');
       load(S.prep.path.inputs,'D') 
   catch
       error('Supply a D structure or path/file of the data')
   end
else
    D= varargin{1};
end

% run though all EEG files in a loop
for d = 1:length(D)
    
    %% Covariates
    if isfield(S.prep.pred,'covariate')
        CVs=S.prep.pred.covariate;
        for ci = 1:length(CVs)
            switch CVs(ci).type
                case 'HGF'
                    pred_out = HGF_covariates(...
                        CVs(ci),...
                        D(d).prep.subname{:},...
                        D(d).prep.dtab.tnums',...
                        S.prep.path.subjects,...
                        0);
                case 'HGF_trialback1'
                    pred_out = HGF_covariates(...
                        CVs(ci),...
                        D(d).prep.subname{:},...
                        D(d).prep.dtab.tnums',...
                        S.prep.path.subjects,...
                        1);
                case 'trials'
                    pred_out = trial_covariate(...
                        CVs(ci),...
                        D(d).prep.subname{:},...
                        D(d).prep.dtab.tnums');
                case 'eegfile'
                    pred_out = eegfile_covariate(...
                        CVs(ci),...
                        D(d).prep.other,...
                        D(d).prep.dtab.tnums');
            end
            repl=find(ismember(pred_out.Properties.VariableNames,D(d).prep.dtab.Properties.VariableNames));
            if ~isempty(repl)
                % add suffixes if labels are replicated
                for rp = 1:length(repl)
                    pred_out.([pred_out.Properties.VariableNames(repl(rp)) 'a']) = pred_out.(pred_out.Properties.VariableNames(repl(rp)));
                    pred_out.(pred_out.Properties.VariableNames(repl(rp)))=[];
                end
            end
            % add pred_out to dtab
            D(d).prep.dtab=[D(d).prep.dtab pred_out];
            % currently unused:
            CVs(ci).level; % subject, trial
        end
    end

    % remove all-NaN predictors
    vt = vartype('numeric');
    numericvar_names=D(d).prep.dtab(:,vt).Properties.VariableNames;
    ip = all(isnan(table2array(D(d).prep.dtab(:,vt))),1);
    D(d).prep.dtab(:,numericvar_names(ip)) = [];
     
    %% continuous covariate predictor transforms
    % Pred data transformation
    if ~isempty(S.prep.calc.pred.transform)
        for pr = 1:length(numericvar_names)
            if strcmp(S.prep.calc.pred.transform,'arcsinh') % need to modify this to apply to only selected predictors
                x=table2array(D(d).prep.dtab(:,numericvar_names));
                D(d).prep.dtab(:,numericvar_names)=array2table(log(x+sqrt(x.^2+1)));
            
            elseif strcmp(S.prep.calc.pred.transform,'rank') % need to modify this to apply to only selected predictors
                x=D(d).prep.dtab(:,numericvar_names);
                [~,D(d).prep.dtab(:,numericvar_names)]=array2table(sort(x));
            end
        end
    end
    
end

if S.prep.output.save
    disp('saving data to disk...')
    save(fullfile(S.prep.path.outputs,[S.prep.sname '.mat']),'D','S','-v7.3')
    disp('...done')
end


%% SUBFUNCTIONS
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

function pred = HGF_covariates(cv,subname,tnums,datfile,trialback)
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
        pred_label={'RT'};
    elseif strcmp(cv.def{ci}{1},'traj')
        % for each traj group... reformat into predictor matrix
        D.HGF = dat;
        S.traj = {cv.def{ci}(2:5)}; 
        [Stemp] = HGF_traj2mat(S,D,cv.def{ci}{6});
        predtemp=Stemp.pred;
        pred_label=Stemp.pred_label;
    end
    
    % optionally use the data from n trials back
    predlabelsuff = '';
    if trialback>0
        predtemp = [nan(trialback,size(predtemp,2)); predtemp(1:end-trialback,:)];
        predlabelsuff = ['_trialback' num2str(trialback)];
    end

    % remove trials not in EEG and put in order of EEG data trials
    for pr = 1:length(pred_label)
        pred.([pred_label{pr} predlabelsuff])=predtemp(tnums,pr);
    end
end

if 0
    % calculate predictive surprise
    u = dat.u(tnums,1);
    pred.HGF_PL_PSurp = -log2(pred.HGF_PL_muhat_1.^u.*(1-pred.HGF_PL_muhat_1).^(1-u));

    % calculate Bayesian surprise
    x = [-3:.1:3];
    for n = 1:length(pred.HGF_PL_muhat_1)
        P = normpdf(x,pred.HGF_PL_muhat_1(n),pred.HGF_PL_sahat_1(n)^(1/2));
        Q = normpdf(x,pred.HGF_PL_mu_1(n),pred.HGF_PL_sa_1(n)^(1/2));
        try
            pred.HGF_PL_BSurp(n)=KLDiv(P,Q);
        catch
            pred.HGF_PL_BSurp(n)=NaN;
        end
    end
    if any(isinf(pred.HGF_PL_BSurp))
        warning('inf values');
        pred.HGF_PL_BSurp(isinf(pred.HGF_PL_BSurp))=NaN;
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
rm_dat = dat<0.000001;
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

function dist=KLDiv(P,Q)
%  dist = KLDiv(P,Q) Kullback-Leibler divergence of two discrete probability
%  distributions
%  P and Q  are automatically normalised to have the sum of one on rows
% have the length of one at each 
% P =  n x nbins
% Q =  1 x nbins or n x nbins(one to one)
% dist = n x 1
if size(P,2)~=size(Q,2)
    error('the number of columns in P and Q should be the same');
end
if sum(~isfinite(P(:))) + sum(~isfinite(Q(:)))
   error('the inputs contain non-finite values!') 
end

%CAB
P(P==0) = realmin;
Q(Q==0) = realmin;
P(P<realmin) = realmin;
Q(Q<realmin) = realmin;

% normalizing the P and Q
if size(Q,1)==1
    Q = Q ./sum(Q);
    P = P ./repmat(sum(P,2),[1 size(P,2)]);
    temp =  P.*log(P./repmat(Q,[size(P,1) 1]));
    temp(isnan(temp))=0;% resolving the case when P(i)==0
    dist = sum(temp,2);
    
    
elseif size(Q,1)==size(P,1)
    
    Q = Q ./repmat(sum(Q,2),[1 size(Q,2)]);
    P = P ./repmat(sum(P,2),[1 size(P,2)]);
    temp =  P.*log(P./Q);
    temp(isnan(temp))=0; % resolving the case when P(i)==0
    dist = sum(temp,2);
end
