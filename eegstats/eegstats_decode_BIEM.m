function D=eegstats_decode_BIEM(S,varargin)
% Bayesian inverse encoding model: requires weights/error from an encoding model

%% find encoding model outcomes
if isempty(varargin)
   % try loading the data if not supplied in varargin
   try
        disp('loading encoding model data')
       load(fullfile(S.biem.path.inputs.encode,'D.mat'),'D') 
   catch
       error('Supply a D structure or path/file of the data')
   end
else
    D= varargin{1};
end

if isempty(S.biem.model.index)
    S.biem.model.index = 1:length(D(1).model);
end

%% load prepped data for predictor values (for prior)
disp('loading prepped data')
Dprep=load(fullfile(S.biem.path.inputs.prep),'D');
Dp=Dprep.D;

% set paths
S.path.code = {
    };
set_paths(S)

if length(D)>1
    save_pref = 'subject_';
else
    save_pref = 'group_';
end

%% Calculate prior covariance
disp('calculating prior covariance for all subjects/models')
for i = S.biem.model.index
    %% Select either single subject or group prior/data, according to S.biem.prior
    if length(D)==1 || strcmp(S.biem.prior,'subject_pred')
        for d = 1:length(D)
            
            % get predictors
            pred = D(d).model(i).fixeddesign;
            trainidx=find(double(Dp(d).prep.dtab.train));
            train_pred = pred(trainidx,:);

            % normalize prior, calculate the covariance of the prior
            % Prior is a multivariate Gaussian with zero mean
            prior = train_pred; % prior: examples x features
            prior_nonan = ~isnan(prior);
            [prior(prior_nonan),prior_mu,prior_sigma] = zscore(double(prior(prior_nonan)));
            priorcov = cov(reshape(prior(prior_nonan),size(prior)));

            % store 
            D(d).model(i).BIEM.pred = pred;
            D(d).model(i).BIEM.prior = {prior, prior_mu, prior_sigma};
            D(d).model(i).BIEM.priorcov = priorcov;
            if 0
                figure; imagesc(R);colormap('hot');
            end
        end
    elseif length(D)>1 && strcmp(S.biem.prior,'group_pred')
        
        % get predictors
        for d = 1:length(D)
            Y1 = Dp(d).prep.Y(1);
            trainidx=find(double(Dp(d).prep.dtab.train));
            if d==1
                pred = D(d).model(i).fixeddesign;
                train_pred = pred(trainidx,:);
            else
                pred = vertcat(pred, D(d).model(i).fixeddesign);
                train_pred = vertcat(train_pred, pred(trainidx,:));
            end
        end

        % normalize prior and calculate covariance. Prior: examples x features
        prior = train_pred; 
        prior_nonan = ~isnan(train_pred);
        [prior(prior_nonan),prior_mu,prior_sigma] = zscore(double(prior(prior_nonan)));
        priorcov = cov(prior(prior_nonan));
        if 0
            figure; imagesc(priorcov);colormap('hot');
        end
        
        % store prior
        for d = 1:length(D)
            D(d).model(i).BIEM.pred = pred;
            D(d).model(i).BIEM.prior = {prior, prior_mu, prior_sigma};
            D(d).model(i).BIEM.priorcov = priorcov;
        end
    elseif strcmp(S.biem.prior,'uniform')
        prior = nan; % unfinished
    end
end
        
%% Loop through subjects/models to obtain data and decode
for d = 1:length(D)
    for i = S.biem.model.index
        disp(['decoding ' save_pref num2str(d) ', model ' num2str(i)])
        
        %% Get the data: keep as Tables for now to select either single subject or whole group later
        % training data and predictors are fed into this function to calculate mean/std of X/Y in
        % training set in order to normalise the test set.
        % Y: brain data (trials x voxels/samples)
        % X: design matrix (trials x features)
        % B: encoding filters (regression weights)
        % Sigma: estimated residual variances used to define the model
        % R: prior covariance matrix (features x features)
        
        % Data for decoding: create array from table
        Y = Dp(d).prep.Y;
        for s = 1:length(Y)
            if s==1
                data = table2array(Y(s).dtab);
            else
                data = horzcat(data,table2array(Y(s).dtab));
            end
        end
        
        % mask the data
        mask = ones(size(data,2),1);
        if ~isempty(S.biem.maskimg)
            maskimg = spm_read_vols(spm_vol(S.biem.maskimg));
            mask = ungrid(maskimg,S)>0;
            mask=mask(:);
            data(:,~mask)=0;
        end
        
        % training predictors
        pred = D(d).model(i).BIEM.pred;
        
        % get fixed and random weights
        if ~isempty(S.biem.randomeffects)
            lmm_rnd = D(d).model(i).randomnames.Level;
            dtab_rnd = cellstr(Dp(d).prep.dtab.(S.biem.randomeffects{:}));
            lmm_rnd_levels = unique(lmm_rnd,'stable');
            dtab_rnd_levels = unique(dtab_rnd,'stable');
            random=zeros(length(Y),length(lmm_rnd_levels),length(D(d).model(i).coeff));
        end
        for b = 1:length(D(d).model(i).coeff)
            % fixed
            fixed(:,b) = D(d).model(i).coeff(b).b(:);
            prednames{b} = D(d).model(i).coeff(b).name;
            if ~isempty(S.biem.randomeffects)
                ridx = strcmp(D(d).model(i).randomnames.Name,prednames{b});
                if any(ridx)
                    random(:,:,b) = reshape(D(d).model(i).random(:,:,ridx),[],sum(ridx));
                end
            end
        end
        
        % sigma: residual variance (e.g. MSE) from encoding. Make sure it's
        % a square matrix with variance on the diagonal.
        sigma = D(d).model(i).s;
        if diff(size(sigma))~=0
           sigma = reshape(sigma,1,[]);
           sigma = diag(sigma);
        end
        
        % prior covariance
        R=D(d).model(i).BIEM.priorcov;
        
        % Decode: ENSURE THE FUNCTION DOES NOT Z-SCORE
        % Add random effects if needed. 
        if ~isempty(S.biem.randomeffects)
            orig=[]; recons=[];
            for rn = 1:length(dtab_rnd_levels)
                lmm_ridx = strcmp(lmm_rnd_levels,dtab_rnd_levels{rn});
                dtab_ridx = strcmp(dtab_rnd,dtab_rnd_levels{rn});
                beta = fixed + squeeze(random(:,lmm_ridx,:));
                trainidx=find(double(Dp(d).prep.dtab.train(dtab_ridx,:)));
                testidx=find(double(Dp(d).prep.dtab.test(dtab_ridx,:)));
                [orig_rnd, recons_rnd] = decode_unimod(data(dtab_ridx,:), pred(dtab_ridx,:), trainidx, testidx, beta', sigma, R);
                orig = [orig;orig_rnd];
                recons = [recons;recons_rnd];
            end
        else
            beta = fixed;
            trainidx=find(double(Dp(d).prep.dtab.train));
            testidx=find(double(Dp(d).prep.dtab.test));
            [orig, recons] = decode_unimod(data, pred, trainidx, testidx, beta', sigma, R);
        end
        
        
        % re-scale
        orig = bsxfun(@plus,bsxfun(@times,orig,D(d).model(i).BIEM.prior{3}),D(d).model(i).BIEM.prior{2});
        recons = bsxfun(@plus,bsxfun(@times,recons,D(d).model(i).BIEM.prior{3}),D(d).model(i).BIEM.prior{2});
        
        % test correlations
        test_corr=[];
        for pn = 1:size(recons,2)
            test_corr(pn) = corr(orig(:,pn),recons(:,pn),'type','Spearman');
            
            for ii = 1:length(recons(:,pn))
                [~,idx]=min(abs(orig(:,pn)-recons(ii,pn)));
                correct(ii,1)=orig(ii,pn)==orig(idx,pn);
            end
            pcorr = round(100*sum(correct)/length(correct));
            
            disp([prednames{pn} ': rho=' num2str(test_corr(pn)),', %corr=' num2str(pcorr)])
        end
        
        % plot orig and reconstruction
        if S.biem.plot_correlations
            fg(d,i)=figure;
            for pn=1:size(recons,2) 
                subplot(1,size(recons,2),pn)
                scatter(orig(:,pn),recons(:,pn));
                title([prednames{pn} ': rho=' num2str(test_corr(pn))]);
            end
        end
        
        % Store results
        D(d).model(i).BIEM.prednames = prednames;
        D(d).model(i).BIEM.test_pred = orig;
        D(d).model(i).BIEM.recons_pred = recons;
        D(d).model(i).BIEM.pred_corr = test_corr;
        
        
    end
end

save(fullfile(S.biem.path.outputs, 'D.mat'),'D');

function set_paths(S)
for p = 1:length(S.path.code)
    if S.path.code{p,1}
        addpath(genpath(S.path.code{p,2}));
    else
        addpath(S.path.code{p,2});
    end
end

function Y = ungrid(YY,S)
% inputs
SPMdata = spm_eeg_load(fullfile(S.biem.path.SPMdata,S.biem.file.SPMdata));
chanind = 1:size(SPMdata,1); % index of chans
ntime = size(YY,3);
n = size(YY,1); % dimension
Cel = spm_eeg_locate_channels(SPMdata, n, chanind);
for ti = 1:ntime
    for ch = chanind
        Y(ch,ti) = YY(Cel(ch,1),Cel(ch,2),ti);
    end
end