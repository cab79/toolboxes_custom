function [D]=eegstats_components_analysis(D,S,varargin)
%- select temporal, spatial, both
%- select PCA/CCA method

% 1st level functions:
% - PCA
% - sPCA
% - factor
% - eigs
% - PLS

% 2nd level functions:
% - CCA

% inputs and outputs:
    % 'D': data structure containing experimental design table dtab, data
    % dimensions, etc.
    % 'grpdata': vertically concatenated EEG data array (over subjects)
    % Inputs are EEG chans/time/trials, outputs are components over trials.
    

% optional Y for PLS
if ~isempty(varargin)
    Y=varargin{1}{1};
    Ynames=varargin{1}{2};
else
    Y={};
    Ydat={};
end

%%
[U,iu,iU]=unique(D.prep.dtab.ID,'stable');
[G,~,iG]=unique(D.prep.dtab.group,'stable');
origgrpdata=D.prep.grpdata{1};
grpdata=origgrpdata;
%D.prep.origgrpdata=origgrpdata;
newdim = D.prep.dim; 

for pt = 1:length(S.type)
    disp(['running PCA: ' S.type{pt}])
    NUM_FAC = S.type_numfac{pt};
    newgrpdata{pt}=[];
    currdim = newdim;

    subdata={};
    dat={};
    switch S.data_in
        case 'conds'
            % Average over trials within each condition, per subject.
            % Currently assumes there are no missing conditions - needs
            % updating to allow missing conditions.
            [C,~,iC]=unique(double(D.prep.dtab.eventTypes));
            nreps = C(end);
            for u=1:length(U)
                for c=1:length(C)
                    dat{u}(c,:,:) = reshape(mean(grpdata(iU==u & iC==c,:),1),1,currdim(1),currdim(2));
                end
                repidx{u} = C;
            end
        case 'trials'
            trialidx = unique(D.prep.dtab.tnums); % repetitions
            nreps = trialidx(end);
            for u=1:length(U)
                dat{u} = reshape(grpdata(iU==u,:),[],currdim(1),currdim(2));
                repidx{u} = D.prep.dtab.tnums(iU==u);
                if ~isempty(Y)
                    for ym = 1:length(Y)
                        Ydat{ym}{u} = Y{ym}(iU==u,:);
                    end
                end
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

    if strcmp(S.type{pt},'spatial')
        obs_dim = 2; % observation dimension
        var_dim = 1; % variable dimension
        % split observations?
        if S.type_obs_sep(pt)
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
                if ~isempty(Y)
                    for ym = 1:length(Ydat)
                        Ydat{ym}{u} = repmat(Ydat{ym}{u},currdim(obs_dim),1);
                    end
                end
            end
        end
    elseif strcmp(S.type{pt},'temporal')
        obs_dim = 1; % observation dimension
        var_dim = 2; % variable dimension
        % split observations?
        if S.type_obs_sep(pt)
            for obs = 1:currdim(1)
                for u=1:length(U)
                    subdata{obs}{u} = reshape(dat{u}(:,obs,:),[],currdim(2));
                end
            end
%                     newdim(obs_dim)=1;
        else
            for u=1:length(U)
                subdata{1}{u} = reshape(dat{u},[],currdim(2)); % trials*chan x time
                if ~isempty(Y)
                    for ym = 1:length(Ydat)
                        Ydat{ym}{u} = repmat(Ydat{ym}{u},currdim(obs_dim),1);
                    end
                end
            end
        end
    elseif strcmp(S.type{pt},'both')
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

    %%% PCA/CCA functions %%%
    O=struct;
    for o = 1:length(subdata)
        if ~obs_dim || length(subdata)>1
            obsdim=1;
        else
            obsdim=newdim(obs_dim);
        end

        % PCA
        [O(o).COEFF,O(o).W,O(o).scores,O(o).mu,O(o).sigma,NUM_FAC,MIN_FAC,varargout{1}] = perform_PCA(S.PCAmethod{1},subdata{o},NUM_FAC,S,rep,nreps*obsdim,Ydat);
        sz=size(O(o).scores);
        O(o).scores_avg = squeeze(mean(reshape(O(o).scores,sz(1),nreps,obsdim,sz(3)),2));
        if ~isempty(varargout)
            if any(strcmp(S.PCAmethod{1},'PCA'))
                O(o).explained = varargout{1}.explained;
                disp(['PCA max var explained after scree: ' num2str(mean(O(o).explained))])

                O(o).Wrand = varargout{1}.Wrand;
                O(o).scoresrand = varargout{1}.PCdata_rand;

                if 0
                    % check recontructed PCA outputs are correlated with
                    % original data
                    for u=1:length(U)
                        % unscaling improves correlation results
                        scores_unscaled = O(o).scores(:,rep{u},u).*repmat(O(o).sigma{u},length(rep{u}),1)';
                        data_pred = scores_unscaled'*pinv(O(o).W(1:NUM_FAC(1),:,u))';
                        % uncentering has no effect on correlation results,
                        % but we do it anyway
                        data_pred_uncentred = data_pred + repmat(O(o).mu{u},length(rep{u}),1);
                        data_test = subdata{o}{u};
                        corr(data_pred_uncentred(:),data_test(:),'type','Spearman')
                        figure; scatter(data_pred_uncentred(:),data_test(:))
                        % highly correlated, rho>0.91
                    end

                end
            elseif any(strcmp(S.PCAmethod{1},'PLS'))
                O(o).PLS = varargout{1};
            elseif any(strcmp(S.PCAmethod{1},'CCA'))
                O(o).CCA = varargout{1};
            end
        end

        % create some randomised data for later
        nrepCCA=S.randomperm;
        scoresrand = nan(size(O(o).scores,1),size(O(o).scores,2),size(O(o).scores,3),nrepCCA);
        for nrp = 1:nrepCCA
            for u=1:length(U)
                for col = 1:size(O(o).scores,1)
                    permidx = rep{u}(randperm(length(rep{u})));
                    scoresrand(col,rep{u},u,nrp) = O(o).scores(col,permidx,u);
                end
            end
        end
        
        try
            NUM_FAC(2) = min(NUM_FAC(2),NUM_FAC(1)); % previously this was inside the next if statement, but caused a crash doign PLS on its own later
        catch
            NUM_FAC(2) = 1;
        end
        if length(S.PCAmethod)>1 && any(strcmp(S.PCAmethod{2},'CCA'))
            disp(['CCA: ' S.type{pt}])
            
            
            % parallel processing
            checkp = gcp('nocreate');
            if isempty(checkp)
                myPool = parpool;
            end
            
            % CCA % https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.23689
            
            % cross-validated r/ncomp: SETTINGS
            reg=1; 
            r= S.cca_reg_weights;
            nr=length(r);
            if S.cca_test_nPCAcomp
                %ncomp =  1 : ceil(size(O(o).scores,1)/min(size(O(o).scores,1),S.cca_test_nPCAcomp)) : size(O(o).scores,1);
                A = 1:1:100; 
                C = unique(round((1+S.cca_test_nPCAcomp/100).^A));
                ncomp = unique([C(C<=size(O(o).scores,1)), size(O(o).scores,1)]);

                % set minimum number
                ncomp = ncomp(ncomp>=MIN_FAC);
            else
                ncomp = size(O(o).scores,1);
                %ncomp = MIN_FAC;
            end
            CV_fold=10;
            CVsample = cvpartition(size(O(o).scores,2),'KFold',CV_fold);
   
            nci=0;
            ncii=[];
            RMSECV=[];
            test_corr = [];
            for nc = 1:length(ncomp) 
                
                % run CV CCAs in parallel for speed
                CV_W = cell(CV_fold,nr);
                CV_W_rand = cell(CV_fold,nr);
                CV_Wn = cell(CV_fold,nr);
                CV_mdata_train = cell(CV_fold,nr);
                CV_mdatan_train = cell(CV_fold,nr);
                CV_mdata_test = cell(CV_fold,nr);
                CV_score_test = cell(CV_fold,nr);
                CV_mdata_train_rand = cell(CV_fold,nr);
                CV_mdatan_train_rand = cell(CV_fold,nr);
                CV_mdata_test_rand = cell(CV_fold,nr);
                CV_score_test_rand = cell(CV_fold,nr);
                nCCA = min([ncomp(nc),NUM_FAC(2)]);
                nCCA_train = nan(CV_fold,nr);
                train_correlations = nan(CV_fold,nr);
                train_correlations_rand = nan(CV_fold,nr);
%                 test_correlations = nan(CV_fold,nr);
                disp(['Estimating CCAs: comp ' num2str(nc) '/' num2str(length(ncomp))])
                parfor cv = 1:CV_fold
                    for ri = 1:nr
                        %disp(['ri ' num2str(ri) ', cv ' num2str(cv)])

                        % train sample
                        CVdata_train = O(o).scores(1:ncomp(nc),training(CVsample,cv),:);
                        [CV_W{cv,ri},CV_mdata_train{cv,ri}, ~, CV_Wn{cv,ri}, CV_mdatan_train{cv,ri}] = mccas(CVdata_train,nCCA,reg,r(ri),O(o).W(1:ncomp(nc),:,:),S);
                        nCCA_train(cv,ri) = size(CV_W{cv,ri},2);

%                         % check outputs: TRAINING
%                         % reconstructed PCs
%                         for u=1:length(U)
%                             PCpred = (squeeze(CV_mdata_train{cv,ri}(1:nCCA_train(cv,ri),ismember(find(training(CVsample,cv)),rep{u}),u))'*pinv(squeeze(CV_W{cv,ri}(:,1:nCCA_train(cv,ri),u))))';
%                             PCtrain = CVdata_train(:,ismember(find(training(CVsample,cv)),rep{u}),u);
%                             cu(u) = corr(PCpred(:),PCtrain(:),'type','Spearman')
%                             % figure; scatter(PCpred(:),PCtest(:))
% 
%                             CV_mdata_train_unscaled = CV_mdata_train{cv,ri}(1:nCCA_train(cv,ri),ismember(find(training(CVsample,cv)),rep{u}),u).*...
%                                 repmat(CV_mdatan_train{cv,ri}(:,u),1,sum(ismember(find(training(CVsample,cv)),rep{u})));
%                             CV_W_unscaled = CV_W{cv,ri}(:,1:nCCA_train(cv,ri),u)*CV_Wn{cv,ri}(:,u);
%                             PCpred_unscaled = (CV_mdata_train_unscaled'*pinv(CV_W_unscaled))';
%                             PCtrain = CVdata_train(:,ismember(find(training(CVsample,cv)),rep{u}),u);
%                             cu_unscaled(u) = corr(PCpred_unscaled(:),PCtrain(:),'type','Spearman')
% 
%                         end
%                         mean(cu)
%                         mean(cu_unscaled)

                        

                        %correlations{cv,ri,nc} = zeros(nCCA_train(cv,ri),2);
                        train_corri = nan(1,nCCA_train(cv,ri));
%                         test_corri = nan(1,nCCA_train(cv,ri));
                        for i = 1:nCCA_train(cv,ri)
                           train_corri(1,i) = mean(mean(corr(squeeze(CV_mdata_train{cv,ri}(i,:,:)),'type','Spearman', 'Rows','pairwise'))); % inter-subject correlations over training data
%                            test_corri(1,i) = mean(mean(corr(squeeze(CV_mdata_test{cv,ri}(i,:,:)),'type','Spearman', 'Rows','pairwise'))); % inter-subject correlations over testing data
                        end
                        train_correlations(cv,ri) = mean(train_corri);
%                         test_correlations(cv,ri) = mean(test_corri);


                        % Test data 
                        CVdata_test = O(o).scores(1:ncomp(nc),test(CVsample,cv),:);
                        if S.cca_select_method == 1

                            % test sample: fix the number of eigenvectors
                            [~,CV_mdata_test{cv,ri},~,~,CV_mdatan_test{cv,ri}] = mccas(CVdata_test,nCCA_train(cv,ri),reg,r(ri),O(o).W(1:ncomp(nc),:,:),S);
                            nCCA_train(cv,ri) = size(CV_mdata_test{cv,ri},1);

                        elseif S.cca_select_method == 2

                            % create CCA scores
                            CV_score_test{cv,ri} = nan(nCCA_train(cv,ri),sum(test(CVsample,cv)),length(U));
                            for u=1:length(U)
                                CV_score_test{cv,ri}(:,ismember(find(test(CVsample,cv)),rep{u}),u) = CV_W{cv,ri}(:,1:nCCA_train(cv,ri),u)' * CVdata_test(:,ismember(find(test(CVsample,cv)),rep{u}),u);
                            end

                            if S.cca_reduce_by_scree
% 
%                                 for nrp = 1:nrepCCA
%                                     for u=1:length(U)
%                                         for col = 1:size(O(o).scores,1)
%                                             scoresrand(col,:,u,nrp);
%                                         end
%                                     end
%                                 end

                                % get training data from randomised PCA
                                % scores
                                CVdata_train_rand = scoresrand(1:ncomp(nc),training(CVsample,cv),:,1);

                                % get test data from randomised PCA
                                % scores
                                CVdata_test_rand = scoresrand(1:ncomp(nc),test(CVsample,cv),:,:);

%                                 CVdata_train_rand = O(o).scoresrand(1:ncomp(nc),training(CVsample,cv),:);
%                                 CVdata_test_rand = O(o).scoresrand(1:ncomp(nc),test(CVsample,cv),:);

                                CV_score_test_rand{cv,ri} = nan(nCCA_train(cv,ri),sum(test(CVsample,cv)),length(U),nrepCCA);
                                
                                if 0
                                    % generate random scores from multiple CCA models
    
                                    % train sample
                                    % THIS NEEDS MODIFICATION: USES WRAND
                                    % FROM PCA FUNCTION CURRENTLY. Need to
                                    % randperm the weights. 
                                    [CV_W_rand{cv,ri},CV_mdata_train_rand{cv,ri}] = mccas(CVdata_train_rand,nCCA,reg,r(ri),O(o).Wrand(1:ncomp(nc),:,:),S);
            
                                    % TRAINING correlations
                                    train_corri_rand = nan(1,nCCA_train(cv,ri));
                                    for i = 1:nCCA_train(cv,ri)
                                       train_corri_rand(1,i) = mean(mean(corr(squeeze(CV_mdata_train_rand{cv,ri}(i,:,:)),'type','Spearman', 'Rows','pairwise'))); % inter-subject correlations over training data
                                    end
                                    train_correlations_rand(cv,ri) = mean(train_corri_rand);
        
                                    for u=1:length(U)
                                        for nrp = 1:nrepCCA
                                            CV_score_test_rand{cv,ri}(:,ismember(find(test(CVsample,cv)),rep{u}),u,nrp) = CV_W_rand{cv,ri}(:,1:nCCA_train(cv,ri),u)' * CVdata_test_rand(:,ismember(find(test(CVsample,cv)),rep{u}),u,nrp);
                                        end
                                    end
                                else
                                    % generate random scores from training
                                    % model (more efficient)
                                    for u=1:length(U)
                                        for nrp = 1:nrepCCA
                                            CV_score_test_rand{cv,ri}(:,ismember(find(test(CVsample,cv)),rep{u}),u,nrp) = CV_W{cv,ri}(:,1:nCCA_train(cv,ri),u)' * CVdata_test_rand(:,ismember(find(test(CVsample,cv)),rep{u}),u,nrp);
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                traincorr_max = max(nanmean(train_correlations,1))
                if S.cca_reduce_by_scree
                    traincorr_max_rand = max(nanmean(train_correlations_rand,1))
                end
%                 testcorr_max = max(nanmean(test_correlations,1))
                
                % how many CCA to test?
                nCCA_train = nCCA_train(:);
                %nCCA_train(nCCA_train==0) = []; % remove zeros
                if S.cca_test_nPCAcomp
                    nCCAtest = 1:min([ncomp(nc),NUM_FAC(2),min(nCCA_train)]);
                elseif S.cca_test_nCCAcomp>0
                    nCCAtest = 1:min([NUM_FAC(2),min(nCCA_train)]);
                else
                    nCCAtest = min([NUM_FAC(2),min(nCCA_train)]);
                end
                
                for ncCCA = nCCAtest % for each CCA component
                    
                    disp(['Reconstructing EEG from PCA comp ' num2str(nc) '/' num2str(length(ncomp)) ', CCA comp ' num2str(ncCCA) '/' num2str(nCCAtest(end))])

                    nci=nci+1;
                    ncii(nci,1)=nc;
                    ncii(nci,2)=ncCCA;

                    % initialise
                    for ri = 1:nr
                        score_test{ri} = nan(ncCCA,size(O(o).scores,2),length(U));
                        if S.cca_reduce_by_scree
                            score_test_rand{ri} = nan(ncCCA,size(O(o).scores,2),length(U),nrepCCA);
                        end
                    end

                    for u=1:length(U)
                        if S.cca_select_method == 1
                            parfor cv = 1:CV_fold
                                for ri = 1:nr
                                    
                                    % reconstructed PCs
                                    PCpred = squeeze(CV_mdata_test{cv,ri}(1:ncCCA,ismember(find(test(CVsample,cv)),rep{u}),u))'*pinv(squeeze(CV_W{cv,ri}(:,1:ncCCA,u)));
    
                                    % test data
                                    testind = test(CVsample,cv);
                                    data_test = subdata{o}{u}(testind(rep{u}),:);
    
                                    % reconstructed EEG data from scaled
                                    % CCA weights and scores.
                                    % No further unscaling from PCA sigma.
    %                                 data_pred = PCpred*pinv(O(o).W(1:ncomp(nc),:,u))';
    
    %                                 % test it
    %                                 corr(data_test(:),data_pred(:),'type','Spearman')
    %                                 sq_diff=bsxfun(@minus,data_test,data_pred).^2;
    %                                 RMSE = squeeze(sqrt(nanmean(sq_diff(:))))
    
                                    % unscale the PCA scores to see if this
                                    % improves it
                                    PCpred_PCunscaled = PCpred.*repmat(O(o).sigma{u}(1,1:ncomp(nc)),sum(ismember(find(test(CVsample,cv)),rep{u})),1);
                                    data_pred_PCunscaled = PCpred_PCunscaled*pinv(O(o).W(1:ncomp(nc),:,u))';
    
    %                                 % test it - improved
    %                                 corr(data_test(:),data_pred_PCunscaled(:),'type','Spearman')
    %                                 sq_diff=bsxfun(@minus,data_test,data_pred_PCunscaled).^2;
    %                                 RMSE = squeeze(sqrt(nanmean(sq_diff(:))))
    
    %                                 % unscale the CCA weights and scores to
    %                                 % generate unscaled PC predictors
    %                                 CV_mdata_test_unscaled = CV_mdata_test{cv,ri}(1:ncCCA,ismember(find(test(CVsample,cv)),rep{u}),u).*...
    %                                     repmat(CV_mdatan_test{cv,ri}(1:ncCCA,u),1,sum(ismember(find(test(CVsample,cv)),rep{u})));
    %                                 CV_W_unscaled = CV_W{cv,ri}(:,1:ncCCA,u)*CV_Wn{cv,ri}(:,u);
    %                                 PCpred_CCAunscaled = (CV_mdata_test_unscaled'*pinv(CV_W_unscaled))';
    % 
    %                                 % test against PCA test data - doesn't
    %                                 % predict it well.
    %                                 CVdata_test = O(o).scores(1:ncomp(nc),test(CVsample,cv),:);
    %                                 PCtest = CVdata_test(:,ismember(find(test(CVsample,cv)),rep{u}),u);
    %                                 corr(PCtest(:),PCpred_CCAunscaled(:),'type','Spearman')
    %                                 %figure; scatter(PCtest(:),PCpred_CCAunscaled(:))
    
    %                                 % test against PCA train data - much better
    %                                 % than for test data. So nothing wrong with
    %                                 % the code.
    %                                 CV_mdata_train_unscaled = CV_mdata_train{cv,ri}(1:ncCCA,ismember(find(training(CVsample,cv)),rep{u}),u).*...
    %                                     repmat(CV_mdatan_train{cv,ri}(1:ncCCA,u),1,sum(ismember(find(training(CVsample,cv)),rep{u})));
    %                                 PCpred_CCAunscaled_train = (CV_mdata_train_unscaled'*pinv(CV_W_unscaled))';
    %                                 CVdata_train = O(o).scores(1:ncomp(nc),training(CVsample,cv),:);
    %                                 PCtrain = CVdata_train(:,ismember(find(training(CVsample,cv)),rep{u}),u);
    %                                 corr(PCtrain(:),PCpred_CCAunscaled_train(:),'type','Spearman')
    %                                 %figure; scatter(PCtrain(:),PCpred_CCAunscaled_train(:))
    
    %                                 % reconstructed EEG data from unscaled
    %                                 % CCA weights and scores AND PC unscaled
    %                                 PCpred_PCunscaled_CCAunscaled = PCpred_CCAunscaled'.*repmat(O(o).sigma{u}(1,1:ncomp(nc)),sum(ismember(find(test(CVsample,cv)),rep{u})),1);
    %                                 data_pred_PCunscaled_CCAunscaled = PCpred_PCunscaled_CCAunscaled*pinv(O(o).W(1:ncomp(nc),:,u))';
    
    %                                 % test it - no improvement in corr - RMSE
    %                                 % is worse
    %                                 corr(data_test(:),data_pred_PCunscaled_CCAunscaled(:),'type','Spearman')
    %                                 sq_diff=bsxfun(@minus,data_test,data_pred_PCunscaled_CCAunscaled).^2;
    %                                 RMSE = squeeze(sqrt(nanmean(sq_diff(:))))
    
                                    % use PC unscaled only
    %                                 RMSECV(ri,cv,nci,u) = 1-corr(data_test(:),data_pred(:),'type','Spearman')^2; % representational similarity
                                    RMSECV(ri,cv,nci,u) = corr(data_test(:),data_pred_PCunscaled(:),'type','Spearman'); % representational similarity
    %                                 RMSECV3(ri,cv,nci,u) = 1-corr(data_test(:),data_pred_PCunscaled_CCAunscaled(:),'type','Spearman')^2; % representational similarity
                                    %PCpred{u} = CV_mdata_test(:,:,u)'*pinv(CV_W(:,:,u));
                                    %sq_diff=bsxfun(@minus,CVdata_test(:,:,u),PCpred{u}').^2;
                                    %RMSECV(ri,cv,u) = squeeze(sqrt(nanmean(sq_diff(:))));
    
                                   
    %                                 % Calculate the total sum of squares (SS_tot)
    %                                 ss_tot = sum((data_test(:)).^2);
    %                                 % Calculate the residual sum of squares (SS_res)
    %                                 ss_res = sum(sq_diff(:));
    %                                 % Calculate the R-squared value
    %                                 r_squared = 1 - (ss_res / ss_tot);
                                end
                            end

                        elseif S.cca_select_method == 2
                            % correlations in test data scores
                            for cv = 1:CV_fold
                                for ri = 1:nr
                                    score_test{ri}(1:ncCCA,intersect(find(test(CVsample,cv)),rep{u}),u) =...
                                        CV_score_test{cv,ri}(1:ncCCA,ismember(find(test(CVsample,cv)),rep{u}),u);
                                    if S.cca_reduce_by_scree
                                        score_test_rand{ri}(1:ncCCA,intersect(find(test(CVsample,cv)),rep{u}),u,:) =...
                                            CV_score_test_rand{cv,ri}(1:ncCCA,ismember(find(test(CVsample,cv)),rep{u}),u,:);
                                    end
                                end
                            end
                        end

                    end

                    if S.cca_select_method == 1
                        meanR = mean(RMSECV(:,:,nci,:),[2 4]);
                        disp(['Max corr: ' num2str(max(meanR))])
                    elseif S.cca_select_method == 2
                        % correlations in test data scores
                        test_corri = nan(nr,ncCCA);
                        for ri = 1:nr
                            for i = 1:ncCCA
                               test_corri(ri,i) = mean(mean(corr(squeeze(score_test{ri}(i,:,:)),'type','Spearman', 'Rows','pairwise'))); % inter-subject correlations over testing data
                            end
                        end
                        %test_corr(:,nci) = nanmean(test_corri,2); % get mean over all CCAs
                        test_corr(:,nci) = test_corri(:,ncCCA); % get last CCA only
                        testcorr_max = max(test_corr(:,nci));
                        disp(['testcorr max: ' num2str(testcorr_max)])

                        if S.cca_reduce_by_scree
                            test_corri_rand = nan(nr,ncCCA);
                            parfor ri = 1:nr
                                for i = 1:ncCCA
                                    for nrp = 1:nrepCCA
                                        test_corri_rand(ri,i,nrp) = mean(mean(corr(squeeze(score_test_rand{ri}(i,:,:,nrp)),'type','Spearman', 'Rows','pairwise'))); % inter-subject correlations over testing data
                                    end
                                end
                            end
                            %test_corr_rand(:,nci) = nanmean(test_corri_rand,2); % get mean over all CCAs
                            for nrp = 1:nrepCCA
                                test_corr_rand(:,nci,nrp) = test_corri_rand(:,ncCCA,nrp); % get last CCA only
                            end
                            [testcorr_max_rand, maxi] = max(mean(test_corr_rand(:,nci,:),3));
                            testcorr_max_rand_std_upper = testcorr_max_rand + 2*std(test_corr_rand(maxi,nci,:));
                            disp(['testcorr rand mean: ' num2str(testcorr_max_rand)])
                            disp(['testcorr rand 2*std upper: ' num2str(testcorr_max_rand_std_upper)])

                        end
                    end

%                     meanR2 = mean(RMSECV2(:,:,nci,:),[2 4]);
%                     disp(['Min error 2: ' num2str(min(meanR2))])
%                     meanR3 = mean(RMSECV3(:,:,nci,:),[2 4]);
%                     disp(['Min error 3: ' num2str(min(meanR3))])
                end
            end
            
            if S.cca_select_method == 1
                % which regularisation parameter?
                RMSECVm = squeeze(mean(RMSECV,[2 4])); % mean over CVs and subjects
                [~,minRi] = max(RMSECVm,[],1); % max reg parameter for every combination of N PCAs and N CCAs
                minR=r(minRi);
                RMSECVr =[];
                for nc = 1:length(minRi)
                    RMSECVr(:,nc,:) = RMSECV(minRi(nc),:,nc,:);
                end
    
%                 % set minimum reg parameter
%                 RMSECVr = RMSECVr(:,minR<1,:);
%                 ncii = ncii(minR<1,:);
%                 minR = minR(minR<1);
                
                % which number of components?
                if isfield(S,'cca_select_nPCA_per_group')
                    select_per_group = S.cca_select_nPCA_per_group;
                else
                    select_per_group = 1;
                end
                figure; hold on
                if select_per_group
                    for g=1:length(G)
                        rmdat = squeeze(mean(RMSECVr(:,:,iG(iu)==g),[1 3]));
                        [mintemp,minYi(g)] = max(rmdat);
                        for nc = 1:length(ncomp)
                           plot(rmdat(ncii(:,1)==nc));
                        end
                        scatter(ncii(minYi(g),2),mintemp);
                    end
                    minYi = max(minYi);
                else
                    rmdat = squeeze(mean(RMSECVr,[1 3]));
                    [mintemp,minYi] = max(rmdat);
                    for nc = 1:length(ncomp)
                       plot(rmdat(ncii(:,1)==nc));
                    end
                    scatter(ncii(minYi,2),mintemp);
                end
                hold off
                O(o).RMSECV = rmdat;
                nCCA = ncii(minYi,2);
                nPCA = ncomp(ncii(minYi,1));

            elseif S.cca_select_method == 2
                % which regularisation parameter?
                [~,minRi] = max(test_corr,[],1); % max reg parameter
                minR=r(minRi);
                testmatr =[];
                testmatrandr =[];
                for nc = 1:length(minRi)
                    testmatr(nc,1) = test_corr(minRi(nc),nc);
                    testmatrandr(nc,:) = test_corr_rand(minRi(nc),nc,:);
                end
    
                figure

                for nc = 1:length(ncomp)
                    hold on
                    plot(testmatr(ncii(:,1)==nc));
                    plot(mean(testmatrandr(ncii(:,1)==nc,:),2),'LineStyle','--');
                    plot(mean(testmatrandr(ncii(:,1)==nc,:),2) + 2*std(testmatrandr(ncii(:,1)==nc,:),[],2),'LineStyle','--');
                    hold off
                end

                % find CCAs that are greater than random instead
                [~,minYi] = max(testmatr); % N PCA with highest value
                testmatrP = testmatr(ncii(:,1)==ncii(minYi,1)); 
                testmatrandrP = testmatrandr(ncii(:,1)==ncii(minYi,1),:); 
                less = find(testmatrP < mean(testmatrandrP,2) + 2*std(testmatrandrP,[],2));

                nCCA = less(1)-1;
                nPCA = ncomp(ncii(minYi,1));

                O(o).testcorr = testmatr;


            end

            disp(['n CCA comp: ' num2str(nCCA) ', n PCA comp: ' num2str(nPCA)])
            disp(['regularisation: ' num2str(minR(minYi))])
            disp(['PCA max var explained after scree: ' num2str(mean(O(o).explained))])
            O(o).CCA_regularisation = minR(minYi);
                        
            % reduce PCs
            O(o).scores = O(o).scores(1:nPCA,:,:);
            O(o).W = O(o).W(1:nPCA,:,:);
            for u=1:length(U)
                O(o).COEFF{u} = O(o).COEFF{u}(:,1:nPCA);
                O(o).sigma{u} = O(o).sigma{u}(:,1:nPCA);
            end
            
            % final CCA solution
            [O(o).W,O(o).mdata] = mccas(O(o).scores,nPCA,reg,minR(minYi),O(o).W,S);
            O(o).W = O(o).W(:,1:nCCA+S.nCCAadd,:);
            O(o).mdata = O(o).mdata(1:nCCA+S.nCCAadd,:,:);
            NUM_FAC(2) = nCCA+S.nCCAadd;
            
            sz=size(O(o).mdata);
            O(o).mdata_avg = squeeze(mean(reshape(O(o).mdata,sz(1),nreps,obsdim,sz(3)),2));
            
        end
    end

%             
    % get PCA and CCA scores for single trials
    for u=1:length(U)
        if strcmp(S.type{pt},'spatial')
            if length(O)==1
                temp = permute(reshape(grpdata(iU==u,:),[],currdim(1),currdim(2)),[1 3 2]);% trials x time x chan
                temp = reshape(temp,[],size(temp,3)); 
            else
                temp = reshape(grpdata(iU==u,:),[],currdim(1),currdim(2));% trials x chan x comp 
            end
        elseif strcmp(S.type{pt},'temporal')
            if length(O)==1
                temp = reshape(grpdata(iU==u,:),[],currdim(1),currdim(2));% trials x chan x time
                temp = reshape(temp,[],size(temp,3)); 
            else
                temp = permute(reshape(grpdata(iU==u,:),[],currdim(1),currdim(2)),[1 3 2]);% trials x time x comp
            end
        elseif strcmp(S.type{pt},'both')
            if length(O)==1
                temp = reshape(grpdata(iU==u,:),[],currdim(1),currdim(2));% trials x chan x time
                temp = reshape(temp,[],size(temp,2)*size(temp,3)); 
            else
                temp = permute(reshape(grpdata(iU==u,:),[],currdim(1),currdim(2)),[1 3 2]);% trials x time x chan
            end
        end

        if var_dim
            newdim(var_dim) = NUM_FAC(2);
        end
        
        
        oind{u}=[];
        PCAs{u} = [];
        for o = 1:length(O)

            if S.centre_output
                tmu = repmat(O(o).mu{u},size(temp,1),1);
                O(o).PCAs{u} = (temp(:,:,o)-tmu)*O(o).COEFF{u};
            else
                O(o).PCAs{u} = temp(:,:,o)*O(o).COEFF{u}; 
            end

            if S.standardise_output && S.centre_output % should onyl standardise if centred
                tsigma = repmat(O(o).sigma{u},size(temp(:,:,o),1),1);
                O(o).PCAs{u} = O(o).PCAs{u}./tsigma;   
            end

%             if S.normalise_output
%                 tsigma = repmat(O(o).sigma{u},size(temp(:,:,o),1),1);
%                 O(o).PCAs{u} = O(o).PCAs{u}./tsigma;   
%             end

            oind{u} = [oind{u}, o*ones(1,size(O(o).PCAs{u},2))];
            PCAs{u} = [PCAs{u}, O(o).PCAs{u}];

        end
        
        if length(S.PCAmethod)>1 && any(strcmp(S.PCAmethod{2},'CCA'))
            CCAs{u} = [];
            for o = 1:length(O)
                CCAs{u} = [CCAs{u}, O(o).PCAs{u}*O(o).W(:,:,u)];
            end
            data_out=CCAs{u};
        else
            data_out=PCAs{u};
        end

        if length(O)>1
            temp = permute(reshape(data_out,length(repidx{u}),NUM_FAC(2),[]),[1 3 2]); % trials x time/chan x cc
        else
            temp = reshape(data_out,length(repidx{u}),[],NUM_FAC(2)); % trials x time/chan x cc
        end

        tempsz = [size(temp,2),size(temp,3)];

        if strcmp(S.type{pt},'spatial')
            temp = reshape(permute(temp,[1 3 2]),[],tempsz(1)*tempsz(2));
        elseif strcmp(S.type{pt},'temporal')
            temp = reshape(temp,[],tempsz(1)*tempsz(2));
        elseif strcmp(S.type{pt},'both')
            temp = squeeze(temp);
        end

        newgrpdata{pt}(iU==u,:) = temp;
    end
    grpdata=newgrpdata{pt};
    if ~obs_dim
        newdim([1 2]) = tempsz;
    end

    % store outputs
    D.prep.PCA(pt).type = S.type{pt};
    D.prep.PCA(pt).O = O;
    D.prep.PCA(pt).oind = oind;
    D.prep.PCA(pt).currdim = currdim;
    D.prep.PCA(pt).options = S;
    D.prep.PCA(pt).PCA = PCAs;
    if length(S.PCAmethod)>1 && any(strcmp(S.PCAmethod{2},'CCA'))
        D.prep.PCA(pt).CCA = CCAs;
    end
    D.prep.grpdata=grpdata;

    % PLOT
    for o = 1:length(O)
        if length(S.PCAmethod)>1 && any(strcmp(S.PCAmethod{2},'CCA'))
            plot_mdata{o}=O(o).mdata;
            plot_mdata_avg{o}=O(o).mdata_avg;
        else
            plot_mdata{o}=O(o).scores;
            plot_mdata_avg{o}=O(o).scores_avg;
        end
    end
    
    cidx=floor(linspace(1,NUM_FAC(2),min(NUM_FAC(2),16)));
    for o = 1:length(O)
        figure('name',[S.type{pt} ', obs ' num2str(o)]); clf;
        cci=0; chandatacorr=[];
        for cc = 1:NUM_FAC(2)
            if strcmp(S.type{pt},'spatial')

                % topo data
                O(o).chandata=[];
                for u=1:length(U)
                    O(o).chandata(cc,:,u) = O(o).COEFF{u}*O(o).W(:,cc,u);
                end

                % PLOTS
                if ismember(cc,cidx)
                    % waveform plot
                    cci=cci+1;
                    mCCA = mean(plot_mdata_avg{o}(cc,:,:),3);
                    if length(mCCA) == D.prep.dim(2)
                        subplot(4,4*2,cci)
                        hold on
                        for u=1:length(U)
                            % plot(CCA{u}(:,cc));
                            plot(plot_mdata_avg{o}(cc,:,u));
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
%                                 s_subdata = reshape(newgrpdata{1}(iU==u,:),[],40,19);
%                                 temp(1,:,u) = corr(plot_mdata{o}(cc,repidx{u},u)',s_subdata(:,:,o)); % 3x60x20 x trials(57)xComp(1083=57x19)
%                             end
%                             chandatacorr(cc,:) = mean(temp,3); % mean over subjects
%                             topo = topotime_3D(chandatacorr(cc,:)',S);
%                             pcolor(topo), shading interp; axis off
%                             title(['comp corr ' num2str(cc)])
                end

            elseif strcmp(S.type{pt},'temporal')

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
                    mCCA = nanmean(plot_mdata_avg{o}(cc,:,:),3)';
                    if length(mCCA) == length(S.img.chansel)
                        topo = topotime_3D(mCCA,S);
                        pcolor(topo), shading interp; axis off
                    end
                    title(['comp ' num2str(cc)])
                end

            end

        end
        D.prep.PCA(pt).O(o).chandatacorr = chandatacorr;
    end


end


% final chan-time data
gavg = mean(origgrpdata,1); % grand average of original EEG data
for o = 1:length(O)
    for cc = 1:NUM_FAC(2)
        temp = [];
        for u=1:length(U)
            temp(1,:,u) = corr(plot_mdata{o}(cc,repidx{u},u)',origgrpdata(iU==u,:),'type','Spearman'); % 3x60x20 x trials(57)xComp(1083=57x19)
        end
        chantimecorrsub(cc,:,:,:,o) = reshape(temp,D.prep.dim(1),D.prep.dim(2),[]); % each subject
        chantimecorr(cc,:,:,o) = squeeze(nanmean(chantimecorrsub(cc,:,:,:,o),4)); % mean over subjects
        chantimecorravg(cc,:,o) = squeeze(mean(temp,3)); 
    end
    
    % correct the sign so that topomaps are consistent with original data
    signcorr = sign(regress(gavg',chantimecorravg'));
    chantimecorrsub = bsxfun(@times,chantimecorrsub,signcorr);
    chantimecorr = bsxfun(@times,chantimecorr,signcorr);
    D.prep.grpdata = bsxfun(@times,grpdata,signcorr');
end

D.prep.PCA(pt).chantimecorrsub = chantimecorrsub;
D.prep.PCA(pt).chantimecorr = chantimecorr;

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

D.prep.dim = newdim;

if ~isempty(Y) && any(strcmp(S.PCAmethod,'PLS'))
    
    % reconstruct Y variables from CCA/PLS/PCA components.
    % currently assumes Y underwent PCA prior to PLS.
    % also assume 'both' option (temporal* spatial).
    
    % model
    ym = D.prep.PCA(1).O(1).PLS.Ymodel;
    % PCA variables (prior to PLS)
    YCOEFF_PCA = D.prep.pred_PCA.Ycoeff{ym}; % (Y x PC) 
    % PLS variables
    YCOEFF = D.prep.PCA(1).O(1).PLS.YCOEFF; % (PC x PLS)
    % CCA variables
    if length(S.PCAmethod)>1 && any(strcmp(S.PCAmethod{2},'CCA'))
        CCA_COEFF = D.prep.PCA(1).O(1).W; % CCA weights/coeffs
    end

    % reconstruct
    for u=1:length(U)
        subdat = grpdata(iU==u,:);
        if length(S.PCAmethod)>1 && any(strcmp(S.PCAmethod{2},'CCA'))
            Yrecon{u} = subdat * CCA_COEFF(:,:,u)' * YCOEFF{u}(:,1:size(CCA_COEFF,1))' * YCOEFF_PCA';
        else
            Yrecon{u} = subdat * YCOEFF{u}' * YCOEFF_PCA';
        end

        % store
        allYrecon(iU==u,:) = Yrecon{u};
    end

    % rescale to original variables over all participants to conserve
    % relative individual differences
    [~,mu,sa] = zscore(D.prep.pred_PCA.Y{ym});
    allYrecon = bsxfun(@plus,bsxfun(@times,zscore(allYrecon),sa),mu);

    for u=1:length(U)
        Yrecon{u} = allYrecon(iU==u,:);
    end

    % store outputs
    D.prep.PCA(1).O(1).PLS.Yrecon = Yrecon;
    D.prep.PCA(1).O(1).PLS.Ynames = Ynames{ym};
    if S.Y_use4output
        D.prep.grpdata=allYrecon;
        D.prep.dim(2) = size(allYrecon,2);
    end

    % cross-correlation heat map of recon Y vs. Y
    f=figure;
    corrmat=corr(Y{ym},allYrecon);
    imagesc(corrmat); % Display correlation matrix as an image
    set(gca, 'XTick', 1:size(allYrecon,2)); % center x-axis ticks on bins
    set(gca, 'YTick', 1:size(Y{ym},2)); % center y-axis ticks on bins
    set(gca, 'XTickLabel', Ynames{ym}); % set x-axis labels
    set(gca, 'YTickLabel', Ynames{ym}); % set y-axis labels
    xlabel('reconstructed Y variables')
    ylabel('predicted Y variables')
    title('Cross-correlation matrix: recon Y vs. Y', 'FontSize', 10); % set title
    colormap('jet'); % Choose jet or any other color scheme
    colorbar 
    caxis([-max(abs(corrmat(:))),max(abs(corrmat(:)))])
        
    
    % cross-correlation heat map of PLS-CCA components vs. Y
    f=figure;
    corrmat=corr(Y{ym},grpdata);
    imagesc(corrmat); % Display correlation matrix as an image
    set(gca, 'XTick', 1:size(grpdata,2)); % center x-axis ticks on bins
    set(gca, 'YTick', 1:size(Y{ym},2)); % center y-axis ticks on bins
%     set(gca, 'XTickLabel', varname); % set x-axis labels
    set(gca, 'YTickLabel', Ynames{ym}); % set y-axis labels
    xlabel('PLS/CCA components')
    ylabel('predicted Y variables')
    title('Cross-correlation matrix: PLS components vs. Y', 'FontSize', 10); % set title
    colormap('jet'); % Choose jet or any other color scheme
    colorbar 
    caxis([-max(abs(corrmat(:))),max(abs(corrmat(:)))])
    
     % cross-correlation between Y components
    f=figure;
    rs = floor(sqrt(length(Y)));
    cs = ceil(length(Y)/rs);
    for ym = 1:length(Y)
        subplot(rs,cs,ym)
        corrmat=corr(Y{ym});
        imagesc(corrmat); % Display correlation matrix as an image
        set(gca, 'XTick', 1:size(Y{ym},2)); % center x-axis ticks on bins
        set(gca, 'YTick', 1:size(Y{ym},2)); % center y-axis ticks on bins
        set(gca, 'XTickLabel', Ynames{ym}); % set x-axis labels
        set(gca, 'YTickLabel', Ynames{ym}); % set y-axis labels
        xlabel('PLS Y variables')
        ylabel('PLS Y variables')
        title('Cross-correlation matrix: PLS Y variables', 'FontSize', 10); % set title
        colormap('jet'); % Choose jet or any other color scheme
        colorbar 
        caxis([-max(abs(corrmat(:))),max(abs(corrmat(:)))])
    end
    
    % cross-correlation heat map of PLS-CCA components
    f=figure;
    corrmat=corr(grpdata);
    imagesc(corrmat); % Display correlation matrix as an image
    set(gca, 'XTick', 1:size(grpdata,2)); % center x-axis ticks on bins
    set(gca, 'YTick', 1:size(grpdata,2)); % center y-axis ticks on bins
%         set(gca, 'XTickLabel', varname); % set x-axis labels
%         set(gca, 'YTickLabel', varname); % set y-axis labels
    xlabel('PLS/CCA components')
    ylabel('PLS/CCA components')
    title('Cross-correlation matrix: PLS components', 'FontSize', 10); % set title
    colormap('jet'); % Choose jet or any other color scheme
    colorbar 
    caxis([-max(abs(corrmat(:))),max(abs(corrmat(:)))])
end

function [COEFF,W,PCdata,mu,sigma,NUM_FAC,MIN_FAC,out_Y] = perform_PCA(PCAmethod,dataAll,NUM_FAC,S,rep,nobs,Y)
% feed in index of 2nd dimension (trials or time)
               
TEMP_FAC=repmat(size(dataAll{1},2),1,2);
if isempty(NUM_FAC)
    NUM_FAC([1,2]) = TEMP_FAC;
else
    NUM_FAC([1,2]) = min([NUM_FAC;TEMP_FAC]);
end

MIN_FAC = 1; % minimum number of factors, default. PCA will estimate this if set.

maxmat = floor((S.maxGig*1e9)^(1/2));
if size(dataAll,2)>maxmat
    error(['downsample the data by ' num2str(size(dataAll,2)/maxmat) ' to enable CCA'])
end

% centre data
%centering the data before PCA ensures that the PCA scores are centered as well, with a mean of zero. This centering step is essential for capturing the covariance structure and interpreting the relative contributions of the principal components.
cdata={};
rand_cdata={};
for i = 1:length(dataAll)
    tmpscore = squeeze(dataAll{i});
    mu{i} = mean(tmpscore);
    cdata{i} = tmpscore - repmat(mu{i},size(tmpscore,1),1);
    rand_cdata{i} = cdata{i};
end

% parallel processing
checkp = gcp('nocreate');
if isempty(checkp)
    myPool = parpool(S.para);
end
 
out_Y=[];

switch PCAmethod

    case 'FA'
        type=S.pca_FA_type;
        % find nfac for each subject
        parfor i = 1:length(cdata)
            disp(['factor analysis on EEG: finding number of factors for subject ' num2str(i)])
            FactorResults = erp_pca(cdata{i},NUM_FAC(1),type);
            randResults = erp_pca(randn(size(cdata{i})),NUM_FAC(1),type);
            randResults.scree = randResults.scree*(sum(FactorResults.scree)/sum(randResults.scree)); %scale to same size as real data
            nfac_temp = find(FactorResults.scree<randResults.scree);
            nfac(i) = min(NUM_FAC(1),nfac_temp(1)-1);
        end
        % take median nfac
        NUM_FAC(1) = floor(median(nfac));
        parfor i = 1:length(cdata)
            disp(['factor analysis on EEG: final solution for subject ' num2str(i)])
            FactorResults = erp_pca(cdata{i},NUM_FAC(1),type);
            COEFF{i} = FactorResults.FacCof;
            W(:,:,i) = COEFF{i}';
            score{i} = FactorResults.FacScr;
        end
    case 'PCA'
       parfor i = 1:length(cdata)
            disp(['PCA on EEG: finding number of factors for subject ' num2str(i)])
            [COEFF{i}, score{i},~,~,explained] = pca(cdata{i},'NumComponents',NUM_FAC(1),'Centered',false,'Algorithm','eig');
            if S.select_ncomp.frac_explained>0
                nfac_exp = find(cumsum(explained) > S.select_ncomp.frac_explained*100);
                if isempty(nfac_exp)
                    nfac_exp = NUM_FAC(1);
                end
                mfac(i) = nfac_exp(1);
            end
            if S.select_ncomp.scree
                %[COEFFrand{i}, scorerand{i},~,~,explainedrand] = pca(randn(size(cdata{i})),'NumComponents',NUM_FAC(1),'Centered',false,'Algorithm','eig');
                
                % The more correct way is to randomly permute the data over
                % trials:
                for col = 1:size(cdata{i},2)
                    rand_cdata{i}(:,col) = cdata{i}(randperm(size(cdata{i},1)),col);
                end
                [COEFFrand{i}, scorerand{i},~,~,explainedrand] = pca(rand_cdata{i},'NumComponents',NUM_FAC(1),'Centered',false,'Algorithm','eig');
                nfac_temp = find(explained<explainedrand);
                nfac(i) = min(NUM_FAC(1),nfac_temp(1)-1);
                explained_by_scree(i) = sum(explained(1:nfac(i)));
            end
       end
        out_Y.explained = explained_by_scree;

        % take median nfac
        MIN_FAC(1) = floor(median(mfac));
        NUM_FAC(1) = floor(median(nfac));
        for i = 1:length(cdata)
            %disp(['PCA on EEG: final solution for subject ' num2str(i)])
            %[COEFF{i}, score{i}] = pca(cdata{i},'NumComponents',NUM_FAC(1),'Centered',false,'Algorithm','eig');
            COEFF{i} = COEFF{i}(:,1:NUM_FAC(1));
            W(:,:,i) = COEFF{i}(:,1:NUM_FAC(1))';
            score{i} = score{i}(:,1:NUM_FAC(1));
            if S.select_ncomp.scree
                COEFFrand{i} = COEFFrand{i}(:,1:NUM_FAC(1));
                Wrand(:,:,i) = COEFFrand{i}(:,1:NUM_FAC(1))';
                scorerand{i} = scorerand{i}(:,1:NUM_FAC(1));
            end
        end
    case 'eigs'
       parfor i = 1:length(cdata)
            [COEFF{i},values]=eigs(cdata{i},NUM_FAC(1));
            explained = diag(values)/sum(diag(values));
            [COEFFrand{i},valuesrand]=eigs(randn(size(cdata{i})),NUM_FAC(1));
            explainedrand = diag(valuesrand)/sum(diag(valuesrand));
            nfac_temp = find(explained<explainedrand);
            nfac(i) = min(NUM_FAC(1),nfac_temp(1)-1);
       end
        % take median nfac
        NUM_FAC(1) = floor(median(nfac));
        parfor i = 1:length(cdata)
            COEFF{i}=eigs(cdata{i},NUM_FAC(1));
            COEFF{i} = COEFF{i}(:,1:NUM_FAC(1));
            W(:,:,i) = COEFF{i}(:,1:NUM_FAC(1))';
            score{i} = cdata{i}*COEFF{i}(:,1:NUM_FAC(1));
        end
    case 'sPCA'
        
        % initial standard PCA to get estimate of number of comps - upper
        % limit on sPCS ncomps
        [~,~,~,~,explained] = pca(cdata{i},'NumComponents',NUM_FAC(1),'Centered',false,'Algorithm','eig');
        [~,~,~,~,explainedrand] = pca(randn(size(cdata{i})),'NumComponents',NUM_FAC(1),'Centered',false,'Algorithm','eig');
        NUM_FAC(1) = find(explained>explainedrand,1,'last');
        
        delta = inf;
        maxiter = 1000;
        convergenceCriterion = 1e-9;
        verbose = false;
        
        parfor i = 1:length(cdata)
            nfac(i)=NUM_FAC(1);
            [COEFF{i},values] = spca(cdata{i}, [], nfac(i), delta, -nfac(i), maxiter, convergenceCriterion, verbose);
            explained = diag(values)/sum(diag(values));
            [COEFFrand{i},valuesrand] = spca(randn(size(cdata{i})), [], nfac(i), delta, -nfac(i), maxiter, convergenceCriterion, verbose);
            explainedrand = diag(valuesrand)/sum(diag(valuesrand));
            nfac_temp = find(explained>explainedrand,1,'last');
            nfac(i) = min(NUM_FAC(1),nfac_temp(1));
        end
        % take median nfac
        NUM_FAC(1) = floor(median(nfac));
        parfor i = 1:length(cdata)
            COEFF{i} = spca(cdata{i}, [], NUM_FAC(1), delta, -NUM_FAC(1), maxiter, convergenceCriterion, verbose);
            COEFF{i} = COEFF{i}(:,1:NUM_FAC(1));
            W(:,:,i) = COEFF{i}(:,1:NUM_FAC(1))';
            score{i} = cdata{i}*COEFF{i}(:,1:NUM_FAC(1));
        end
    case 'PLS'
        for ym = 1:length(Y) % for each Y model
            
            if S.select_ncomp.boot
                % find number of PLS components using bootstrapping
                
                nboot=S.select_ncomp.boot;
                for i = 1:length(cdata)
                    
                    % original sample
                    [~,Los,~,~,Wos] = plsregress(cdata{i},Y{ym}{i},NUM_FAC(1));
                    
                    % bootstrap
                    [n,p]=size(Y{ym}{i});
                    Bsample=floor(unifrnd(0,n,n,nboot))+1;
                    Wb=[];
                    Lb=[];
                    for b=1:nboot
                        index=Bsample(:,b);
                        if mod(b,nboot/10)==0
                            disp(['Model ' num2str(ym) ', subject ' num2str(i) ', Bootstrap replications N. ' num2str(b)]);
                        end
                        [~,Lb(:,:,b),~,~,Wb(:,:,b)] = plsregress(cdata{i}(index,:),Y{ym}{i}(index,:),NUM_FAC(1));
                    end
                    
                    % p values
                    ptype = 'var'; % var: over variables. points: over points.
                    if strcmp(ptype,'var')
                        Wp = mean(Wb,2);
                        Lp = mean(Lb,2);
                        Wop = mean(Wos,2);
                        Lop = mean(Los,2);
                    end
                    stdW=std(Wp,[],3);
                    stdL=std(Lp,[],3);
                    tW{ym}(:,i)=Wop./stdW;
                    tL{ym}(:,i)=Lop./stdL;
                    pW{ym}(:,i)=1-tcdf(abs(tW{ym}(:,i)),n-1);
                    pL{ym}(:,i)=1-tcdf(abs(tL{ym}(:,i)),n-1);
                    grp_nfac{ym}(i) = sum(pL{ym}(:,i)<0.05);
                    
                    if 0
                        % CIs: settings
                        ci_type='basic';
                        conf_level=0.95;
                        bias=false;
                        tail_inf=(1-conf_level)/2*100;
                        tail_sup=100-tail_inf;

    %                     if strcmp(ptype,'var')
    %                         Wb = reshape(Wb,[],1,size(Wb,3));
    %                         Lb = reshape(Lb,[],1,size(Lb,3));
    %                         Los = reshape(Los,[],1,size(Los,3));
    %                         Wos = reshape(Wos,[],1,size(Wos,3));
    %                         if bias==true
    %                             biasW=repmat(mean(Wb,3)-Wos,1,1,2);
    %                             biasL=repmat(mean(Lb,3)-Los,1,1,2);
    %                         else
    %                             biasW=zeros(size(Wb,1),size(Wb,2),1,2);
    %                             biasL=zeros(size(Lb,1),size(Lb,2),1,2);
    %                         end
    %                     else
                            if bias==true
                                biasW=repmat(mean(Wb,3)-Wos,1,2);
                                biasL=repmat(mean(Lb,3)-Los,1,2);
                            else
                                biasW=zeros(size(Wb,1),size(Wb,2),2);
                                biasL=zeros(size(Lb,1),size(Lb,2),2);
                            end
    %                     end

                        switch ci_type
                            case 'basic'
                                CI.W=cat(3,2*Wos-prctile(Wb,tail_sup,3), 2*Wos-prctile(Wb,tail_inf,3))-biasW;
                                CI.L=cat(3,2*Los-prctile(Lb,tail_sup,3), 2*Los-prctile(Lb,tail_inf,3))-biasL;

                            case 'percentile'
                                CI.W=cat(3,prctile(Wb',tail_inf)', prctile(Wb',tail_sup)')-biasW;
                                CI.L=cat(3,prctile(Lb',tail_inf)', prctile(Lb',tail_sup)')-biasL;
                            case 'studentized'
                                tW{ym}=(Wb-repmat(Wos,1,1,nboot))./repmat(stdW,1,1,nboot);
                                tL{ym}=(Lb-repmat(Los,1,1,nboot))./repmat(stdL,1,1,nboot);

                                CI.W=cat(3,Wos-stdW.*prctile(tW{ym}',tail_sup)', Wos-stdW.*prctile(tW{ym}',tail_inf)')-biasW;
                                CI.L=cat(3,Los-stdL.*prctile(tL{ym}',tail_sup)', Los-stdL.*prctile(tL{ym}',tail_inf)')-biasL;

                            otherwise
                                % if ci_type doesn't match with predefined typologies - set 'basic'
                                CI.W=cat(3,2*Wos-prctile(Wb',tail_sup)', 2*Wos-prctile(Wb',tail_inf)')-biasW;
                                CI.L=cat(3,2*Los-prctile(Lb',tail_sup)', 2*Los-prctile(Lb',tail_inf)')-biasL;
                        end

                        %CI.W=cat(3,Wos, mean(Wb,3), Wos-mean(Wb,3), tW, pW, CI.W);
                        %CI.L=cat(3,Los, mean(Lb,3), Los-mean(Lb,3), tL, pL, CI.L);
                        CI.W(CI.W>100)=inf;
                        CI.L(CI.L>100)=inf;
                    end
                    
                end
                % number of factors
                grp_nfac{ym} = max(1,median(grp_nfac{ym}));
                
            elseif S.select_ncomp.MSE_CV || S.select_ncomp.frac_explained
                % find number of PLS components using comparison of CV-MSE
                % between actual and random data; or use fraction explained
                clear explained explainedrand
                for i = 1:length(cdata)
                    [~,~,~,~,BETA{ym}{i},explained{ym}(:,:,i),MSEcv{ym}(:,:,i)] = plsregress(cdata{i},Y{ym}{i},NUM_FAC(1),'cv',5);
                    [~,~,~,~,~,~,MSE{ym}(:,:,i)] = plsregress(cdata{i},Y{ym}{i},NUM_FAC(1));

                    % scale random data
                    si = std(cdata{i}(:));
                    randcdata=randn(size(cdata{i}))*si;
                    [~,~,~,~,BETArand{ym}{i},explainedrand{ym}(:,:,i),MSEcvrand{ym}(:,:,i)] = plsregress(randcdata,Y{ym}{i},NUM_FAC(1),'cv',5);
                    [~,~,~,~,~,~,MSErand{ym}(:,:,i)] = plsregress(randcdata,Y{ym}{i},NUM_FAC(1));

                end
                nfac_msecv = find(mean(MSEcv{ym}(2,:,:),3)<mean(MSEcvrand{ym}(2,:,:),3));

                figure('name',['PLS MSE and explained variance: Model ' num2str(ym)])
                subplot(2,3,1); hold on
                    plot(0:NUM_FAC(1),mean(MSE{ym}(1,:,:),3)'); xlabel('n comp'), ylabel('X CV MSE'); %gcaExpandable;
                    plot(0:NUM_FAC(1),mean(MSErand{ym}(1,:,:),3)','r--'); 
                subplot(2,3,2); hold on
                    plot(0:NUM_FAC(1),mean(MSE{ym}(2,:,:),3)'); xlabel('n comp'), ylabel('Y MSE'); %gcaExpandable;
                    plot(0:NUM_FAC(1),mean(MSErand{ym}(2,:,:),3)','r--'); 
                subplot(2,3,3); hold on
                    plot(0:NUM_FAC(1),mean(MSEcv{ym}(2,:,:),3)'); xlabel('n comp'), ylabel('Y CV MSE'); %gcaExpandable;
                    plot(0:NUM_FAC(1),mean(MSEcvrand{ym}(2,:,:),3)','r--');
                    plot([nfac_msecv(end),nfac_msecv(end)],[0 max(mean(MSEcv{ym}(2,:,:),3))],'k'); 
                subplot(2,3,4); hold on
                    plot(1:NUM_FAC(1),cumsum(mean(explained{ym}(1,:,:),3)')); xlabel('n comp'), ylabel('X explained variance'); %gcaExpandable;
                    plot(1:NUM_FAC(1),cumsum(mean(explainedrand{ym}(1,:,:),3))','r--'); 
                subplot(2,3,5); hold on
                    plot(1:NUM_FAC(1),cumsum(mean(explained{ym}(2,:,:),3)')); xlabel('n comp'), ylabel('Y explained variance'); %gcaExpandable;
                    plot(1:NUM_FAC(1),cumsum(mean(explainedrand{ym}(2,:,:),3))','r--'); 

                nr=ceil(sqrt(length(cdata)));
                figure('name',['PLS fitted responses: Model ' num2str(ym)])
                for i = 1:length(cdata)
                    subplot(nr,nr,i); plot(Y{ym}{i},[ones(size(cdata{i},1),1) cdata{i}]*BETA{ym}{i},'bo');
                        xlabel('Observed Response');
                        ylabel('Fitted Response');
                end
                figure('name',['PLS fitted responses: random: Model ' num2str(ym)])
                for i = 1:length(cdata)
                    subplot(nr,nr,i); plot(Y{ym}{i},[ones(size(cdata{i},1),1) cdata{i}]*BETArand{ym}{i},'bo');
                        xlabel('Observed Response');
                        ylabel('Fitted Response');
                end
                
                if S.select_ncomp.MSE_CV
                    grp_nfac{ym} = nfac_msecv(end);
                elseif S.select_ncomp.frac_explained
                    nfac_exp = find(cumsum(mean(explained{ym}(2,:,:),3)) > S.select_ncomp.frac_explained);
                    if isempty(nfac_exp)
                    	nfac_exp = NUM_FAC(1);
                    end
                    grp_nfac{ym} = nfac_exp(1);
                end
            else
                grp_nfac{ym} = NUM_FAC(1);
            end
            
        end
        
        if 0
            % select model based on best t value
            for ym = 1:length(Y) % for each Y model
                meantL(ym) = mean(abs(tL{ym}(:)));
            end
            [~,minidx] = find(max(meantL));
            
        elseif 0
            % select model based on minimum MSE
            for ym = 1:length(Y) % for each Y model
                MSEcvym(ym) = mean(MSEcv{ym}(2,grp_nfac{ym},:),3)/size(Y{ym}{1},2); 
            end
            [minMSE, minidx] = min(MSEcvym);
            figure; hold on
                    plot(1:length(MSEcvym),MSEcvym); xlabel('model'), ylabel('mean MSE per variable'); %gcaExpandable;
                    scatter(minidx,minMSE,'k'); 
        else
            minidx = 1;
        end
        
        % take grp nfac
        NUM_FAC(1) = grp_nfac{minidx};
        clear explained MSE MSE_CV W
        for i = 1:length(cdata)
            [XCOEFF,YCOEFF,Xscore,Yscore,BETA{i},explained(:,:,i),MSE(:,:,i),stats] = plsregress(cdata{i},Y{minidx}{i},NUM_FAC(1));
            [~,~,~,~,~,~,MSE_CV(:,:,i),~] = plsregress(cdata{i},Y{minidx}{i},NUM_FAC(1),'cv',5);
            
            
            COEFF{i} = XCOEFF;
            W(:,:,i) = COEFF{i}';
            score{i} = Xscore;
            out_Y.YCOEFF{i}=YCOEFF;
            out_Y.Yscore{i}=Yscore;
            out_Y.BETA{i}=BETA{i};
            out_Y.MSE{i}=MSE(:,:,i);
            out_Y.MSE_CV{i}=MSE_CV(:,:,i);
            out_Y.explained{i}=explained(:,:,i);
            out_Y.W{i}=stats.W;
            out_Y.Yfitted{i} = [ones(size(cdata{i},1),1) cdata{i}]*BETA{i};
            
            
            
            if 0
                % manual cross-validation: training R2 and testing Q2
                CV_fold=5;
                CVsample = cvpartition(size(cdata{i},1),'KFold',CV_fold);

                for cv = 1:CV_fold
                    [~,~,~,~,Be] = plsregress(cdata{i}(training(CVsample,cv),:),Y{minidx}{i}(training(CVsample,cv),:),NUM_FAC(1));
                    % training
                    Yp = [ones(size(cdata{i}(training(CVsample,cv),:),1),1) cdata{i}(training(CVsample,cv),:)] * Be; % predicted Y
                    Yt = Y{minidx}{i}(training(CVsample,cv),:); % training Y
                    TSS = sum((Yt-mean(Yt)).^2);
                    RSS = sum((Yp-Yt).^2);
                    R2_train = 1 - RSS/TSS;
                    % testing
                    Yp = [ones(size(cdata{i}(test(CVsample,cv),:),1),1) cdata{i}(test(CVsample,cv),:)] * Be; % predicted Y
                    Yt = Y{minidx}{i}(test(CVsample,cv),:); % test Y
                    TSS = sum((Yt-mean(Yt)).^2);
                    RSS = sum((Yp-Yt).^2);
                    MSE_CV_temp(cv)=RSS/length(Yt);
                    Q2_test = 1 - RSS/TSS;
                end
                MSE_CV2 = mean(MSE_CV_temp);
                out_Y.MSE_CV2{i}=MSE_CV2;
                out_Y.R2_CV_train{i}=R2_train;
                out_Y.Q2_CV_test{i}=Q2_test;
            end
        end
        out_Y.Ymodel=minidx;
%         if exist('YCOEFFPCA','var')
%             out_Y.YCOEFF_PCA=YCOEFFPCA{minidx};
%             out_Y.Yscore_PCA=YscorePCA{minidx};
%         end
        out_Y.grp_nfac = grp_nfac;
        
        figure('name',['PLS MSE and explained variance: Model ' num2str(minidx)])
            subplot(2,3,1);
                boxplot(squeeze(MSE(1,:,:))'); xlabel('n comp'), ylabel('X CV MSE'); 
            subplot(2,3,2); 
                boxplot(squeeze(MSE(2,:,:))'); xlabel('n comp'), ylabel('Y MSE'); %gcaExpandable;
            subplot(2,3,3); 
                boxplot(squeeze(MSE_CV(2,:,:))'); xlabel('n comp'), ylabel('Y CV MSE'); %gcaExpandable;
            subplot(2,3,4); 
                boxplot(squeeze(explained(1,:,:))');  xlabel('n comp'), ylabel('X explained variance'); %gcaExpandable;
            subplot(2,3,5); 
                boxplot(squeeze(explained(2,:,:))'); xlabel('n comp'), ylabel('Y explained variance'); %gcaExpandable;
        
        figure('name',['PLS fitted responses: first ' num2str(NUM_FAC(1)) ' components'])
        nr=ceil(sqrt(length(cdata)));
        for i = 1:length(cdata)
            subplot(nr,nr,i); plot(Y{minidx}{i},out_Y.Yfitted{i},'bo');
                xlabel('Observed Response');
                ylabel('Fitted Response');
        end
        i
    case 'CCA'
        
        ncomp=inf;
        PCA=1;
        if PCA
            disp('CCA per subject: PCA on EEG')
            % PCA of X
            nfac=[];
            parfor i = 1:length(cdata)
                [COEFF{i}, score{i},~,~,explained] = pca(cdata{i},'NumComponents',NUM_FAC(1),'Centered',false,'Algorithm','eig');
                [COEFFrand{i}, scorerand{i},~,~,explainedrand] = pca(randn(size(cdata{i})),'NumComponents',NUM_FAC(1),'Centered',false,'Algorithm','eig');
                nfac_temp = find(explained<explainedrand);
                nfac(i) = min(NUM_FAC(1),nfac_temp(1)-1);
            end
            % take median nfac
            NUM_FAC(1) = floor(median(nfac));
            parfor i = 1:length(cdata)
                [XCOEFF{i}, Xscore{i}] = pca(cdata{i},'NumComponents',NUM_FAC(1),'Centered',false,'Algorithm','eig');
                XCOEFF{i} = XCOEFF{i}(:,1:NUM_FAC(1));
                XW(:,:,i) = XCOEFF{i}(:,1:NUM_FAC(1))';
                Xscore{i} = Xscore{i}(:,1:NUM_FAC(1));
            end
        else
            Yscore = Y;
            Xscore = cdata;
        end
        
        disp('CCA per subject: PCA on Y and CCA')
        for ym = 1:length(Y) % for each Y model
            clear COEFF score COEFFrand scorerand YCOEFF
            
            if PCA
                % PCA of Y
                nfac=[];
                ncomp=size(Y{ym}{i},2);
                for i = 1:length(Y{ym})
                    [COEFF{i}, score{i},~,~,explained] = pca(Y{ym}{i},'NumComponents',ncomp,'Centered',false,'Algorithm','eig');
                    [COEFFrand{i}, scorerand{i},~,~,explainedrand] = pca(randn(size(Y{ym}{i})),'NumComponents',ncomp,'Centered',false,'Algorithm','eig');
                    nfac_temp = find(explained<explainedrand);
                    nfac(i) = nfac_temp(1);
                end
                % take median nfac
                ncomp = floor(median(nfac));
                for i = 1:length(Y{ym})
                    [YCOEFF{i}, Yscore{ym}{i}] = pca(Y{ym}{i},'NumComponents',ncomp,'Centered',false,'Algorithm','eig');
                    YCOEFF{i} = YCOEFF{i};
                    %YW(:,:,i) = YCOEFF{i}';
                end
            end
            
%             if 0
%                 % CCA step: CV on ncomp using Matlab CCA function
%                 for i = 1:length(cdata)
% 
%                     % cross-validate predictions of Yscore and Xscore
%                     CV_fold=5;
%                     CVsample = cvpartition(size(Y{ym}{i},1),'KFold',CV_fold);
%                     for cv = 1:CV_fold
%                         trainidx=training(CVsample,cv);
%                         testidx=test(CVsample,cv);
% 
%                         % mean centre X and Y first?
%                         [Atrain,Btrain,] = canoncorr(Xscore{i}(trainidx,:),Yscore{ym}{i}(trainidx,:));
%                         [Atest,Btest,~,Utest,Vtest,~] = canoncorr(Xscore{i}(testidx,:),Yscore{ym}{i}(testidx,:));
% 
%                         ncomp = size(Atrain,2);
%                         for nc = 1:ncomp
%                             Xpred = Utest(:,1:nc)*pinv(Atrain(:,1:nc));
%                             Ypred = Vtest(:,1:nc)*pinv(Btrain(:,1:nc));
%                             RMSECV{ym}(1,cv,nc)=sqrt(mean(bsxfun(@minus,Xscore{i}(testidx,:),Xpred).^2,'all'));
%                             RMSECV{ym}(2,cv,nc)=sqrt(mean(bsxfun(@minus,Yscore{ym}{i}(testidx,:),Ypred).^2,'all'));
%                         end
%                     end
%                     [minX,minXi] = min(mean(RMSECV{ym}(1,:,:),2));
%                     [minY,minYi] = min(mean(RMSECV{ym}(2,:,:),2));
% 
%                 end
%                 grp_nfac{ym} = minYi;
%                 
%                 figure('name',['CCA MSE: Model ' num2str(ym)])
%                 subplot(2,2,1); hold on
%                     plot(1:ncomp,squeeze(mean(RMSECV{ym}(1,:,:),2))); xlabel('n comp'), ylabel('X CV MSE'); %gcaExpandable;
%     %                 plot([minXi,minXi],[0 max(squeeze(mean(RMSECV{ym}(1,:,:),2)))],'k'); 
%                 subplot(2,2,2); hold on
%                     plot(1:ncomp,squeeze(mean(RMSECV{ym}(2,:,:),2))); xlabel('n comp'), ylabel('Y CV MSE'); %gcaExpandable;
%     %                 plot([grp_nfac{ym},grp_nfac{ym}],[0 max(squeeze(mean(RMSECV{ym}(2,:,:),2)))],'k'); 
%     
%     
%                 % select model based on minimum MSE
%                 for ym = 1:length(Y) % for each Y model
%                     MSEym(ym) = squeeze(mean(RMSECV{ym}(2,:,grp_nfac{minYi}),2))/size(Yscore{ym}{1},2); 
%                 end
%                 [minMSE, minidx] = min(MSEym);
%                 figure; hold on
%                         plot(1:length(MSEym),MSEym); xlabel('model'), ylabel('mean MSE per variable'); %gcaExpandable;
%                         scatter(minidx,minMSE,'k'); 
% 
%                 % take grp nfac
%                 NUM_FAC(1) = grp_nfac{minidx};
%                 for i = 1:length(cdata)
%                     [Atrain,Btrain,~,Ascore,Bscore] = canoncorr(Xscore{i},Yscore{minidx}{i});
%                     COEFF{i} = Atrain(:,1:NUM_FAC(1));
%                     W(:,:,i) = COEFF{i}';
%                     score{i} = Ascore(:,1:NUM_FAC(1));
%                     out_Y.YCOEFF{i}=Btrain(:,1:NUM_FAC(1));
%                     out_Y.Yscore{i}=Bscore(:,1:NUM_FAC(1));
%                     out_Y.MSE{i}=MSEym(minidx);
%                 end
%                 out_Y.grp_nfac = grp_nfac;
%             end
            
            ncomp=min(NUM_FAC(1),ncomp);
            % CCA step: CV on ncomp using Matlab CCA function
            for i = 1:length(cdata)
                
                % cross-validate predictions of Yscore and Xscore
                CV_fold=5;
                CVsample = cvpartition(size(Y{ym}{i},1),'KFold',CV_fold);
                
                reg=1; r= [0:0.05:1];
                for ri = 1:length(r) 
                    for cv = 1:CV_fold
                        trainidx=training(CVsample,cv);
                        testidx=test(CVsample,cv);

                        [~,~, CV_Xweights,CV_Yweights] = dccas(Xscore{i}(trainidx,:)',Yscore{ym}{i}(trainidx,:)',...
                            ncomp,reg,r(ri),S);
                        [CV_mdataX,CV_mdataY,~,~] = dccas(Xscore{i}(testidx,:)',Yscore{ym}{i}(testidx,:)',...
                            ncomp,reg,r(ri),S);

                        for nc = 1:size(CV_mdataX,1) % for each component
                            Xpred = CV_mdataX(1:nc,:)'*pinv(CV_Xweights(:,1:nc));
                            Ypred = CV_mdataY(1:nc,:)'*pinv(CV_Yweights(:,1:nc));
                            RMSECV{ym}(ri,nc,cv,1,i)=sqrt(mean(bsxfun(@minus,Xscore{i}(testidx,:),Xpred).^2,'all'));
                            RMSECV{ym}(ri,nc,cv,2,i)=sqrt(mean(bsxfun(@minus,Yscore{ym}{i}(testidx,:),Ypred).^2,'all'));
                        end
                    end
                end
            end
            % which regularisation parameter?
            [~,minRix{ym}] = min(mean(RMSECV{ym},[2 3 4 5]));
            minRi{ym}=r(minRix{ym});
            % which number of components?
            [~,minYi{ym}] = min(squeeze(mean(RMSECV{ym}(minRix{ym},:,:,2,:),[3 5])));
            grp_nfac{ym} = minYi{ym};

        end
        
        % select model based on minimum MSE
        for ym = 1:length(Y) % for each Y model
            MSEym(ym) = squeeze(mean(RMSECV{ym}(minRix{ym},grp_nfac{ym},:,:,:),[3 4 5]))/size(Yscore{ym}{1},2); 
            MSEym_var(ym) = squeeze(mean(RMSECV{ym}(minRix{ym},grp_nfac{ym},:,:,:),[3 4 5])); 
        end
        [minMSE, minidx] = min(MSEym);
        figure; hold on
                yyaxis left
                plot(1:length(MSEym),MSEym); xlabel('model'), ylabel('mean MSE per variable'); %gcaExpandable;
                scatter(minidx,minMSE,'k'); 
                yyaxis right
                plot(1:length(MSEym_var),MSEym_var); xlabel('model'), ylabel('mean MSE'); %gcaExpandable;
        
        % take grp nfac
        ncomp = grp_nfac{minidx};
        NUM_FAC(1) = ncomp;
        for i = 1:length(cdata)
            [mdataX,mdataY,Xweights,Yweights] = dccas(Xscore{i}',Yscore{minidx}{i}',...
                            ncomp,reg,minRi{minidx},S);
%             [Atrain,Btrain,~,Ascore,Bscore] = canoncorr(Xscore{i},Yscore{minidx}{i});
            COEFF{i} = Xweights(:,1:NUM_FAC(1));
            W(:,:,i) = COEFF{i}';
            score{i} = mdataX';
            out_Y.YCOEFF{i}=Yweights(:,1:NUM_FAC(1));
            out_Y.Yscore{i}=mdataY';
            out_Y.MSE{i}=MSEym(minidx);
        end
        out_Y.grp_nfac = ncomp;
end
% obtain subject-specific PCs
PCdata = nan(NUM_FAC(1),nobs,length(cdata));
sigma = {};
for i = 1:length(cdata)
    sigma{i} = sqrt(var(score{i})); 
    if length(S.PCAmethod)>1 && any(strcmp(S.PCAmethod{2},'CCA')) % standardise scores if they will be used for CCA
       % when the data is centered before performing Principal Component
       % Analysis (PCA), the resulting PCA scores will also be centered, so
       % it is ok to standardise
       score{i} = score{i}./repmat(sqrt(var(score{i})),size(score{i},1),1); % standardise so that each PC from each subject contributes equally to CCA solution
    end
    PCdata(:,rep{i},i) = score{i}';
end

if S.select_ncomp.scree
    PCdata_rand = nan(NUM_FAC(1),nobs);
    sigma_rand = {};
    for i = 1:length(scorerand)
        sigma_rand{i} = sqrt(var(scorerand{i})); 
        if length(S.PCAmethod)>1 && any(strcmp(S.PCAmethod{2},'CCA')) % standardise scores if they will be used for CCA
           % when the data is centered before performing Principal Component
           % Analysis (PCA), the resulting PCA scores will also be centered, so
           % it is ok to standardise
           scorerand{i} = scorerand{i}./repmat(sqrt(var(scorerand{i})),size(scorerand{i},1),1); % standardise so that each PC from each subject contributes equally to CCA solution
        end
        PCdata_rand(:,rep{i},i) = scorerand{i}';
    end
    out_Y.Wrand = Wrand;
    out_Y.PCdata_rand = PCdata_rand;
end
 

function [W, mdata, mWeights, Wn, mdatan] = mccas(data,K,reg,r,weights,Sx)
%% regularized-multiset CCA
%   Date           Programmers               Description of change
%   ====        =================            =====================
%  09/10/2016     Qiong Zhang                 Original code
%% Citation
%  Zhang,Q., Borst,J., Kass, R.E., & Anderson, J.A. (2016) Between-Subject
%  Alignment of MEG Datasets at the Neural Representational Space. 

%% INPUT
% data - input training data to CCA (samples * sensors * subjects), or (samples * PCcomponents * datasets)
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
allsub=[];
for i=1:sub
    allsub=[allsub;data(:,:,i)]; % concatenate components
end
R=cov(allsub.','partialrows'); % cab: added partialrows to allow missing data (nan)
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
    for i = 1:dim2*sub2
        tmp = ceil(i/dim2);
        S2((tmp-1)*dim2+1:tmp*dim2,i) = R2((tmp-1)*dim2+1:tmp*dim2,i); 
    end
    R = R + r*R2;
    S = S + r*S2;
end

if ~isfield(Sx,'cca_method')
    Sx.cca_method = 'eigs';
end

% obtain CCA solutions 
switch Sx.cca_method
    case 'eigs'
        [tempW,values]=eigs((R-S),S,K);
%         if cca_reduce_by_scree
%             
%             allsubrand = randn(size(allsub'));
%             Rrand = cov(allsubrand);
%             Srand = zeros(size(Rrand));
%             for i = 1:dim*sub
%                 tmp = ceil(i/dim);
%                 Srand((tmp-1)*dim+1:tmp*dim,i) = Rrand((tmp-1)*dim+1:tmp*dim,i); 
%             end
%             
%             explained = diag(values)/sum(diag(values));
%             [~,valuesrand]=eigs((Rrand-Srand),Srand,K);
%             explainedrand = diag(valuesrand)/sum(diag(valuesrand));
%             nfac_temp = find(explained<explainedrand);
%             if ~isempty(nfac_temp)
%                 nfac = min(nfac_temp(1)-1, K);
%                 tempW = tempW(:,1:nfac);
%             end
%         end
    case 'FA'
        type = Sx.cca_FA_type;
        FactorResults = erp_pca(allsub',K,type,(R-S),S);
        tempW = FactorResults.FacCof;
%         if cca_reduce_by_scree
%             
%             allsubrand = randn(size(allsub'));
%             Rrand = cov(allsubrand);
%             Srand = zeros(size(Rrand));
%             for i = 1:dim*sub
%                 tmp = ceil(i/dim);
%                 Srand((tmp-1)*dim+1:tmp*dim,i) = Rrand((tmp-1)*dim+1:tmp*dim,i); 
%             end
%             
%             FactorResultsRand = erp_pca(allsubrand,K,type,(Rrand-Srand),Srand);
%             FactorResultsRand.scree = FactorResultsRand.scree*(sum(FactorResults.scree)/sum(FactorResultsRand.scree)); %scale to same size as real data
%             nfac_temp = find(FactorResults.scree<FactorResultsRand.scree);
%             if ~isempty(nfac_temp)
%                 nfac = min(nfac_temp(1)-1, K);
%                 tempW = tempW(:,1:nfac);
%                 
%             end
%         end
end

K = size(tempW,2);
W = zeros(dim,K,sub);
Wn = zeros(1,sub);
for i=1:sub
    if Sx.normalise_cca_weights 
        Wn(1,i) = norm(tempW((i-1)*dim+1:i*dim,:));
        W(:,:,i)=tempW((i-1)*dim+1:i*dim,:)./Wn(1,i);
%         for j = 1:K
%             W(:,j,i)=tempW((i-1)*dim+1:i*dim,j)./norm(tempW((i-1)*dim+1:i*dim,j));
%         end
    else
        W(:,:,i)=tempW((i-1)*dim+1:i*dim,:);
    end
end

% projected data
mdata = nan(K,num,sub);
mdatan = nan(K,sub);
for i = 1:sub
    nandat=isnan(data(:,:,i));
    data_nonan = data(:,:,i);
    if any(nandat,'all')
        data_nonan(nandat)=[];
        data_nonan = reshape(data_nonan,size(data,1),[]);
    end
    mdatatemp = W(:,:,i)'*data_nonan;
    for j = 1:K
        if Sx.normalise_cca_weights
            mdatan(j,i) = norm(mdatatemp(j,:));
            mdatatemp(j,:) = mdatatemp(j,:)./mdatan(j,i);
        else
            mdatatemp(j,:) = mdatatemp(j,:);
        end
    end
    mdata(:,~any(nandat,1),i) = mdatatemp;
end

% sign of the data - does it need changing?

% projection weights
mWeights = 0;
if(reg==1)
    mWeights = zeros(K,num2,sub2);
    for i = 1:sub2
        mWeights(:,:,i) = W(:,:,i)'*weights(:,:,i);
    end
end


function [mdataX, mdataY, Xweights,Yweights] = dccas(X,Y,K,reg,r,Sx)
%% regularized-dual CCA

%% INPUT
% X,Y - input training data to CCA
% K - number of CCAs to retain
% reg - add regularization (1) or not (0)
% r - degree of regularization added
% weights - PCA weight maps for each subject

%% OUTPUT
% W - CCA weights that transform PCAs to CCAs for each subject (PCAs * CCAs * subjects)
% mdata - transformed averaged data (CCAs * sapmles * subjects)
% mWeights - projection weights from sensor to CCAs (CCAs * sensors * subjects)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dimX = size(X,1);
dimY = size(Y,1);
if K>min(dimX,dimY)
    K = min(dimX,dimY);
end
% temp=[];
% for i=1:sub
%     temp=[temp;data(:,:,i)];
% end
XY = [X;Y]; % concatenate components
R=cov(XY','partialrows'); % cab: added partialrows to allow missing data (nan)
S = zeros(size(R));
% for i = 1:dim*sub % subjects*components
%     tmp = ceil(i/dim); % 
%     S((tmp-1)*dim+1:tmp*dim,i) = R((tmp-1)*dim+1:tmp*dim,i); % x: all samples for the subject; y: sample index
% end
S(1:dimX,1:dimX)=R(1:dimX,1:dimX);
S(dimX+1:dimX+dimY,dimX+1:dimX+dimY)=R(dimX+1:dimX+dimY,dimX+1:dimX+dimY);

% add regularization 
if(reg==1)
    R = R + r*eye(length(R));
    S = S + r*eye(length(S));
end

% obtain CCA solutions 
[tempW,values]=eigs((R-S),S,K);

if Sx.normalise_cca_weights 
    Xweights=tempW(1:dimX,:)./norm(tempW(1:dimX,:));
    Yweights=tempW(dimX+1:dimX+dimY,:)./norm(tempW(dimX+1:dimX+dimY,:));
else
    Xweights=tempW(1:dimX,:);
    Yweights=tempW(dimX+1:dimX+dimY,:);
end

% projected data
nandatX=isnan(X);
X_nonan = X;
if any(nandatX,'all')
    X_nonan(nandatX)=[];
    X_nonan = reshape(X_nonan,size(X,1),[]);
end
mdataXtemp = Xweights'*X_nonan;

nandatY=isnan(Y);
Y_nonan = Y;
if any(nandatY,'all')
    Y_nonan(nandatY)=[];
    Y_nonan = reshape(Y_nonan,size(Y,1),[]);
end
mdataYtemp = Yweights'*Y_nonan;
    
for j = 1:K
    if Sx.normalise_cca_weights
        mdataXtemp(j,:) = mdataXtemp(j,:)/norm(mdataXtemp(j,:));
        mdataYtemp(j,:) = mdataYtemp(j,:)/norm(mdataYtemp(j,:));
    end
end
mdataX=nan(K,size(X,2));
mdataY=nan(K,size(Y,2));
mdataX(:,~any(nandatX,1)) = mdataXtemp;
mdataY(:,~any(nandatY,1)) = mdataYtemp;
