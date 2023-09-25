function D = eegstats_dataprep_predictorTransform(S,varargin)
% inputs: S is a settings structure (see example script)
% outputs: D is the prepared EEG data and predictors 

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

%% colormap
cmp=cbrewer('div', 'RdBu', 100, 'pchip');

%% transformations: must take place on grouped data (over subjects) if analysis will be grouped
for d = 1:length(D)
    disp(['data transformations: ' num2str(d) '/' num2str(length(D))])
    
%     if isfield(S.prep.calc.pred,'cubert') && ~isempty(S.prep.calc.pred.cubert) 
%         col_idx=find(contains(D(d).prep.dtab.Properties.VariableNames,S.prep.calc.pred.cubert));
%         for pr = 1:length(col_idx)
%             tempdat=table2array(D(d).prep.dtab(:,col_idx(pr)));
%             if any(tempdat<0)
%                 temp = (abs(tempdat).^(1/3)).*sign(tempdat); 
%             else
%                 temp = (tempdat).^(1/3);
%             end
%             D(d).prep.dtab(:,col_idx(pr)) = array2table(temp);
%         end
%     end
    if isfield(S.prep.calc.pred,'inv') && ~isempty(S.prep.calc.pred.inv) 
        col_idx=find(contains(D(d).prep.dtab.Properties.VariableNames,S.prep.calc.pred.inv));
        for pr = 1:length(col_idx)
            tempdat=table2array(D(d).prep.dtab(:,col_idx(pr)));
            if any(tempdat<0); tempdat = tempdat-min(tempdat); end
            if min(tempdat(tempdat~=0)) ~= max(tempdat(tempdat~=0))
                tempdat(tempdat==0) = min(tempdat(tempdat~=0));
            else
                tempdat(tempdat==0) = realmin;
            end
            temp = 1./tempdat;
            D(d).prep.dtab(:,col_idx(pr)) = array2table(temp);
        end
    end
    if isfield(S.prep.calc.pred,'log') && ~isempty(S.prep.calc.pred.log) 
        col_idx=find(contains(D(d).prep.dtab.Properties.VariableNames,S.prep.calc.pred.log));
        for pr = 1:length(col_idx)
            tempdat=table2array(D(d).prep.dtab(:,col_idx(pr)));
            if any(tempdat<0); tempdat = tempdat-min(tempdat); end
            if min(tempdat(tempdat~=0)) ~= max(tempdat(tempdat~=0))
                tempdat(tempdat==0) = min(tempdat(tempdat~=0));
            else
                tempdat(tempdat==0) = realmin;
            end
            temp = log(tempdat);
            D(d).prep.dtab(:,col_idx(pr)) = array2table(temp);
        end
    end
    if isfield(S.prep.calc.pred,'rootn') && ~isempty(S.prep.calc.pred.rootn)
        for ri = 1:size(S.prep.calc.pred.rootn,1)
            col_idx=find(contains(D(d).prep.dtab.Properties.VariableNames,S.prep.calc.pred.rootn{ri,1}));
            val=S.prep.calc.pred.rootn{ri,2};
            val=val{:};
            for pr = 1:length(col_idx)
                tempdat=table2array(D(d).prep.dtab(:,col_idx(pr)));
                if any(tempdat<0)
                    temp = (abs(tempdat).^(1/val)).*sign(tempdat); 
                else
                    temp = (tempdat).^(1/val);
                end
                D(d).prep.dtab(:,col_idx(pr)) = array2table(temp);
            end
        end
    end
    if isfield(S.prep.calc.pred,'exp') && ~isempty(S.prep.calc.pred.exp) 
        col_idx=find(contains(D(d).prep.dtab.Properties.VariableNames,S.prep.calc.pred.exp));
        for pr = 1:length(col_idx)
            temp = exp(table2array(D(d).prep.dtab(:,col_idx(pr))));
            D(d).prep.dtab(:,col_idx(pr)) = array2table(temp);
        end
    end
    if isfield(S.prep.calc.pred,'power10') && ~isempty(S.prep.calc.pred.power10) 
        col_idx=find(contains(D(d).prep.dtab.Properties.VariableNames,S.prep.calc.pred.power10));
        for pr = 1:length(col_idx)
            temp = 10.^(table2array(D(d).prep.dtab(:,col_idx(pr))));
            D(d).prep.dtab(:,col_idx(pr)) = array2table(temp);
        end
    end
    if isfield(S.prep.calc.pred,'newvarmult') && ~isempty(S.prep.calc.pred.newvarmult) 
        for pr = 1:size(S.prep.calc.pred.newvarmult)
            try
                if strcmp(S.prep.calc.pred.newvarmult{pr,4},'*')
                    temp = D(d).prep.dtab.(S.prep.calc.pred.newvarmult{pr,1}).*D(d).prep.dtab.(S.prep.calc.pred.newvarmult{pr,2});
                elseif strcmp(S.prep.calc.pred.newvarmult{pr,4},'*c') % multiply centered variables
                    % centre per ppt
                    ID = unique(D(d).prep.dtab.ID);
                    temp = nan(height(D(d).prep.dtab),1);
                    for id = 1:length(ID)
                        idx = ismember(D(d).prep.dtab.ID,ID(id));
                        temp(idx) = D(d).prep.dtab.(S.prep.calc.pred.newvarmult{pr,1})(idx)-nanmean(D(d).prep.dtab.(S.prep.calc.pred.newvarmult{pr,1})(idx))...
                            .* D(d).prep.dtab.(S.prep.calc.pred.newvarmult{pr,2})(idx)-nanmean(D(d).prep.dtab.(S.prep.calc.pred.newvarmult{pr,2})(idx));
                    end
                elseif strcmp(S.prep.calc.pred.newvarmult{pr,4},'/')
                    if isnumeric(S.prep.calc.pred.newvarmult{pr,1})
                        % this is the inverse transform
                        vals = D(d).prep.dtab.(S.prep.calc.pred.newvarmult{pr,2});
                        negvalidx = vals<0;
                        temp = S.prep.calc.pred.newvarmult{pr,1}./abs(vals);
                        temp(negvalidx) = -temp(negvalidx);
                    else
                        temp = D(d).prep.dtab.(S.prep.calc.pred.newvarmult{pr,1})./D(d).prep.dtab.(S.prep.calc.pred.newvarmult{pr,2});
                    end
                elseif strcmp(S.prep.calc.pred.newvarmult{pr,4},'-')
                    temp = D(d).prep.dtab.(S.prep.calc.pred.newvarmult{pr,1})-D(d).prep.dtab.(S.prep.calc.pred.newvarmult{pr,2});
%                 elseif strcmp(S.prep.calc.pred.newvarmult{pr,4},'k-') % previous trial
%                     tidx = 1:length(D(d).prep.dtab.(S.prep.calc.pred.newvarmult{pr,2})) S.prep.calc.pred.newvarmult{pr,1};
%                     temp = D(d).prep.dtab.(S.prep.calc.pred.newvarmult{pr,2})(tidx);
                elseif strcmp(S.prep.calc.pred.newvarmult{pr,4},'+')
                    temp = D(d).prep.dtab.(S.prep.calc.pred.newvarmult{pr,1})+D(d).prep.dtab.(S.prep.calc.pred.newvarmult{pr,2});
                end
                D(d).prep.dtab.(S.prep.calc.pred.newvarmult{pr,3}) = temp;
            end
        end
    end
    if isfield(S.prep.calc.pred,'log_2nd') && ~isempty(S.prep.calc.pred.log_2nd) 
        col_idx=find(contains(D(d).prep.dtab.Properties.VariableNames,S.prep.calc.pred.log_2nd));
        for pr = 1:length(col_idx)
            temp = log(table2array(D(d).prep.dtab(:,col_idx(pr))));
            D(d).prep.dtab(:,col_idx(pr)) = array2table(temp);
        end
    end
    
    if isfield(S.prep.calc.pred,'rootn_2nd') && ~isempty(S.prep.calc.pred.rootn_2nd)
        for ri = 1:size(S.prep.calc.pred.rootn_2nd,1)
            col_idx=find(contains(D(d).prep.dtab.Properties.VariableNames,S.prep.calc.pred.rootn_2nd{ri,1}));
            val=S.prep.calc.pred.rootn_2nd{ri,2};
            val=val{:};
            for pr = 1:length(col_idx)
                tempdat=table2array(D(d).prep.dtab(:,col_idx(pr)));
                if any(tempdat<0)
                    temp = (abs(tempdat).^(1/val)).*sign(tempdat); 
                else
                    temp = (tempdat).^(1/val);
                end
                D(d).prep.dtab(:,col_idx(pr)) = array2table(temp);
            end
        end
    end
    % z-scoring on continuous predictors, but not if doing PLS on EEG data
    if S.prep.calc.pred.zscore
        %if S.prep.calc.eeg.pca.on == 1 && any(strcmp(S.prep.calc.eeg.pca.PCAmethod,'PLS'))
        %    disp('no z-scoring to properly rescale Y variables for PLS')
        %else
        if S.prep.calc.pred.zscore_perID
            [U,~,iU]=unique(D.prep.dtab.ID,'stable');
        else
            U={'all'}; iU=ones(height(D.prep.dtab),1);
        end

        for u=1:length(U)
            vt = vartype('numeric');
            numericvar_names=D(d).prep.dtab(:,vt).Properties.VariableNames;
            for pr = 1:length(numericvar_names)
                if ismember(numericvar_names{pr}, S.prep.calc.pred.zscore_exclude); continue; end
                % remove inf values
                dat=table2array(D(d).prep.dtab(iU==u,numericvar_names{pr}));
                if any(isinf(dat))
                    dat(isinf(dat))= max(dat(~isinf(dat)));
                end
                [zdata, D(d).prep.pred_means(pr,u), D(d).prep.pred_stds(pr,u)] = zscore(dat);
                D(d).prep.dtab(iU==u,numericvar_names{pr}) = array2table(zdata);
            end
            
        end
        
        % additionally, convert any categorical variables that are due
        % for PCA to numeric and z-score them
        if S.prep.calc.pred.PCA_on
            ct = vartype('categorical');
            catvar_names=D(d).prep.dtab(:,ct).Properties.VariableNames;
            for a2t = S.prep.calc.pred.PCA_model_add2dtab
                catvar_pca_idx = find(endsWith(catvar_names,S.prep.calc.pred.PCA_cov{a2t}));
                if numel(catvar_pca_idx)>0
                    for pr = 1:length(catvar_pca_idx)
                        if ismember(catvar_names{catvar_pca_idx(pr)}, S.prep.calc.pred.zscore_exclude); continue; end
                        [zdata, D(d).prep.pred_means(pr,1), D(d).prep.pred_stds(pr,1)] = zscore(double(table2array(D(d).prep.dtab(:,catvar_names{catvar_pca_idx(pr)}))));
                        D(d).prep.dtab(:,catvar_names{catvar_pca_idx(pr)}) = [];
                        D(d).prep.dtab(:,catvar_names{catvar_pca_idx(pr)}) = array2table(zdata);
                    end
                end
            end
        end    
            
        %end
    %elseif S.prep.calc.pred.PCA_on == 1
    %    error('must z-score data prior to PCA')
    end
    
    
    
    if S.prep.calc.pred.PCA_on

        if S.prep.calc.pred.PCA_LMM
            [corrmat,pmat,var_names] = LMM_correlation_matrix(D(d),S,0);
            S.LMM_corrmat = corrmat;
            S.LMM_varnames = var_names;
        end
        D(d).prep.pred_PCA = predictor_PCA(D(d).prep.dtab,S);
        if isempty(D(d).prep.pred_PCA.Y)
            return;
        end
        
        if strcmp(S.prep.calc.pred.PCA_type,'FA') || (strcmp(S.prep.calc.pred.PCA_type,'PLS') && S.prep.calc.pred.pls.rotate_factors)
            pattern = D(d).prep.pred_PCA.Ypattern{1};
            pattern_thresh = D(d).prep.pred_PCA.Ypattern{1}.*double(D(d).prep.pred_PCA.Ypattern{1}>=S.prep.calc.pred.PCA_model_minloading | D(d).prep.pred_PCA.Ypattern{1}<=-S.prep.calc.pred.PCA_model_minloading);
            structure = D(d).prep.pred_PCA.Ystruct{1};
            plot_names = D.prep.pred_PCA.col_names{1};
            fig=figure;
            ax=subplot(1,3,1);
            imagesc(pattern);colorbar;colormap(cmp); caxis([-1 1])
            xticks(1:size(pattern,2)); 
            yticks(1:length(plot_names)); yticklabels(plot_names);
            title('pattern matrix (unique var - regress coeffs)')
            set(ax, 'XAxisLocation', 'top');
            ax=subplot(1,3,2);
            imagesc(pattern_thresh);colorbar;colormap(cmp); caxis([-1 1])
            xticks(1:size(pattern,2)); 
            yticks(1:length(plot_names)); yticklabels(plot_names);
            title('pattern matrix, thresholded')
            set(ax, 'XAxisLocation', 'top');
            ax=subplot(1,3,3);
            imagesc(structure);colorbar;colormap(cmp); caxis([-1 1])
            xticks(1:size(pattern,2)); 
            yticks(1:length(plot_names)); yticklabels(plot_names);
            title('structure matrix (corr between vars and facs)')
            set(ax, 'XAxisLocation', 'top');

            idx = eye(size(D.prep.pred_PCA.Ycorr{1}));
            avg_btw_fac_corr = mean(abs(D.prep.pred_PCA.Ycorr{1}(~idx)));
            avg_wth_fac_corr = mean(abs(structure),[1 2]);
            prop_wth_btw = (avg_wth_fac_corr/(avg_btw_fac_corr+avg_wth_fac_corr))*100;

            vardata=table2array(D(d).prep.dtab(:,endsWith(D(d).prep.dtab.Properties.VariableNames,S.prep.calc.pred.PCA_cov{1})));
            varnames=D(d).prep.dtab.Properties.VariableNames(endsWith(D(d).prep.dtab.Properties.VariableNames,S.prep.calc.pred.PCA_cov{1}));
            
            % Cronbach
            ym=S.prep.calc.pred.PCA_model_add2dtab;
            for nf = 1:size(D(d).prep.pred_PCA.Yscore{ym},2)
                if sum(pattern_thresh(:,nf)~=0)>1
                    crona(nf)=cronbach(vardata(:,pattern_thresh(:,nf)~=0));
                else
                    crona(nf)=NaN;
                end
            end

            dim = [0 .7 .3 .3];
            str = {
                ['TotVar = ' num2str(D.prep.pred_PCA.Yfacvartotal{1})],...
                ['FacVar = ' num2str(D.prep.pred_PCA.Yfacvar{1})],...
                ['FacVarUniq = ' num2str(D.prep.pred_PCA.Yfacvaruniq{1})],...
                ['prop_wth_var = ' num2str(prop_wth_btw) '%'],...
                ['cronbach = ' num2str(crona)]
                };
            annotation('textbox',dim,'String',str,'FitBoxToText','on');

            savefig(fig,fullfile(S.prep.path.outputs,[S.prep.sname 'FApattern.fig']));
        end
        
        if strcmp(S.prep.calc.pred.PCA_type,'FA') || (strcmp(S.prep.calc.pred.PCA_type,'PLS') && S.prep.calc.pred.pls.rotate_factors)
           
            % add FA components to data table
            if S.prep.calc.pred.PCA_model_add2dtab
                for k = 1:size(D(d).prep.pred_PCA.Yscore{ym},2)
                    if any(pattern_thresh(:,k)>0)
                        pcname=[S.prep.calc.pred.PCA_type num2str(k)];
                        %D(d).prep.dtab.(pcname) = zscore((transpose(pattern_thresh(:,k))*transpose(vardata))'); % this is incorrect
                        D(d).prep.dtab.(pcname) = zscore(D(d).prep.pred_PCA.Yscore{ym}(:,k));
                    end
                end

                % trialback data
    %             if ~isempty(D(d).prep.pred_PCA.Yscore_tb)
    %                 for k = 1:size(D(d).prep.pred_PCA.Yscore_tb{ym},2)
    %                     pcname=[S.prep.calc.pred.PCA_type '_trialback_' num2str(k)];
    %                     D(d).prep.dtab.(pcname) = D(d).prep.pred_PCA.Yscore_tb{ym}(:,k);
    %                 end
    %             end
    
                % remove factor variance from original predictors
                if S.prep.calc.pred.PCA_model_rmFacVar && S.prep.calc.pred.PCA_model_rmPredVar == 0
                    for nf = 1:size(D(d).prep.pred_PCA.Yscore{ym},2)
                        vnme = varnames(pattern_thresh(:,nf)~=0);
                        pcname=[S.prep.calc.pred.PCA_type num2str(nf)];
                        for v = 1:length(vnme)
                            formula = [vnme{v} '~' pcname '+(1|ID)'];
                            lme = fitlme(D(d).prep.dtab,formula);
                            D(d).prep.dtab.(vnme{v}) = residuals(lme);
                        end
                    end
                end
                
                % remove other predictor variance
                if S.prep.calc.pred.PCA_model_rmPredVar == 1 % grouped by factor
                    dtab = D(d).prep.dtab;
                    for nf = 1:size(D(d).prep.pred_PCA.Yscore{ym},2)
                        vnme = varnames(pattern_thresh(:,nf)~=0);
                        pcname=[S.prep.calc.pred.PCA_type num2str(nf)];
                        for v = 1:length(vnme)
                            formula = [vnme{v} '~'];
                            for v2 = 1:length(vnme)
                                if v2==v; continue; end
                                formula = [formula vnme{v2} '+'];
                            end
                            if S.prep.calc.pred.PCA_model_rmFacVar
                                formula = [formula pcname '+'];
                            end
                            formula = [formula '(1|ID)'];
                            lme = fitlme(D(d).prep.dtab,formula);
                            dtab.(vnme{v}) = residuals(lme);
                        end
                    end
                    D(d).prep.dtab = dtab;
                elseif S.prep.calc.pred.PCA_model_rmPredVar == 2 % all predictors
                    dtab = D(d).prep.dtab;
                    vnme = varnames;
                    for v = 1:length(vnme)
                        formula = [vnme{v} '~'];
                        for v2 = 1:length(vnme)
                            if v2==v; continue; end
                            formula = [formula vnme{v2} '+'];
                        end
                        if S.prep.calc.pred.PCA_model_rmFacVar
                            for nf = 1:size(D(d).prep.pred_PCA.Yscore{ym},2)
                                pcname=[S.prep.calc.pred.PCA_type num2str(nf)];
                                formula = [formula pcname '+'];
                            end
                        end
                        formula = [formula '(1|ID)'];
                        lme = fitlme(D(d).prep.dtab,formula);
                        dtab.(vnme{v}) = residuals(lme);
                    end
                    D(d).prep.dtab = dtab;
                end
                
                % zscore again
                vnme = varnames;
                for v = 1:length(vnme)
                    dat=table2array(D(d).prep.dtab(:,vnme{v}));
                    [zdata, D(d).prep.pred_means(v), D(d).prep.pred_stds(v)] = zscore(dat);
                    D(d).prep.dtab(:,vnme{v}) = array2table(zdata);
                end
            end
            
            
        else
            
            % add PCA components to data table
            if S.prep.calc.pred.PCA_model_add2dtab
                for ym=S.prep.calc.pred.PCA_model_add2dtab
                    for k = 1:size(D(d).prep.pred_PCA.Yscore{ym},2)
                        pcname=[S.prep.calc.pred.PCA_type num2str(ym) '_' num2str(k)];
                        D(d).prep.dtab.(pcname) = zscore(D(d).prep.pred_PCA.Yscore{ym}(:,k));
                    end
                end
            end
        end
    end
    
    if S.prep.calc.pred.PCA_fit
        % Obtain good of fit of whole model and of each factor excluded,
        % for each participant
        
        PCA_fit=table;
        
        [U,iu,iU]=unique(D(d).prep.dtab.ID,'stable');
        for u=1:length(U)
        
            % data covariance
            datacov = tril(cov(D(d).prep.pred_PCA.Y{1}(iU==u,:)),-1);

            % model covariance (whole model)
            modelcov = tril(cov(D(d).prep.pred_PCA.Yscore{1}(iU==u,:)*D(d).prep.pred_PCA.Ycoeff{1}'),-1);

            % goodness of fit (whole model)
            PCA_fit.whole(u) = goodnessOfFit({modelcov(modelcov~=0)},{datacov(datacov~=0)},'NRMSE');

            % for each factor excluded:
            fac_all = 1:size(D(d).prep.pred_PCA.Yscore{1},2);
            for k = fac_all
                fac_include = setdiff(fac_all,k);

                % model covariance (factor excluded)
                modelcov = tril(cov(D(d).prep.pred_PCA.Yscore{1}(iU==u,fac_include)*D(d).prep.pred_PCA.Ycoeff{1}(:,fac_include)'),-1);

                % goodness of fit (factor excluded)
                if length(fac_all)>1
                    PCA_fit.(['noFA' num2str(k)])(u) = goodnessOfFit({modelcov(modelcov~=0)},{datacov(datacov~=0)},'NRMSE');
                end
            end
        end
        
        D(d).prep.pred_PCA.PCA_fit = PCA_fit;
        
        figure;
        means = mean(table2array(PCA_fit),1);
        stds = std(table2array(PCA_fit),[],1);
        if length(fac_all)>1
            errorbar([0 fac_all],means,stds)
            xticks([0 fac_all]);
            xticklabels(PCA_fit.Properties.VariableNames)
            xlim([-0.5 fac_all(end)+0.5])
        else
            errorbar(0,means,stds)
            xticks([0]);
            xticklabels(PCA_fit.Properties.VariableNames)
            xlim([-0.5 +0.5])
        end
        title('NRMSE of PCA/FA models')
        
    end
    
    if S.prep.calc.pred.factor_regressions
        
        % variables
        % var_names = D(d).prep.dtab.Properties.VariableNames;
        var_names = {};
        for pv = 1:length(S.prep.calc.pred.plot)
            Index = regexp(D(d).prep.dtab.Properties.VariableNames, regexptranslate('wildcard', S.prep.calc.pred.plot{pv}));
            Index = find(not(cellfun('isempty',Index)));
            var_names = [var_names D(d).prep.dtab.Properties.VariableNames(Index)];
        end
        if ~isempty(S.prep.calc.pred.plotx)
            xvar_names = {};
            for pv = 1:length(S.prep.calc.pred.plotx)
                Index = regexp(var_names, regexptranslate('wildcard', S.prep.calc.pred.plotx{pv}));
                Index = find(not(cellfun('isempty',Index)));
                xvar_names = [xvar_names var_names(Index)];
            end
        else
            xvar_names = var_names;
        end
        xvar_index = ismember(var_names,xvar_names);
        %var_names = setdiff(var_names,{'ID','group','test','train','eventTypes'},'stable');
        vt = vartype('numeric');
        numericvar_names=D(d).prep.dtab(:,vt).Properties.VariableNames;
        numericvar_ind = find(ismember(var_names,numericvar_names));
        nonnumericvar_ind = find(~ismember(var_names,numericvar_names));
        all_ind = [numericvar_ind,nonnumericvar_ind]; % must be in this order
        
        % combinations
        comb = nchoosek(all_ind,2);   
        comb = comb(ismember(comb(:,1),numericvar_ind),:); 
        corrmat = diag(zeros(1,length(var_names)));
        pmat = diag(zeros(1,length(var_names)));
        for ci = 1:length(comb)
            formula = [var_names{comb(ci,1)} '~' var_names{comb(ci,2)} '+(1|ID)'];
            lme = fitlme(D(d).prep.dtab,formula);
            if strcmp(S.prep.calc.pred.output_metric,'r')
                corrmat(comb(ci,1),comb(ci,2)) = lme.Coefficients(2,2);
            elseif strcmp(S.prep.calc.pred.output_metric,'r2')
                corrmat(comb(ci,1),comb(ci,2)) = sign(double(lme.Coefficients(2,2)))*lme.Rsquared.Ordinary;
            end
            pmat(comb(ci,1),comb(ci,2)) = lme.Coefficients(2,6);
        end
        corrmat = tril(corrmat + corrmat',-1);
        pmat = tril(pmat + pmat',-1);
        corrmat = corrmat.*double(pmat<S.prep.calc.pred.sig);
        if strcmp(S.prep.calc.pred.output_metric,'r')
            corrmat = corrmat.*double(corrmat>=S.prep.calc.pred.R_min | corrmat<-S.prep.calc.pred.R_min);
        elseif strcmp(S.prep.calc.pred.output_metric,'r2')
            corrmat = corrmat.*double(corrmat>=(S.prep.calc.pred.R_min)^2 | corrmat<-(S.prep.calc.pred.R_min)^2);
        end
        figure;
        ax=subplot(1,1,1);
        imagesc(corrmat(:,xvar_index));colorbar;colormap(cmp); caxis([-1 1])
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
        xticks(1:sum(xvar_index)); xticklabels(plot_names(xvar_index)); %xtickangle(90);
        yticks(1:length(plot_names)); yticklabels(plot_names);
        title('factor regressions (concat)')
        set(ax, 'XAxisLocation', 'top');
        
%         save(fullfile(S.prep.path.outputs,'predictor_correlations.mat'),'corrmat','var_names');
%         coli=find(corrmat>S.prep.calc.pred.test_collinearity);
%         [row,col] = ind2sub(size(corrmat),coli);
%         for ci = 1:length(coli)
%             var1=var_names{row(ci)};
%             var2=var_names{col(ci)};
%             disp(['collinear predictors: ' var1 ', ' var2])
%         end
%         if S.prep.calc.pred.allow_collinearity==false && ~isempty(coli)
%             error('collinear predictors - see above')
%         end
        
    end
    
    
    if S.prep.calc.pred.test_collinearity
        
        [corrmat,pmat,var_names] = LMM_correlation_matrix(D,S,S.prep.calc.pred.test_collinearity_PCAcov_only);

        if length(var_names)<2
            return
        end

        corrmat = corrmat.*double(pmat<S.prep.calc.pred.sig);
        if strcmp(S.prep.calc.pred.output_metric,'r')
            corrmat = corrmat.*double(corrmat>=S.prep.calc.pred.R_min | corrmat<-S.prep.calc.pred.R_min);
        elseif strcmp(S.prep.calc.pred.output_metric,'r2')
            corrmat = corrmat.*double(corrmat>=(S.prep.calc.pred.R_min)^2 | corrmat<-(S.prep.calc.pred.R_min)^2);
        end
        
        % PLOT
        fig2=figure;
        ax=subplot(1,2,1);
        imagesc(corrmat);colorbar;colormap(cmp); caxis([-1 1])
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
        xticks(1:length(plot_names)); xticklabels(plot_names); xtickangle(90);
        yticks(1:length(plot_names)); yticklabels(plot_names);
        title('collinearity of predictors (LMM)')
        set(ax, 'XAxisLocation', 'top');
        save(fullfile(S.prep.path.outputs,'predictor_correlations.mat'),'corrmat','var_names');
        savefig(fig2,fullfile(S.prep.path.outputs,[S.prep.sname 'collinearity.fig']));
        
        % multicollinaerity
        coli=find(abs(corrmat)>S.prep.calc.pred.test_collinearity);
        [row,col] = ind2sub(size(corrmat),coli);
        ax2=subplot(1,2,2);
        imagesc(abs(corrmat)>S.prep.calc.pred.test_collinearity); caxis([0 1])
        xticks(1:length(plot_names)); xticklabels(plot_names); xtickangle(90);
        yticks(1:length(plot_names)); yticklabels(plot_names);
        title('multicollinearity of predictors (LMM)')
        set(ax2, 'XAxisLocation', 'top');
        for ci = 1:length(coli)
            var1=var_names{row(ci)};
            var2=var_names{col(ci)};
            disp(['collinear predictors (LMM): ' var1 ', ' var2])
        end
        if S.prep.calc.pred.allow_collinearity==false && ~isempty(coli)
            error('collinear predictors - see above')
        end
        
        
        % correlations on concatenated data
        vt = vartype('numeric');
        numericvar_names=D.prep.dtab(:,vt).Properties.VariableNames;
        numericvar_ind = find(ismember(var_names,numericvar_names));
        [corrmat,pmat] = corr(table2array(D(d).prep.dtab(:,var_names(numericvar_ind))));
%         VIF = diag(inv(corrmat))'; % VIF
        if strcmp(S.prep.calc.pred.output_metric,'r2')
            corrmat=corrmat.^2;
        end
        corrmat = tril(corrmat,-1);
        pmat = tril(pmat + pmat',-1);
        corrmat = corrmat.*double(pmat<S.prep.calc.pred.sig);
        if strcmp(S.prep.calc.pred.output_metric,'r')
            corrmat = corrmat.*double(corrmat>=S.prep.calc.pred.R_min | corrmat<-S.prep.calc.pred.R_min);
        elseif strcmp(S.prep.calc.pred.output_metric,'r2')
            corrmat = corrmat.*double(corrmat>=(S.prep.calc.pred.R_min)^2 | corrmat<-(S.prep.calc.pred.R_min)^2);
        end
        figure;
        ax=subplot(1,2,1);
        imagesc(corrmat);colorbar;colormap(cmp); caxis([-1 1])
        plot_names = var_names(numericvar_ind);
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
        xticks(1:length(plot_names)); xticklabels(plot_names); xtickangle(90);
        yticks(1:length(plot_names)); yticklabels(plot_names);
        title('collinearity of predictors (concat)')
        set(ax, 'XAxisLocation', 'top');
        %save(fullfile(S.prep.path.outputs,'predictor_correlations.mat'),'corrmat','var_names');
        
        % multicollinaerity
        coli=find(abs(corrmat)>S.prep.calc.pred.test_collinearity);
        [row,col] = ind2sub(size(corrmat),coli);
        ax2=subplot(1,2,2);
        imagesc(abs(corrmat)>S.prep.calc.pred.test_collinearity); caxis([0 1])
        xticks(1:length(plot_names)); xticklabels(plot_names); xtickangle(90);
        yticks(1:length(plot_names)); yticklabels(plot_names);
        title('multicollinearity of predictors (concat)')
        set(ax2, 'XAxisLocation', 'top');
        
%         % multicollinaerity VIF
%         figure
%         ax3=subplot(1,1,1);
%         if strcmp(S.prep.calc.pred.output_metric,'r2')
%             VIF = 1./(1-corrmat);
%         else
%             VIF = 1./(1-corrmat^2);
%         end
%         imagesc(VIF); caxis([1 10]); colorbar;
%         xticks(1:length(plot_names)); xticklabels(plot_names); xtickangle(90);
%         yticks(1:length(plot_names)); yticklabels(plot_names);
%         title('multicollinearity of predictors (concat - VIF)')
%         set(ax3, 'XAxisLocation', 'top');
        
        % scatterplot matrix - observe linearity
        fig3=figure;
        plot_names = var_names(numericvar_ind);
        [~,axpm] = plotmatrix(table2array(D(d).prep.dtab(:,plot_names)));
        for pn = 1:length(plot_names)
            ylabel(axpm(pn,1),plot_names{pn},'Rotation',0,'HorizontalAlignment','right')
            xlabel(axpm(end,pn),plot_names{pn},'Rotation',90,'HorizontalAlignment','right')
        end
        title('linearity of predictors')
        savefig(fig3,fullfile(S.prep.path.outputs,[S.prep.sname 'linearity.fig']));

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
end

% save data to disk
if S.prep.output.save
    disp('saving data to disk...')
    save(fullfile(S.prep.path.outputs,[S.prep.sname '.mat']),'D','S','-v7.3')
    disp('...done')
end


function out = predictor_PCA(dtab,S)
% dtab: all predictors in table format
% out: PCA coefficients etc.
% S: which predictors to include, which PCA method, etc.

% select predictors Y
for ym = 1:length(S.prep.calc.pred.PCA_cov) % for each Y model

    % find columns with relevant predictors
    col_idx=endsWith(dtab.Properties.VariableNames,S.prep.calc.pred.PCA_cov{ym});
    
    % get an index of those colidx overlapping with trialback
    trialback_idx{ym}=contains(dtab.Properties.VariableNames,"trialback") & col_idx;
    
    % remove any with "trialback" from col_idx
    col_idx(col_idx==trialback_idx{ym}) = 0;
%     
%     % remove any multicollinear variables >0.95
%     vt = vartype('numeric');
%     toRemove = {};
%     correlationMatrix = corrcoef(table2array(dtab(:,vt)));
%     for i = 2:size(correlationMatrix, 1)
%         if any(abs(correlationMatrix(i, i+1:end)) > 0.99)
%             i
%             toRemove = [toRemove, dtab(:,vt).Properties.VariableNames{i}];
%         end
%     end
%     col_idx(ismember(dtab.Properties.VariableNames,toRemove)) = 0;


    out.col_names{ym} = dtab.Properties.VariableNames(col_idx);

    for nc = 1:length(out.col_names{ym})
        Y{ym}(:,nc)=dtab.(out.col_names{ym}{nc});
    end
end


if exist('Y','var')
    out.Y = Y;
else
    out.Y=[];
    return
end


for ym = 1:length(S.prep.calc.pred.PCA_cov) % for each Y model

    % remove variables that don't exist
    varidx = endsWith(S.LMM_varnames,S.prep.calc.pred.PCA_cov{ym});
    PCA_cov{ym} = S.LMM_varnames(varidx);

    % remove collinear variables fro stability
    if S.prep.calc.pred.PCA_LMM
        varidx = endsWith(S.LMM_varnames,PCA_cov{ym});
        CM = S.LMM_corrmat(varidx,varidx)';
    else
        CM = cov(Y{ym},'partialrows'); 
    end
    rm = any(triu(CM.*~eye(size(CM)))>0.99,1);
    Y{ym} = Y{ym}(:,~rm);
    PCA_cov{ym} = PCA_cov{ym}(:,~rm);
    out.col_names{ym} = out.col_names{ym}(~rm);
end

% Do the same for the correlation matrix if needed
if ym==1 && S.prep.calc.pred.PCA_LMM
    varidx = endsWith(S.LMM_varnames,PCA_cov{ym});
    S.LMM_corrmat = S.LMM_corrmat(varidx,varidx);
    S.LMM_corrmat = S.LMM_corrmat+S.LMM_corrmat'+eye(length(S.LMM_corrmat));
end

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
                nfac_temp = find(cumsum(explained)/sum(explained) > S.prep.calc.pred.PCA_minvar);
                if isempty(nfac_temp)
                    nfac_temp = ncomp;
                end
                nfac = nfac_temp(1);
            else
                for col = 1:size(Y{ym},2)
                    rand_data(:,col) = Y{ym}(randperm(size(Y{ym},1)),col);
                end
                [~,~,~,~,explainedrand] = pca(rand_data,'NumComponents',ncomp,'Centered',false,'Algorithm','eig');
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
                nfac_temp = find(cumsum(explained)/sum(explained) > S.prep.calc.pred.PCA_minvar);
                if isempty(nfac_temp)
                    nfac_temp = ncomp;
                end
                nfac = nfac_temp(1);
            else
                for col = 1:size(Y{ym},2)
                    rand_data(:,col) = Y{ym}(randperm(size(Y{ym},1)),col);
                end
                [~, values] = spca(rand_data, [], ncomp, delta, -ncomp, maxiter, convergenceCriterion, verbose);
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

            % KMO test
            disp('model 1, KMO test')
            [A,B] = kmo(Y{ym});
            
            % inital num components
            ncomp=size(Y{ym},2); 
            
            % PCA
            if S.prep.calc.pred.PCA_LMM && ~isempty(S.LMM_corrmat)
                % feed in correlation matrix from LMM
                FactorResults = erp_pca(Y{ym}, ncomp, type, S.LMM_corrmat);
            else
                FactorResults = erp_pca(Y{ym}, ncomp, type);
            end
            explained = FactorResults.facVar;
            
            % select num components
            if S.prep.calc.pred.PCA_minvar
                nfac_temp = find(cumsum(explained)/sum(explained) > S.prep.calc.pred.PCA_minvar);
                if isempty(nfac_temp)
                    nfac_temp = ncomp;
                end
                nfac1 = nfac_temp(1);
            else 
                nfac1 = 0;
            end
            
            for col = 1:size(Y{ym},2)
                rand_data(:,col) = Y{ym}(randperm(size(Y{ym},1)),col);
            end
            if S.prep.calc.pred.PCA_LMM && ~isempty(S.LMM_corrmat)
                % feed in correlation matrix from LMM
                randResults = erp_pca(rand_data, ncomp, type, S.LMM_corrmat);
            else
                randResults = erp_pca(rand_data, ncomp, type);
            end
            explainedrand = randResults.facVar;
            %explainedrand = explainedrand*(sum(explained)/sum(explainedrand)); %scale to same size as real data
            nfac_temp = find(explained<explainedrand);
            if isempty(nfac_temp)
                nfac2=1;
            else
                nfac2 = max(1,nfac_temp(1)-1);
            end

            nfac = max([nfac1 nfac2]);

            if S.prep.calc.pred.PCA_LMM && ~isempty(S.LMM_corrmat)
                % feed in correlation matrix from LMM
                FactorResults = erp_pca(Y{ym}, nfac, type, S.LMM_corrmat);
            else
                FactorResults = erp_pca(Y{ym}, nfac, type);
            end
            out.Ycoeff{ym} = FactorResults.FacCof;
            out.Ycorr{ym} = FactorResults.FacCor;
            out.Ypattern{ym} = FactorResults.FacPat;
            out.Ystruct{ym} = FactorResults.FacStr;
            out.Yscore{ym} = FactorResults.FacScr;
            out.Yfacvar{ym} = FactorResults.facVar;
            out.Yfacvaruniq{ym} = FactorResults.facVarQ;
            out.Yfacvartotal{ym} = FactorResults.facVarTot;
            
        end
    case 'PLS'
        % here we consider Y as X, i.e. predictors
        Xall=Y;
        clear Y;
        % and the data to be predicted is some EEG data components,
        % centered (e.g. CCA components). cdata{i} is one subject's data
        % consisting of trials (rows) by components (cols)
        if isfield(S.prep.calc.pred.pls,'eeg_data') && ~isempty(S.prep.calc.pred.pls.eeg_data)
            disp('loading EEG data')
            Dy = load(S.prep.calc.pred.pls.eeg_data,'D');
            for nY = 1:length(Dy.D.prep.Y)
                Y(:,nY) = Dy.D.prep.Y(nY).dtab.data;
            end
        else
            error('supply EEG data')
        end
        
        % concatenate over subjects? Or treat each separately
        if S.prep.calc.pred.pls.concat_over_subjects
            cdata{1} = Y; 
            for xm = 1:length(Xall)
                X{xm}{1} = Xall{xm};
            end
        else
            [U,~,iU]=unique(Dy.D.prep.dtab.ID,'stable');
            for i=1:length(U)
                cdata{i} = Y(iU==i,:);
                for xm = 1:length(Xall)
                    X{xm}{i} = Xall{xm}(iU==i,:);
                end
            end
        end
        clear Dy
        
        for xm = 1:length(X) % for each X model, i.e. predictor model
            
            if S.prep.calc.pred.pls.select_ncomp.boot

            elseif S.prep.calc.pred.pls.select_ncomp.MSE_CV || S.prep.calc.pred.pls.select_ncomp.frac_explained || S.prep.calc.pred.pls.select_ncomp.PRESS
                
                clear explained explainedrand
                NUM_FAC = size(X{xm}{1},2);
                for i = 1:length(cdata)
                    [~,~,~,~,BETA{xm}{i},explained{xm}(:,:,i),MSEcv{xm}(:,:,i)] = plsregress(X{xm}{i},cdata{i},NUM_FAC,'cv',5);
                    [~,~,~,~,~,~,MSE{xm}(:,:,i)] = plsregress(X{xm}{i},cdata{i},NUM_FAC);

                    % scale random data
                    si = std(cdata{i}(:));
                    randcdata=randn(size(cdata{i}))*si;
                    [~,~,~,~,BETArand{xm}{i},explainedrand{xm}(:,:,i),MSEcvrand{xm}(:,:,i)] = plsregress(X{xm}{i},randcdata,NUM_FAC,'cv',5);
                    [~,~,~,~,~,~,MSErand{xm}(:,:,i)] = plsregress(X{xm}{i},randcdata,NUM_FAC);

                end
                
                % find number of PLS components using comparison of CV-MSE
                % between actual and random data; or use fraction explained
                nfac_msecv = find(mean(MSEcv{xm}(2,:,:),3)<mean(MSEcvrand{xm}(2,:,:),3));

                % find number of PLS components using PRESS
                % https://stats.idre.ucla.edu/wp-content/uploads/2016/02/pls.pdf
                % http://facweb.cs.depaul.edu/sjost/csc423/documents/f-test-reg.htm
                % https://sites.duke.edu/bossbackup/files/2013/02/FTestTutorial.pdf
                
                nsamples = [];
                for i = 1:length(cdata)
                    nsamples(i) = size(cdata{i},1);
                end
                mean_nsamples = mean(nsamples);
                
                % PRESS compared to previous model: Y
                Fy=[];Py=[];
                for n = 1:NUM_FAC
                    df1 = mean_nsamples - (n-1);
                    df2 = mean_nsamples - n;
                    SSR1 = mean(MSEcv{xm}(2,n,:),3)*df1;
                    SSR2 = mean(MSEcv{xm}(2,n+1,:),3)*df2;
                    %F = SSR1/SSR2;
                    Fy(n) = ((SSR1-SSR2)/(df1-df2)) / (SSR2/df2);
                    Py(n)=1-fcdf(Fy(n),df1,df2);
                end
                % PRESS compared to previous model: X
                Fx=[];Px=[];
                for n = 1:NUM_FAC
                    df1 = mean_nsamples - (n-1);
                    df2 = mean_nsamples - n;
                    SSR1 = mean(MSEcv{xm}(1,n,:),3)*df1;
                    SSR2 = mean(MSEcv{xm}(1,n+1,:),3)*df2;
                    %F = SSR1/SSR2;
                    Fx(n) = ((SSR1-SSR2)/(df1-df2)) / (SSR2/df2);
                    Px(n)=1-fcdf(Fx(n),df1,df2);
                end
                
%                 % PRESS compared to random: Y
%                 F=[];P=[];
%                 for n = 1:NUM_FAC
%                     df1 = mean_nsamples - n;
%                     df2 = mean_nsamples - n;
%                     SSR1 = mean(MSEcvrand{xm}(2,n,:),3)*df1;
%                     SSR2 = mean(MSEcv{xm}(2,n+1,:),3)*df2;
%                     F(n) = SSR1/SSR2;
% %                     F(n) = ((SSR1-SSR2)/(df1-df2)) / (SSR2/df2);
%                     P(n)=1-fcdf(F(n),df1,df2);
%                 end
%                 % PRESS compared to random: X
%                 F=[];P=[];
%                 for n = 1:NUM_FAC
%                     df1 = mean_nsamples - n;
%                     df2 = mean_nsamples - n;
%                     SSR1 = mean(MSEcvrand{xm}(1,n,:),3)*df1;
%                     SSR2 = mean(MSEcv{xm}(1,n+1,:),3)*df2;
%                     F(n) = SSR1/SSR2;
% %                     F(n) = ((SSR1-SSR2)/(df1-df2)) / (SSR2/df2);
%                     P(n)=1-fcdf(F(n),df1,df2);
%                 end
                
                figure('name',['PLS MSE and explained variance: Model ' num2str(xm)])
                subplot(2,4,1); hold on
                    plot(0:NUM_FAC,mean(MSE{xm}(1,:,:),3)'); xlabel('n comp'), ylabel('X MSE'); %gcaExpandable;
                    plot(0:NUM_FAC,mean(MSErand{xm}(1,:,:),3)','r--'); 
                subplot(2,4,5); hold on
                    plot(0:NUM_FAC,mean(MSE{xm}(2,:,:),3)'); xlabel('n comp'), ylabel('Y MSE'); %gcaExpandable;
                    plot(0:NUM_FAC,mean(MSErand{xm}(2,:,:),3)','r--'); 
                subplot(2,4,2); hold on
                    plot(0:NUM_FAC,mean(MSEcv{xm}(1,:,:),3)'); xlabel('n comp'), ylabel('X CV MSE'); %gcaExpandable;
                    plot(0:NUM_FAC,mean(MSEcvrand{xm}(1,:,:),3)','r--');
                    plot([nfac_msecv(end),nfac_msecv(end)],[0 max(mean(MSEcv{xm}(1,:,:),3))],'k'); 
                subplot(2,4,6); hold on
                    plot(0:NUM_FAC,mean(MSEcv{xm}(2,:,:),3)'); xlabel('n comp'), ylabel('Y CV MSE'); %gcaExpandable;
                    plot(0:NUM_FAC,mean(MSEcvrand{xm}(2,:,:),3)','r--');
                    plot([nfac_msecv(end),nfac_msecv(end)],[0 max(mean(MSEcv{xm}(2,:,:),3))],'k'); 
                subplot(2,4,3); hold on
                    plot(1:NUM_FAC,cumsum(mean(explained{xm}(1,:,:),3)')); xlabel('n comp'), ylabel('X explained variance'); %gcaExpandable;
                    plot(1:NUM_FAC,cumsum(mean(explainedrand{xm}(1,:,:),3))','r--'); 
                subplot(2,4,7); hold on
                    plot(1:NUM_FAC,cumsum(mean(explained{xm}(2,:,:),3)')); xlabel('n comp'), ylabel('Y explained variance'); %gcaExpandable;
                    plot(1:NUM_FAC,cumsum(mean(explainedrand{xm}(2,:,:),3))','r--'); 
                fi = 1:NUM_FAC;
                subplot(2,4,4); hold on
                    plot(1:NUM_FAC,Fx); xlabel('n comp'), ylabel('X F value vs previous'); %gcaExpandable;
                    plot(fi(Px<0.05),Fx(Px<0.05),'ro'); 
                subplot(2,4,8); hold on
                    plot(1:NUM_FAC,Fy); xlabel('n comp'), ylabel('Y F value vs previous'); %gcaExpandable;
                    plot(fi(Py<0.05),Fy(Py<0.05),'ro'); 

                nr=ceil(sqrt(length(cdata)));
                figure('name',['PLS fitted responses: Model ' num2str(xm)])
                for i = 1:length(cdata)
                    subplot(nr,nr,i); plot(cdata{i},[ones(size(X{xm}{i},1),1) X{xm}{i}]*BETA{xm}{i},'bo');
                        xlabel('Observed Response');
                        ylabel('Fitted Response');
                end
                figure('name',['PLS fitted responses: random: Model ' num2str(xm)])
                for i = 1:length(cdata)
                    subplot(nr,nr,i); plot(cdata{i},[ones(size(X{xm}{i},1),1) X{xm}{i}]*BETArand{xm}{i},'bo');
                        xlabel('Observed Response');
                        ylabel('Fitted Response');
                end
                
                if S.prep.calc.pred.pls.select_ncomp.MSE_CV
                    grp_nfac{xm} = nfac_msecv(end)-1;
                elseif S.prep.calc.pred.pls.select_ncomp.frac_explained
                    nfac_exp = find(cumsum(mean(explained{xm}(2,:,:),3)) > S.prep.calc.pred.pls.select_ncomp.frac_explained);
                    if isempty(nfac_exp)
                    	nfac_exp = NUM_FAC;
                    end
                    grp_nfac{xm} = nfac_exp(1);
                else 
                    grp_nfac{xm} = str2double(inputdlg('number of factors?'));
                end
            else
                grp_nfac{xm} = NUM_FAC;
            end
            
        end
        
        if 0
            % select model based on best t value
            for xm = 1:length(X) % for each Y model
                meantL(xm) = mean(abs(tL{xm}(:)));
            end
            [~,minidx] = find(max(meantL));
            
        elseif 0
            % select model based on minimum MSE
            for xm = 1:length(X) % for each Y model
                MSEcvxm(xm) = mean(MSEcv{xm}(2,grp_nfac{xm},:),3)/size(X{xm}{1},2); 
            end
            [minMSE, minidx] = min(MSEcvxm);
            figure; hold on
                    plot(1:length(MSEcvxm),MSEcvxm); xlabel('model'), ylabel('mean MSE per variable'); %gcaExpandable;
                    scatter(minidx,minMSE,'k'); 
        else
            minidx = 1:length(X);
        end
        
        % take grp nfac
        for xm = minidx
            NUM_FAC = grp_nfac{xm};
            clear explained MSE MSE_CV W
            for i = 1:length(cdata)
                [XCOEFF,YCOEFF,Xscore,Yscore,BETA{i},explained(:,:,i),MSE(:,:,i),stats] = plsregress(X{xm}{i},cdata{i},NUM_FAC);
                [~,~,~,~,~,~,MSE_CV(:,:,i),~] = plsregress(X{xm}{i},cdata{i},NUM_FAC,'cv',5);

                out.Ycoeff{xm}{i}=YCOEFF;
                out.Yscore{xm}{i}=Yscore;
                out.Xcoeff{xm}{i}=XCOEFF;
                out.Xscore{xm}{i}=Xscore;
                out.BETA{xm}{i}=BETA{i};
                out.MSE{xm}{i}=MSE(:,:,i);
                out.MSE_CV{xm}{i}=MSE_CV(:,:,i);
                out.explained{xm}{i}=explained(:,:,i);
                out.W{xm}{i}=stats.W;
                out.Yfitted{xm}{i} = [ones(size(X{xm}{i},1),1) X{xm}{i}]*BETA{i};

               
            end
            out.Ymodel=xm;
    %         if exist('YCOEFFPCA','var')
    %             out.YCOEFF_PCA=YCOEFFPCA{xm};
    %             out.Yscore_PCA=YscorePCA{xm};
    %         end
            out.grp_nfac = grp_nfac;

            % only plot if there are multiple subjects (not concat)
            if ndims(explained)==3
                figure('name',['PLS MSE and explained variance: Model ' num2str(xm)])
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
            end
            figure('name',['PLS fitted responses: first ' num2str(NUM_FAC) ' components'])
            nr=ceil(sqrt(length(cdata)));
            for i = 1:length(cdata)
                subplot(nr,nr,i); plot(cdata{i},out.Yfitted{xm}{i},'bo');
                    xlabel('Observed Response');
                    ylabel('Fitted Response');
            end
            
            if S.prep.calc.pred.pls.rotate_factors
%                 % Rotate away from the principal components
%                 [out.Ycoeff{xm}{i},T] = rotatefactors(out.Ycoeff{xm}{i},'Method',S.prep.calc.pred.FA_type.rotation_option,'Power',S.prep.calc.pred.FA_type.degree);
%                 % Variances of the rotated scores, if we have S
%                 evalues = diag(cov(out.Yscore{xm}{i}*T))';
%                 out.explained{xm}{i} = diag(evalues)/sum(diag(evalues));
%                 %FactorResults = erp_pca(Y{ym},nfac,type);
%                 out.Ycorr{ym} = FactorResults.FacCor;
%                 out.Ypattern{ym} = FactorResults.FacPat;
%                 out.Ystruct{ym} = FactorResults.FacStr;
%                 out.Yscore{ym} = FactorResults.FacScr;
%                 out.Yfacvar{ym} = FactorResults.facVar;
%                 out.Yfacvaruniq{ym} = FactorResults.facVarQ;
%                 out.Yfacvartotal{ym} = FactorResults.facVarTot;
                for i = 1:length(cdata)
                    
                    data=X{xm}{i};
                    R=cov(data,'partialrows');
                    R=cov2corr(R);
                    goodVars=find(std(data) ~= 0);
                    Sc = cov(data,'partialrows'); 
                    Sd = diag(sqrt(diag(Sc)));
                    
                    V=out.Xcoeff{xm}{i};
                    ScrDiag = diag(nanstd(out.Xscore{xm}{i}));
                    A = inv(Sd) * (V * ScrDiag);  %unrotated factor loading matrix
                    
                    type.rotation_option = S.prep.calc.pred.FA_type.rotation_option;
                    type.degree = S.prep.calc.pred.FA_type.degree;
                    if strcmp(type.rotation_option,'promax')
                        [FacPat, FacCor] = rotatefactors(A,...
                            'Method',type.rotation_option,... 
                            'Power',type.degree,... % Promax only
                            'Normalize','on',...
                            'Reltol',.00001,...
                            'Maxit',1000);
                        FacStr = FacPat * FacCor;	%factor structure matrix (Harman, eq. 12.19, p. 268)
                    else
                        [FacPat, FacCor] = rotatefactors(A,...
                            'Method',type.rotation_option,... % 'Method' is 'orthomax', 'varimax', 'quartimax', 'equamax', or 'parsimax'
                            'Normalize','on',...
                            'Reltol',.00001,...
                            'Maxit',1000);
                        FacStr = FacPat;
                    end
                    LargestLoading=max(max(abs(FacStr))); %factor pattern loadings can go over 1 for oblique rotations
                    LargestCom=max(max(abs(sum(FacPat.*FacStr,2))));


                    invR=pinv(R);
                    FacCof=invR*FacStr;
                    FacScr=(data)*FacCof;

                    facVar=[];
                    facVarQ=[];
                    Comm=[];

                    if LargestCom ~=0 && LargestLoading~=0

                        FacScr=(FacScr)*inv(diag(std(FacScr))); %Standardize factor scores, not mean corrected.
                        Var=Sd.^2;

                        Comm=sum((Var*FacPat.*FacStr),2)/sum(diag(Var)); % p. 100

                        facVar = sum((Var*FacPat.*FacStr))/sum(diag(Var));

                        facVarQ=sum(Var*(FacPat*inv(diag(sqrt(diag(inv(FacCor)))))).^2)/sum(diag(Var));

                        %sort factors in order of size
                        [dummy index]=sort(facVar); 
                        index = fliplr(index);

                        FacPat=FacPat(:,index);
                        FacStr=FacStr(:,index);
                        FacCof=FacCof(:,index);
                        FacScr=FacScr(:,index);
                        FacCor=FacCor(index,index);
                        facVar=facVar(index);
                        facVarQ=facVarQ(index);

                        for f1 = 1:NUM_FAC
                            if sum(FacPat(:,f1)) < 0	%flip factor loadings so mostly positive
                                FacPat(:,f1) = FacPat(:,f1) .* (-1);
                                FacStr(:,f1) = FacStr(:,f1) .* (-1);
                                FacCof(:,f1) = FacCof(:,f1) .* (-1);
                                FacScr(:,f1) = FacScr(:,f1) .* (-1); %if the loading is flipped, the scores must be flipped too.
                                FacCor(:,f1) = FacCor(:,f1) .* (-1);
                                FacCor(f1,:) = FacCor(f1,:) .* (-1);
                            end
                        end
                    end

                    facNames=[];
                    facTypes=[];

                    factorDigits=length(num2str(NUM_FAC));
                    for f=1:NUM_FAC
                        facNames{f,1}=['TF' sprintf(['%0' num2str(factorDigits) 'd'],f)];
                        facTypes{f,1}='SGL';
                    end

                    out.Ycoeff{xm}{i}(goodVars,:) = FacCof;
                    out.Ycorr{xm}{i} = FacCor;
                    out.Ypattern{xm}{i}(goodVars,:) = FacPat;
                    out.Ystruct{xm}{i}(goodVars,:) = FacStr;
                    out.Yscore{xm}{i} = FacScr;
                    out.Yfacvar{xm}{i} = facVar;
                    out.Yfacvaruniq{xm}{i} = facVarQ;
                    out.Yfacvartotal{xm}{i} = sum(Comm,1);
                    out.numFacs{xm}{i}=NUM_FAC;
                    out.facNames{xm}{i}=facNames;
                    out.facTypes{xm}{i}=facTypes;
                end
            end
            
        end
        
end

if S.prep.calc.pred.CCA_over_subjects

    % trial indices
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
    % obtain subject-specific PCs
    PCdata = nan(NUM_FAC(1),nobs);
    sigma = {};
    for i = 1:length(cdata)
        sigma{i} = sqrt(var(score{i})); 
        if length(S.PCAmethod)>1 && any(strcmp(S.PCAmethod{2},'CCA')) % standardise scores if they will be used for CCA
           score{i} = score{i}./repmat(sqrt(var(score{i})),size(score{i},1),1); % standardise so that each PC from each subject contributes equally to CCA solution
        end
        PCdata(:,rep{i},i) = score{i}';
    end

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
        ncomp =  1 : ceil(size(O(o).scores,1)/min(size(O(o).scores,1),S.cca_test_nPCAcomp)) : size(O(o).scores,1);
    else
        ncomp = size(O(o).scores,1);
    end
    CV_fold=5;
    CVsample = cvpartition(size(O(o).scores,2),'KFold',CV_fold);

    nci=0;
    ncii=[];
    RMSECV=[];
    cca_reduce_by_scree = S.cca_reduce_by_scree;
    for nc = 1:length(ncomp) 

        % run CV CCAs in parallel for speed
        CV_W = cell(CV_fold,nr);
        CV_mdata_test = cell(CV_fold,nr);
        nCCA = min([ncomp(nc),NUM_FAC(2)]);
        nCCA_train = nan(CV_fold,nr);
        disp(['Estimating CCAs: comp ' num2str(nc) '/' num2str(length(ncomp))])
        parfor ri = 1:nr
            for cv = 1:CV_fold
                %disp(['ri ' num2str(ri) ', cv ' num2str(cv)])
                CVdata_train = O(o).scores(1:ncomp(nc),training(CVsample,cv),:);
                CVdata_test = O(o).scores(1:ncomp(nc),test(CVsample,cv),:);

                % train sample
                [CV_W{cv,ri},~] = mccas(CVdata_train,nCCA,reg,r(ri),O(o).W(1:ncomp(nc),:,:),S,cca_reduce_by_scree);
                nCCA_train(cv,ri) = size(CV_W{cv,ri},2);

                % test sample: fix the number of eigenvectors
                [~,CV_mdata_test{cv,ri}] = mccas(CVdata_test,nCCA_train(cv,ri),reg,r(ri),O(o).W(1:ncomp(nc),:,:),S,0);
                nCCA_train(cv,ri) = size(CV_mdata_test{cv,ri},1);
            end
        end

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
            for u=1:length(U)
                for cv = 1:CV_fold
                    for ri = 1:nr

                        % reconstructed PCs
                        PCpred = squeeze(CV_mdata_test{cv,ri}(1:ncCCA,:,u))'*pinv(squeeze(CV_W{cv,ri}(:,1:ncCCA,u)));
                        % reconstructed EEG data
                        testind = test(CVsample,cv);
                        data_test = subdata{o}{u}(testind(rep{u}),:);
                        data_pred = PCpred(ismember(find(testind),rep{u}),:)*pinv(O(o).W(1:ncomp(nc),:,u))';
                        sq_diff=bsxfun(@minus,data_test,data_pred).^2;
                        RMSECV(ri,cv,nci,u) = squeeze(sqrt(nanmean(sq_diff(:))));
                        %PCpred{u} = CV_mdata_test(:,:,u)'*pinv(CV_W(:,:,u));
                        %sq_diff=bsxfun(@minus,CVdata_test(:,:,u),PCpred{u}').^2;
                        %RMSECV(ri,cv,u) = squeeze(sqrt(nanmean(sq_diff(:))));
                    end
                end
            end
        end
    end

    % which regularisation parameter?
    RMSECVm = squeeze(mean(RMSECV,[2 4]));
    [~,minRi] = min(RMSECVm,[],1);
    minR=r(minRi);
    RMSECVr =[];
    for nc = 1:length(minRi)
        RMSECVr(:,nc,:) = RMSECV(minRi(nc),:,nc,:);
    end

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
            [mintemp,minYi(g)] = min(rmdat);
            for nc = 1:length(ncomp)
               plot(rmdat(ncii(:,1)==nc));
            end
            scatter(ncii(minYi(g),2),mintemp);
        end
        minYi = max(minYi);
    else
        rmdat = squeeze(mean(RMSECVr,[1 3]));
        [mintemp,minYi] = min(rmdat);
        for nc = 1:length(ncomp)
           plot(rmdat(ncii(:,1)==nc));
        end
        scatter(ncii(minYi,2),mintemp);
    end
    hold off

    nCCA = ncii(minYi,2);
    nPCA = ncomp(ncii(minYi,1));
    disp(['n CCA comp: ' num2str(nCCA) ', n PCA comp: ' num2str(nPCA)])
    O(o).CCA_regularisation = minR(minYi);
    O(o).RMSECV = rmdat;

    % reduce PCs
    O(o).scores = O(o).scores(1:nPCA,:,:);
    O(o).W = O(o).W(1:nPCA,:,:);
    for u=1:length(U)
        O(o).COEFF{u} = O(o).COEFF{u}(:,1:nPCA);
        O(o).sigma{u} = O(o).sigma{u}(:,1:nPCA);
    end

    % final CCA solution
    [O(o).W,O(o).mdata] = mccas(O(o).scores,nPCA,reg,minR(minYi),O(o).W,S,0);
    O(o).W = O(o).W(:,1:nCCA,:);
    O(o).mdata = O(o).mdata(1:nCCA,:,:);
    NUM_FAC(2) = nCCA;

    sz=size(O(o).mdata);
    O(o).mdata_avg = squeeze(mean(reshape(O(o).mdata,sz(1),nreps,obsdim,sz(3)),2));
elseif S.prep.calc.pred.concat_over_subjects && iscell(out.Ycoeff{1})
    for xm = 1:length(out.Ycoeff)
        % concatenate over subjects
        if strcmp(S.prep.calc.pred.PCA_type,'FA') || (strcmp(S.prep.calc.pred.PCA_type,'PLS') && S.prep.calc.pred.pls.rotate_factors)
            out.Ycoeff{xm} = out.Ycoeff{xm}{1};
            out.Ycorr{xm} = out.Ycorr{xm}{1} ;
            out.Ypattern{xm} = out.Ypattern{xm}{1};
            out.Ystruct{xm} = out.Ystruct{xm}{1};
            out.Yscore{xm} = out.Yscore{xm}{1};
            out.Yfacvar{xm} = out.Yfacvar{xm}{1};
            out.Yfacvaruniq{xm} = out.Yfacvaruniq{xm}{1};
            out.Yfacvartotal{xm} = out.Yfacvartotal{xm}{1};
            out.numFacs{xm} =  out.numFacs{xm}{1};
            out.facNames{xm} = out.facNames{xm}{1};
            out.facTypes{xm} = out.facTypes{xm}{1};
        elseif strcmp(S.prep.calc.pred.PCA_type,'PLS')
            out.Yscore{xm} = out.Yscore{xm}{1};
        end
    end
end

% apply weights to trialback variables
if ~isempty(trialback_idx{1}) && any(trialback_idx{1})
    Yscore_tb = {};
    for ym = 1:length(S.prep.calc.pred.PCA_cov) % for each Y model
    
        out.col_names_trialback{ym} = dtab.Properties.VariableNames(trialback_idx{ym});

        for nc = 1:length(out.col_names_trialback{ym})
            Yscore_tb{ym}(:,nc) = dtab.(out.col_names_trialback{ym}{nc});
        end
        Yscore_tb{ym} = Yscore_tb{ym} * out.Ycoeff{ym};
    end
    out.Yscore_tb = Yscore_tb;
end

function [A,B] = kmo(X)
%KMO Kaiser-Meyer-Olkin Measure of Sampling Adequacy.
% Factor analysis can be used as a guide to how coherently a set of variables
% relate to a hypothesized underlying dimension that they are all being used
% to measure. External validity analysis assesses whether the scale that has
% been constructed performs as theoretically expected in correlation with
% other variables to which it is expected to be related.
% There are some assumptions about the characteristics of factors that are
% extracted and defined that are unobserved common dimensions that may be
% listed to account for the correlations among observed variables. Sampling
% adequacy predicts if data are likely to factor well, based on correlation
% and partial correlation. Is used to assess which variables to drop from the
% model because they are too multicollinear.
% It has been suggested that inv(R) should be a near-diagonal matrix in order
% to successfully fit a factor analysis model.  To assess how close inv(R)
% is to a diagonal matrix, Kaiser (1970) proposed a measure of sampling
% adequacy, now called KMO (Kaiser-Meyer-Olkin) index. The common part, called
% the image of a variable, is defined as that part which is predictable by
% regressing each variable on all other variables.
% The anti-image is the specific part of the variable that cannot be predicted.
% Examining the anti-image of the correlation matrix. That is the negative of the
% partial correlations, partialling out all other variables.
% There is a KMO statistic for each individual variable and their sum is
% the overall statistic. If it is not > 0.6 drop the indicator variables with
% the lowest individual statistic value until the overall one rises above 0.6:
% factors which is meritorious. The diagonal elements on the Anti-image 
% correlation matrix are the KMO individual statistics for each variable. A KMO
% index <= 0.5 indicates the correlation matrix is not suitable for factor
% analysis.
%
% Syntax: function kmo(X) 
%      
%     Input:
%          X - Input matrix can be a data matrix (size n-data x p-variables)
%     Output(s):
%            - Kaiser-Meyer-Olkin Index.
%            - Degree of Common Variance Report (shared by a set of variables
%              and thus assesses the degree to which they measure a common
%              underlying factor).
%        optional(s):
%            - Anti-image Covariance Matrix.
%            - Anti-image Correlation Matrix
%
%  Example: From the example given on the web page
%  http://www.ncl.ac.uk/iss/statistics/docs/factoranalysis.html
%  We are interested to calculate the Kaiser-Meyer-Olkin measure of sampling
%  adequacy in order to see if proceeds a satisfactory factor analysis to
%  investigate the reasons why customers buy a product such as a particular
%  brand of soft drink (e.g. coca cola). Several variables were identified
%  which influence customer to buy coca cola. Some of the variables identified
%  as being influential include availability of product (X1), cost of product
%  (X2), experience with product (X3), popularity of product (X4), prestige
%  attached to product (X5), quality of product (X6), quantity of product (X7),
%  and respectability of product (X8). From this, you designed a questionnaire
%  to solicit customers' view on a seven point scale, where 1 = not important
%  and 7 = very important. The results from your questionnaire are show on the
%  table below. Only the first twelve respondents (cases) are used in this
%  example. 
%
%  Table 1: Customer survey
%  --------------------------------------------------
%     X1    X2    X3    X4    X5    X6    X7    X8   
%  --------------------------------------------------
%      4     1     4     5     2     3     6     7
%      4     2     7     6     6     3     3     4
%      6     4     3     4     2     5     7     7
%      5     3     5     4     3     4     6     7
%      5     2     4     5     2     5     5     6
%      6     3     3     5     4     4     7     7
%      6     2     4     4     4     3     4     5
%      4     1     3     4     3     3     5     6
%      5     3     4     3     4     3     6     6
%      5     4     3     4     4     4     6     7
%      6     2     4     4     4     3     7     5
%      5     2     3     3     3     3     7     6  
%  --------------------------------------------------
%
% Data matrix must be:
%  X=[4 1 4 5 2 3 6 7;4 2 7 6 6 3 3 4;6 4 3 4 2 5 7 7;5 3 5 4 3 4 6 7;
%  5 2 4 5 2 5 5 6;6 3 3 5 4 4 7 7;6 2 4 4 4 3 4 5;4 1 3 4 3 3 5 6;
%  5 3 4 3 4 3 6 6;5 4 3 4 4 4 6 7;6 2 4 4 4 3 7 5;5 2 3 3 3 3 7 6];
%
%  Calling on Matlab the function: 
%            kmo(X)
%
%  Answer is:
%
%  Kaiser-Meyer-Olkin Measure of Sampling Adequacy: 0.4172
%  The KMO test yields a degree of common variance unacceptable (Don't Factor).
%
%
%  Created by A. Trujillo-Ortiz, R. Hernandez-Walls, A. Castro-Perez, 
%             K. Barba-Rojo and A. Otero-Limon
%             Facultad de Ciencias Marinas
%             Universidad Autonoma de Baja California
%             Apdo. Postal 453
%             Ensenada, Baja California
%             Mexico.
%             atrujo@uabc.mx
%
%  Copyright. October 10, 2006.
%
% To cite this file, this would be an appropriate format:
% Trujillo-Ortiz, A., R. Hernandez-Walls, A. Castro-Perez, K. Barba-Rojo 
%   and A. Otero-Limon (2006). kmo:Kaiser-Meyer-Olkin Measure of Sampling
%   Adequacy. A MATLAB file. [WWW document]. URL http://www.mathworks.com/
%   matlabcentral/fileexchange/loadFile.do?objectId=12736
%
%  References:
%  Rencher, A. C. (2002), Methods of Multivariate Analysis. 2nd. ed.
%            New-Jersey:John Wiley & Sons. Chapter 13 (pp. 408-450).
%
error(nargchk(1,1,nargin));
msg = nargoutchk(1, 2, nargout);
X = corrcoef(X);
iX = inv(X);
S2 = diag(diag((iX.^-1)));
AIS = S2*iX*S2; %anti-image covariance matrix
IS = X+AIS-2*S2; %image covariance matrix
Dai = diag(diag(sqrt(AIS)));
IR = inv(Dai)*IS*inv(Dai); %image correlation matrix
AIR = inv(Dai)*AIS*inv(Dai); %anti-image correlation matrix
a = sum((AIR - diag(diag(AIR))).^2);
AA = sum(a);
b = sum((X - eye(size(X))).^2);
BB = sum(b);
MSA = b./(b+a); %measures of sampling adequacy
AIR = AIR-eye(size(AIR))+diag(MSA);
%Examine the anti-image of the correlation matrix. That is the negative of the partial correlations,
%partialling out all other variables.
N = BB;
D = AA+BB;
kmo = N/D;
disp(' ')
fprintf('Kaiser-Meyer-Olkin Measure of Sampling Adequacy: %3.4f\n', kmo);
if (kmo >= 0.00 && kmo < 0.50);
    disp('The KMO test yields a degree of common variance unacceptable (Don''t Factor).')
elseif (kmo >= 0.50 && kmo < 0.60);
    disp('The KMO test yields a degree of common variance miserable.')
elseif (kmo >= 0.60 && kmo < 0.70);
    disp('The KMO test yields a degree of common variance mediocre.')
elseif (kmo >= 0.70 && kmo < 0.80);
    disp('The KMO test yields a degree of common variance middling.')
elseif (kmo >= 0.80 && kmo < 0.90);
    disp('The KMO test yields a degree of common variance meritorious.')
elseif (kmo >= 0.90 && kmo <= 1.00);
    disp('The KMO test yields a degree of common variance marvelous.')
else
    disp('The KMO test result is out of range.')
end
if nargout == 1;
    disp(' ')
    disp('A = Anti-image covariance matrix.');    
    A = AIS;
elseif nargout > 1;
    disp(' ')
    disp('A = Anti-image covariance matrix.');
    A = AIS;
    disp('B = Anti-image correlation matrix.');
    B = AIR;
end

function [a,R,N]=cronbach(X)
% CRONBACH Cronbach's Alpha
%   [a,R,N]=cronbach(X) calculates the Cronbach's alpha of the data set X.
%
%   a is the Cronbach's alpha.
%   R is the upper triangular Spearman inter-correlation matrix
%   N is the number of valid items
%   X is the data set. Columns are the items.
%
%
%   Reference:
%   Cronbach L J (1951): Coefficient alpha and the internal structure of
%   tests. Psychometrika 16:297-333
%
%
%   If there are items with zero variance, these itemns are excluded from
%   the calculation. There is a warning in this case; alpha might be
%   guessed to high.
%
%
% Frederik Nagel
% Institute of Music Physiology and Musicians' Medicine
% Hanover University of Music and Drama 
% Hannover
% Germany
%
% e-mail: frederik.nagel@hmt-hannover.de
% homepage: http://www.immm.hmt-hannover.de
%
% April 24, 2006.
if nargin<1 | isempty(X)
   error('You shoud provide a data set.');
else
   % X must be at least a 2 dimensional matrix
   if size(X,2)<2
      error('Invalid data set.');
   end
end
% Items
N=size(X,2);
% Entries of the upper triangular matrix
e=(N*(N-1)/2);
% Spearman's correlation coefficient
R = corr(X,'rows','pairwise','type','spearman');
% Coefficients from upper triangular matrix
R = abs(triu(R,1));   
% Mean of correlation coefficients
r = sum(sum(triu(R,1)))/e;
% If there are columns with zero variance, these have to be excluded.
if(isnan(r))
    disp('There are columns with zero variance!');
    disp('These columns have been excluded from the calculation of alpha!');        
    disp([num2str(sum(sum(isnan(R)))) ' coefficients of ' num2str(N*N) ' have been excluded.']);    
    % Correct # of items 
    e = e-sum(sum(isnan(R)));
    
    % corrected mean of correlation coefficients 
    r = nansum(nansum(R))/e;    
    
    % corrected number of items
    N = N - sum(isnan(R(1,:)));
end
% Formular for alpha (Cronbach 1951)
a=(N*r)/(1+(N-1)*r);

function [corrmat,pmat,var_names] = LMM_correlation_matrix(D,S,test_collinearity_PCAcov_only)
% variables
if test_collinearity_PCAcov_only
    var_names=D.prep.dtab.Properties.VariableNames(endsWith(D.prep.dtab.Properties.VariableNames,S.prep.calc.pred.PCA_cov{1}));
else
    var_names = D.prep.dtab.Properties.VariableNames;
    var_names = setdiff(var_names,{'ID','group','test','train','eventTypes'},'stable');
end
vt = vartype('numeric');
numericvar_names=D.prep.dtab(:,vt).Properties.VariableNames;
numericvar_ind = find(ismember(var_names,numericvar_names));
nonnumericvar_ind = find(~ismember(var_names,numericvar_names));
all_ind = [numericvar_ind,nonnumericvar_ind]; % must be in this order

if length(var_names)<2
    corrmat=0;
    pmat=0;
    return
end

% combinations
comb = nchoosek(all_ind,2);   
comb = comb(ismember(comb(:,1),numericvar_ind),:); 
corrmat = diag(zeros(1,length(var_names)));
pmat = diag(zeros(1,length(var_names)));
for ci = 1:size(comb,1)
    formula = [var_names{comb(ci,1)} '~' var_names{comb(ci,2)} '+(1|ID)'];
    lme = fitlme(D.prep.dtab,formula);
    if strcmp(S.prep.calc.pred.output_metric,'r')
        corrmat(comb(ci,1),comb(ci,2)) = lme.Coefficients(2,2);
    elseif strcmp(S.prep.calc.pred.output_metric,'r2')
        corrmat(comb(ci,1),comb(ci,2)) = sign(double(lme.Coefficients(2,2)))*lme.Rsquared.Ordinary;
    end
    pmat(comb(ci,1),comb(ci,2)) = lme.Coefficients(2,6);
end
corrmat = tril(corrmat + corrmat',-1);
pmat = tril(pmat + pmat',-1);


function [W, mdata, mWeights] = mccas(data,K,reg,r,weights,Sx,cca_reduce_by_scree)
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
        if cca_reduce_by_scree
            
            allsubrand = randn(size(allsub'));
            Rrand = cov(allsubrand);
            Srand = zeros(size(Rrand));
            for i = 1:dim*sub
                tmp = ceil(i/dim);
                Srand((tmp-1)*dim+1:tmp*dim,i) = Rrand((tmp-1)*dim+1:tmp*dim,i); 
            end
            
            explained = diag(values)/sum(diag(values));
            [~,valuesrand]=eigs((Rrand-Srand),Srand,K);
            explainedrand = diag(valuesrand)/sum(diag(valuesrand));
            nfac_temp = find(explained<explainedrand);
            if ~isempty(nfac_temp)
                nfac = min(nfac_temp(1)-1, K);
                tempW = tempW(:,1:nfac);
            end
        end
    case 'FA'
        type = Sx.cca_FA_type;
        FactorResults = erp_pca(allsub',K,type,(R-S),S);
        tempW = out.FacCof;
        if cca_reduce_by_scree
            
            allsubrand = randn(size(allsub'));
            Rrand = cov(allsubrand);
            Srand = zeros(size(Rrand));
            for i = 1:dim*sub
                tmp = ceil(i/dim);
                Srand((tmp-1)*dim+1:tmp*dim,i) = Rrand((tmp-1)*dim+1:tmp*dim,i); 
            end
            
            FactorResultsRand = erp_pca(allsubrand,K,type,(Rrand-Srand),Srand);
            FactorResultsRand.scree = FactorResultsRand.scree*(sum(out.scree)/sum(FactorResultsRand.scree)); %scale to same size as real data
            nfac_temp = find(out.scree<FactorResultsRand.scree);
            if ~isempty(nfac_temp)
                nfac = min(nfac_temp(1)-1, K);
                tempW = tempW(:,1:nfac);
                
            end
        end
end

K = size(tempW,2);
W = zeros(dim,K,sub);
for i=1:sub
    if Sx.normalise_cca_weights 
        W(:,:,i)=tempW((i-1)*dim+1:i*dim,:)./norm(tempW((i-1)*dim+1:i*dim,:));
    else
        W(:,:,i)=tempW((i-1)*dim+1:i*dim,:);
    end
end

% projected data
mdata = nan(K,num,sub);
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
            mdatatemp(j,:) = mdatatemp(j,:)/norm(mdatatemp(j,:));
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
