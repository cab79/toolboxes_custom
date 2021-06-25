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
    if isfield(S.prep.calc.pred,'log') && ~isempty(S.prep.calc.pred.log) 
        col_idx=find(contains(D(d).prep.dtab.Properties.VariableNames,S.prep.calc.pred.log));
        for pr = 1:length(col_idx)
            tempdat=table2array(D(d).prep.dtab(:,col_idx(pr)));
            if any(tempdat<0); tempdat = tempdat-min(tempdat); end
            tempdat(tempdat==0) = min(tempdat(tempdat~=0));
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
                elseif strcmp(S.prep.calc.pred.newvarmult{pr,4},'/')
                    temp = D(d).prep.dtab.(S.prep.calc.pred.newvarmult{pr,1})./D(d).prep.dtab.(S.prep.calc.pred.newvarmult{pr,2});
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
            
            % additionally, convert any categorical variables that are due
            % for PCA to numeric and z-score them
            if S.prep.calc.pred.PCA_on
                ct = vartype('categorical');
                catvar_names=D(d).prep.dtab(:,ct).Properties.VariableNames;
                for a2t = S.prep.calc.pred.PCA_model_add2dtab
                    catvar_pca_idx = find(ismember(catvar_names,S.prep.calc.pred.PCA_cov{a2t}));
                    if numel(catvar_pca_idx)>0
                        for pr = 1:length(catvar_pca_idx)
                            if ismember(catvar_names{catvar_pca_idx(pr)}, S.prep.calc.pred.zscore_exclude); continue; end
                            [zdata, D(d).prep.pred_means(pr,u), D(d).prep.pred_stds(pr,u)] = zscore(double(table2array(D(d).prep.dtab(iU==u,catvar_names{catvar_pca_idx(pr)}))));
                            D(d).prep.dtab(iU==u,catvar_names{catvar_pca_idx(pr)}) = [];
                            D(d).prep.dtab(iU==u,catvar_names{catvar_pca_idx(pr)}) = array2table(zdata);
                        end
                    end
                end
            end
        end
            
            
        %end
    %elseif S.prep.calc.pred.PCA_on == 1
    %    error('must z-score data prior to PCA')
    end
    
    
    
    if S.prep.calc.pred.PCA_on
            
        D(d).prep.pred_PCA = predictor_PCA(D(d).prep.dtab,S);
        
        if strcmp(S.prep.calc.pred.PCA_type,'FA')
            pattern = D(d).prep.pred_PCA.Ypattern{1};
            pattern_thresh = D(d).prep.pred_PCA.Ypattern{1}.*double(D(d).prep.pred_PCA.Ypattern{1}>=S.prep.calc.pred.PCA_model_minloading | D(d).prep.pred_PCA.Ypattern{1}<=-S.prep.calc.pred.PCA_model_minloading);
            structure = D(d).prep.pred_PCA.Ystruct{1};
            plot_names = D.prep.pred_PCA.col_names{1};
            figure;
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

            vardata=table2array(D(d).prep.dtab(:,contains(D(d).prep.dtab.Properties.VariableNames,S.prep.calc.pred.PCA_cov{1})));
            varnames=D(d).prep.dtab.Properties.VariableNames(contains(D(d).prep.dtab.Properties.VariableNames,S.prep.calc.pred.PCA_cov{1}));
            
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
            
            % add FA components to data table
            if S.prep.calc.pred.PCA_model_add2dtab
                for k = 1:size(D(d).prep.pred_PCA.Yscore{ym},2)
                    pcname=[S.prep.calc.pred.PCA_type num2str(k)];
                    D(d).prep.dtab.(pcname) = zscore((transpose(pattern_thresh(:,k))*transpose(vardata))');
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
                        D(d).prep.dtab.(pcname) = D(d).prep.pred_PCA.Yscore{ym}(:,k);
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
        title('factor regressions')
        set(ax, 'XAxisLocation', 'top');
        
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
        [corrmat,pmat] = corr(table2array(D(d).prep.dtab(:,var_names(numericvar_ind))));
        if strcmp(S.prep.calc.pred.output_metric,'r2')
            corrmat=corrmat.^2;
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
        
        % scatterplot matrix - observe linearity
        figure
        plot_names = var_names(numericvar_ind);
        [~,axpm] = plotmatrix(table2array(D(d).prep.dtab(:,plot_names)));
        for pn = 1:length(plot_names)
            ylabel(axpm(pn,1),plot_names{pn},'Rotation',0,'HorizontalAlignment','right')
            xlabel(axpm(end,pn),plot_names{pn},'Rotation',90,'HorizontalAlignment','right')
        end
        title('linearity of predictors')

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
    col_idx=contains(dtab.Properties.VariableNames,S.prep.calc.pred.PCA_cov{ym});
    
    % get an index of those colidx overlapping with trialback
    trialback_idx{ym}=contains(dtab.Properties.VariableNames,"trialback") & col_idx;
    
    % remove any with "trialback" from col_idx
    col_idx(col_idx==trialback_idx{ym}) = 0;
    
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
                nfac = nfac_temp(1);
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
           
            % KMO test
            disp('model 1, KMO test')
            [A,B] = kmo(Y{ym});
            
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
            out.Ycorr{ym} = FactorResults.FacCor;
            out.Ypattern{ym} = FactorResults.FacPat;
            out.Ystruct{ym} = FactorResults.FacStr;
            out.Yscore{ym} = FactorResults.FacScr;
            out.Yfacvar{ym} = FactorResults.facVar;
            out.Yfacvaruniq{ym} = FactorResults.facVarQ;
            out.Yfacvartotal{ym} = FactorResults.facVarTot;
            
        end
        

end

% apply weights to trialback variables
Yscore_tb = {};
for ym = 1:length(S.prep.calc.pred.PCA_cov) % for each Y model
    if ~isempty(trialback_idx{ym}) && any(trialback_idx{ym})
        out.col_names_trialback{ym} = dtab.Properties.VariableNames(trialback_idx{ym});

        for nc = 1:length(out.col_names_trialback{ym})
            Yscore_tb{ym}(:,nc) = dtab.(out.col_names_trialback{ym}{nc});
        end
        Yscore_tb{ym} = Yscore_tb{ym} * out.Ycoeff{ym};
    end
end
out.Yscore_tb = Yscore_tb;

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