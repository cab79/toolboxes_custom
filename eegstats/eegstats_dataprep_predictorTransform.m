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

%% transformations: must take place on grouped data (over subjects) if analysis will be grouped
for d = 1:length(D)
    disp(['data transformations: ' num2str(d) '/' num2str(length(D))])
    
    % z-scoring on continuous predictors, but not if doing PLS on EEG data
    if S.prep.calc.pred.zscore
        %if S.prep.calc.eeg.pca.on == 1 && any(strcmp(S.prep.calc.eeg.pca.PCAmethod,'PLS'))
        %    disp('no z-scoring to properly rescale Y variables for PLS')
        %else
            vt = vartype('numeric');
            numericvar_names=D(d).prep.dtab(:,vt).Properties.VariableNames;
            for pr = 1:length(numericvar_names)
                if ismember(numericvar_names{pr}, S.prep.calc.pred.zscore_exclude); continue; end
                [zdata, D(d).prep.pred_means(pr), D(d).prep.pred_stds(pr)] = zscore(table2array(D(d).prep.dtab(:,numericvar_names{pr})));
                D(d).prep.dtab(:,numericvar_names{pr}) = array2table(zdata);
            end
            
            % additionally, convert any categorical variables that are due
            % for PCA to numeric and z-score them
            if S.prep.calc.pred.PCA_on
                ct = vartype('categorical');
                catvar_names=D(d).prep.dtab(:,ct).Properties.VariableNames;
                catvar_pca_idx = find(ismember(catvar_names,S.prep.calc.pred.PCA_cov{S.prep.calc.pred.PCA_model_add2dtab}));
                if numel(catvar_pca_idx)>0
                    for pr = 1:length(catvar_pca_idx)
                        if ismember(catvar_names{catvar_pca_idx(pr)}, S.prep.calc.pred.zscore_exclude); continue; end
                        [zdata, D(d).prep.pred_means(pr), D(d).prep.pred_stds(pr)] = zscore(double(table2array(D(d).prep.dtab(:,catvar_names{catvar_pca_idx(pr)}))));
                        D(d).prep.dtab(:,catvar_names{catvar_pca_idx(pr)}) = [];
                        D(d).prep.dtab(:,catvar_names{catvar_pca_idx(pr)}) = array2table(zdata);
                    end
                end
            end
        %end
    elseif S.prep.calc.pred.PCA_on == 1
        error('must z-score data prior to PCA')
    end
    
    if S.prep.calc.pred.PCA_on
            
        D(d).prep.pred_PCA = predictor_PCA(D(d).prep.dtab,S);

        % add PCA components to data table
        if S.prep.calc.pred.PCA_model_add2dtab
            ym=S.prep.calc.pred.PCA_model_add2dtab;
            for k = 1:size(D(d).prep.pred_PCA.Yscore{ym},2)
                pcname=[S.prep.calc.pred.PCA_type num2str(k)];
                D(d).prep.dtab.(pcname) = D(d).prep.pred_PCA.Yscore{ym}(:,k);
            end
            
            % trialback data
            if ~isempty(D(d).prep.pred_PCA.Yscore_tb)
                for k = 1:size(D(d).prep.pred_PCA.Yscore_tb{ym},2)
                    pcname=[S.prep.calc.pred.PCA_type '_trialback_' num2str(k)];
                    D(d).prep.dtab.(pcname) = D(d).prep.pred_PCA.Yscore_tb{ym}(:,k);
                end
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
