function [results,corrmat,pmat]=eegstats_covariatecorrelations(S,varargin)
dbstop if error

%% find D saved from createimages or MCC functions
if isempty(varargin)
   % try loading the data if not supplied in varargin
   try
       load(fullfile(S.covcorr.path.inputs,'D.mat'),'D') 
   catch
       error('Supply a D structure or path/file of the data')
   end
else
    D= varargin{1};
end

% set paths
S.path.code = {
    };
set_paths(S)

if length(D)>1
    save_pref = 'subject_';
else
    save_pref = 'group_';
end

for d=1:length(D)
    
    types = S.covcorr.summary_types;
    
    % get model and contrast indices
    if isempty(S.covcorr.model.index)
        S.covcorr.model.index = 1:length(D(d).model);
    end
    if isempty(S.covcorr.model.contrast) && isfield(D(d).model(1),'con')
        for i = S.covcorr.model.index
            S.covcorr.model.contrast{i} = 1:length(D(d).model(i).con);
        end
    end
    
    % participant IDs
    ID=table;
    ID.ID = unique({D(d).model(1).random(:).Level},'stable')';

    % results table
    if ~isempty(S.covcorr.path.pred)
        results=readtable(S.covcorr.path.pred);
        npred = width(results);
    else
        results=table;
    end
    
    if S.covcorr.model.index
        for i = S.covcorr.model.index
            for c = S.covcorr.model.contrast{i}
                nc=numel(D(d).model(i).con(c).vox);
                for ci=1:nc
                    % input
                    if any(strcmp(types,'median'))
                        if ismember('input',S.covcorr.summary_data)
                            results.(['model' num2str(i) '-con' num2str(c) '-clus' num2str(ci) '-inputmedian']) = D(d).model(i).con(c).clus(ci).input_median_median';
                        end
                        if ismember('fitted',S.covcorr.summary_data)
                            results.(['model' num2str(i) '-con' num2str(c) '-clus' num2str(ci) '-fittedmedian']) = D(d).model(i).con(c).clus(ci).fitted_median_median';
                        end
                    end

                    if any(strcmp(types,'mean'))
                        if ismember('input',S.covcorr.summary_data)
                            results.(['model' num2str(i) '-con' num2str(c) '-clus' num2str(ci) '-inputmean']) = D(d).model(i).con(c).clus(ci).input_mean_mean';
                        end
                        if ismember('fitted',S.covcorr.summary_data)
                            results.(['model' num2str(i) '-con' num2str(c) '-clus' num2str(ci) '-fittedmean']) = D(d).model(i).con(c).clus(ci).fitted_mean_mean';
                        end
                    end

                    if any(strcmp(types,'eig'))
                        if ismember('input',S.covcorr.summary_data)
                            results.(['model' num2str(i) '-con' num2str(c) '-clus' num2str(ci) '-inputeig']) = D(d).model(i).con(c).clus(ci).input_eig_mean';
                        end
                        if ismember('fitted',S.covcorr.summary_data)
                            results.(['model' num2str(i) '-con' num2str(c) '-clus' num2str(ci) '-fittedeig']) = D(d).model(i).con(c).clus(ci).fitted_eig_mean';
                        end
                    end
                end
            end
            
            if ismember('coeffs',S.covcorr.summary_data)
                
                for co = 1:length(D(d).model(i).coeff)
                    inp_ind = find(strcmp(S.covcorr.summary_coeffs(:,2),D(d).model(i).coeff(co).name)); % input index
                    if isempty(inp_ind); continue; end
                    con_name = S.covcorr.summary_coeffs{inp_ind,1};
                    c = find(strcmp({D(d).model(i).con(:).term},con_name));
                    nc=numel(D(d).model(i).con(c).vox);
                    
                    for ci=1:nc
                        results.([strrep(con_name,':','-') '-clus' num2str(ci)]) = D(d).model(i).coeff(co).clus(ci).coeff_mean;
                    end
                end
            end
        end
    end
end

if S.covcorr.save_table
    IDresults = [ID results];
    writetable(IDresults,fullfile(S.covcorr.path.inputs,'covariates.xls'));
end
% experimental: ratios
% results.con45ratio1 = D.model.con(4).clus(3).input_median_median' - D.model.con(5).clus(1).input_median_median';


[corrmat, pmat] = corr(table2array(results),'type',S.covcorr.type);
plot_names = results.Properties.VariableNames;
if (S.covcorr.trim_matrix && ~isempty(S.covcorr.path.pred))
    figure
    imagesc(corrmat(1:npred,npred+1:end));colorbar;colormap('jet'); caxis([-1 1])
    xticks(1:width(results)-npred); xticklabels(plot_names(npred+1:width(results))); xtickangle(90);
    yticks(1:npred); yticklabels(plot_names(1:npred));
    title('predictor correlations - Rho values')
    figure
    imagesc(pmat(1:npred,npred+1:end)*S.covcorr.bonferroni);colorbar;colormap('jet'); caxis([0 0.05])
    xticks(1:width(results)); xticklabels(plot_names(npred+1:width(results))); xtickangle(90);
    yticks(1:npred); yticklabels(plot_names(1:npred));
    title('predictor correlations - p values')
else
    figure
    imagesc(corrmat);colorbar;colormap('jet'); caxis([-1 1])
    xticks(1:width(results)); xticklabels(plot_names); xtickangle(90);
    yticks(1:width(results)); yticklabels(plot_names);
    title('predictor correlations - Rho values')
    figure
    imagesc(pmat*S.covcorr.bonferroni);colorbar;colormap('jet'); caxis([0 0.05])
    xticks(1:width(results)); xticklabels(plot_names); xtickangle(90);
    yticks(1:width(results)); yticklabels(plot_names);
    title('predictor correlations - p values')
end

if 0
    %% multiple regression / lasso
    lm_results = table2array(results);
    lm_var = lm_results(:,npred+1:end);
    % collinearity of predictors
    thresh = 0.8;
    [corrmat, pmat] = corr(lm_var,'type',S.covcorr.type);
    plot_names = results(:,npred+1:end).Properties.VariableNames;
    figure
    imagesc(corrmat.*(corrmat>thresh));colorbar;colormap('jet'); caxis([-1 1])
    xticks(1:width(results)); xticklabels(plot_names); xtickangle(90);
    yticks(1:width(results)); yticklabels(plot_names);
    title(['predictor collinearity - Rho values > ' num2str(thresh)])
    % remove collinear predictors
    removetype = 2; % 1=rows, 2=cols
    lm_var(:,any(triu(corrmat>thresh,1),removetype))=[];
    plot_names(:,any(triu(corrmat>thresh,1),removetype))=[];
    % regress
    for p=1%:npred
        lm_pred = lm_results(:,p); 
%         [lm{p},fit{p}] = lasso(lm_var,lm_pred,'Alpha',0.5);
        lm{p} = fitlm(lm_var,lm_pred);
    end
end

function set_paths(S)
for p = 1:length(S.path.code)
    if S.path.code{p,1}
        addpath(genpath(S.path.code{p,2}));
    else
        addpath(S.path.code{p,2});
    end
end
