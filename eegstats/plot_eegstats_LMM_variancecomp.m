function plot_eegstats_LMM_variancecomp(S)
% For calculation of ICCs/VPCs/R2.
% These can be interpreted as the correlation of two randomly drawn units 
% from the same grouping or as the amount of variance explained by those 
% groupings, similar to an R2.
% R2: called this when fixed effect variance is on numerator
% ICC/VPC: called this when random effect variance is on numerator
% https://stats.stackexchange.com/questions/178911/icc-in-a-multi-level-model-with-two-random-effects?rq=1

% GENERAL INFO:
% 1. ICC is calculated within-model. E.g. variance explained by ID relative
% to total variance (btw-ID + error). ICC is for empty model (but OK to control
% for between-subject covariates if needed. If including within-subject fixed
% effects and/or random slopes, this is no longer ICC but VPC.
% 2. VPC formula: https://www.bristol.ac.uk/media-library/sites/cmm/migrated/documents/variance-partitioning.pdf
% Currently the below code only allows one slope with two levels (0 and 1).
% 3. Can also calculate PCVs between models using only the change in random
% effect variance after adding fixed effects to the model. https://jech.bmj.com/content/jech/59/9/729.full.pdf
% e.g. change in between-ID variance after adding ODD. MUST HAVE THE SAME RANDOM EFFECTS STRUCTURE IN THe TWO MODELS 
% 4. Can also compare ID variances for models with different random effects
% (e.g intercept only vs. intercept+slope): the relative increase in the random
% effect covariance indicate how important the slope is to understanding individual
% differences.
% 5. To calculate R2 for fixed effects: https://stats.stackexchange.com/questions/128911/partitioning-explained-variance-to-fixed-effects-by-comparing-r-squared-r2-bet

dbstop if error

% close all
% set paths and file names 
S.path.code = {
   1, 'Q:\MATLAB\toolboxes_external\cbrewer'
   1, 'Q:\MATLAB\toolboxes_external\spm12'
    };

% set paths
set_paths(S)

% load cluster data
load(fullfile(S.path.stats_load,'D.mat'),'D');
cd(S.path.stats_load)
% 
% if ~isempty(S.model.def)
%     D.model(S.model.index).def = S.model.def{:};
% end

disp('loading data table')
Ddtab = load(S.path.dtab_inputs);
dtab = Ddtab.D.prep.dtab;

% find con
i = S.model.index;
try D.model(S.model.index)
catch
    error(['this model number does not exist'])
end
if isempty(S.model.clus)
    S.model.clus=1:100;
end
if ~isempty(S.model.contrast_term)
    if strcmp(S.model.contrast_term,'all')
        S.model.contrast_term = {D.model(S.model.index).con(:).term};
    elseif ~any(strcmp({D.model(S.model.index).con(:).term},S.model.contrast_term))
        error([S.model.contrast_term ' is not a constrast term in this model'])
    end
    cs = find(strcmp({D.model(S.model.index).con(:).term},S.model.contrast_term));
elseif ~isempty(S.model.contrastUC_term)
    if strcmp(S.model.contrastUC_term,'all')
        S.model.contrastUC_term = {D.model(S.model.index).conUC(:).term};
    elseif ~any(strcmp({D.model(S.model.index).conUC(:).term},S.model.contrastUC_term))
        error([S.model.contrastUC_term ' is not a constrast term in this model'])
    end
    cs = find(strcmp({D.model(S.model.index).conUC(:).term},S.model.contrastUC_term));
end

for c = cs

    if ~isempty(S.model.contrast_term)
        term = D.model(S.model.index).con(c).term;
        S.model.clus(S.model.clus>length(D.model(i).con(c).clus))=[];
    elseif ~isempty(S.model.contrastUC_term)
        term = D.model(S.model.index).conUC(c).term;
        S.model.clus(S.model.clus>length(D.model(i).conUC(c).clus))=[];
    end

    % per cluster
    for cl = S.model.clus

        if ~isempty(S.model.contrast_term)
            input=D.model(i).con(c).clus(cl).input_median';
        elseif ~isempty(S.model.contrastUC_term)
            input=D.model(i).conUC(c).clus(cl).input_median';
        end
        dtab.data=input;

        % estimate LMMs to obtain fixed and random effect variances
        for li = 1:length(S.model.def)
            disp(['estimating LMM, cluster ' num2str(cl)])
            lme=fitlme(dtab,S.model.def{li},'FitMethod',S.fitmethod,'DummyVarCoding',S.coding,'CheckHessian', true);

            disp(lme)
            anova(lme)

            % fixed effect coefficient variance
            fcov=sum(diag(lme.CoefficientCovariance));

            % calculate ICC (btw variance relative to sum of btw+error, with or without accounting for fixed effect variances)
            [psi,mse] = covarianceParameters(lme);
            if size(psi{1},1)<=2 % currently only works for slope with 2 levels max
                for p = 1:length(psi{1})
                    if p==1
                        model(li).bV{cl}(p) = psi{1}(1,1);
                    else
                        model(li).bV{cl}(p) = psi{1}(1,1) + 2*psi{1}(1,p)*(p-1) + psi{1}(p,p)*((p-1)^2);
                    end
                end
                % https://royalsocietypublishing.org/doi/10.1098/rsif.2017.0213
                model(li).ICC{cl} = model(li).bV{cl}./(model(li).bV{cl}+fcov+mse); % should include fixed effect on denominator unless calculating the adjusted ICC

                % calculate marginal R2 of fixed effect
                if strcmp(S.model.R2_var,'all')
                    model(li).R2{cl} = fcov./(model(li).bV{cl}+fcov+mse);
                else
                    cnames = lme.CoefficientNames;
                    ci=strcmp(cnames,S.model.R2_var);
                    pcov=sum(diag(lme.CoefficientCovariance(ci,ci)));
                    model(li).R2{cl} = pcov./(model(li).bV{cl}+fcov+mse);
                end
            end

        end

        for pc = 1:length(S.paired_comp)

            ii=0;
            for li = S.paired_comp{pc}
                ii=ii+1;
                cvar{ii} = model(li).bV{cl}; % more than one bV value - what to do?

                % ADD IN PCV FOR TRIAL VARIANCE (MSE)
    %             evar(i) = mse;

            end

            % The equation for the proportional change in
            % variance (PCV) is:
            model_comp(pc).PCVc{cl} = (cvar{1} - cvar{2}) ./ cvar{1};
    %         PCVe(pc) = (evar(1) - evar(2)) / evar(1);
        end
    %     % add nans for plotting
    %     lth=cellfun(@length,PCVc);
    %     ind=find(lth<max(lth));
    %     for ii = 1:length(ind)
    %         PCVc{ind(ii)} = [PCVc{ind(ii)},NaN];
    %     end
    %     lth=cellfun(@length,ICC);
    %     ind=find(lth<max(lth));
    %     for ii = 1:length(ind)
    %         ICC{ind(ii)} = [ICC{ind(ii)},NaN];
    %         R2{ind(ii)} = [R2{ind(ii)},NaN];
    %     end

    %     subplot(1,2,2); bar(PCVe); % just really tells you what the fixed effect variance 
    %     if exist('ICC','var')
    %         title(['ICC = ' num2str(ICC)]);
    %     end

    %     
    %     R2.fix.lme.total = sum(diag(lme.CoefficientCovariance),'all')/(sum(diag(lme.CoefficientCovariance),'all')+lme.MSE);
    %     R2.fixran.lme.total = (sum(diag(lme.CoefficientCovariance),'all')+sum(diag(rp{1}),'all'))/(sum(diag(lme.CoefficientCovariance),'all')+lme.MSE+sum(diag(rp{1}),'all'));
    %     for f = 1:length(lme.CoefficientCovariance)
    %         R2.fix.lme.fixed(f) = lme.CoefficientCovariance(f,f)/(sum(diag(lme.CoefficientCovariance),'all')+lme.MSE);
    %         R2.fixran.lme.fixed(f) = lme.CoefficientCovariance(f,f)/(sum(diag(lme.CoefficientCovariance),'all')+lme.MSE+sum(diag(rp{1}),'all'));
    %     end
    %     for r = 1:length(rp{1})
    %         R2.fix.lme.random(r) = rp{1}(r,r)/(sum(diag(lme.CoefficientCovariance),'all')+lme.MSE);
    %         R2.fixran.lme.random(r) = rp{1}(r,r)/(sum(diag(lme.CoefficientCovariance),'all')+lme.MSE+sum(diag(rp{1}),'all'));
    %     end

        % save
    %     save(lme_fname,'lme')

    end
    figure('name',[term ': PCV']); 
    rs=floor(sqrt(length(model_comp)));
    cs=ceil(length(model_comp)/rs);
    for pc = 1:length(model_comp)
        ax1(pc)=subplot(rs,cs,pc); boxplot(vertcat(model_comp(pc).PCVc{:})); title(['Pairs compared: ' num2str(S.paired_comp{pc})])
    end
    % linkaxes(ax1)

    figure('name',[term ': ICC,R2']); 
    a=0;
    rs=floor(sqrt(length(model_comp)*2));
    cs=ceil(length(model_comp)*2/rs);
    for li = 1:length(model)
        a=a+1; ax2(a)=subplot(length(model),2,a); boxplot(vertcat(model(li).ICC{:})); title(['Model ' num2str(li)])
        a=a+1; ax2(a)=subplot(length(model),2,a); boxplot(vertcat(model(li).R2{:})); title(['Model ' num2str(li)])
    end
    % linkaxes(ax2)
    % subplot(3,1,2); bar(vertcat(ICC{:})); title('ICC per model')
    % subplot(3,1,3); bar(vertcat(R2{:})); title('R2 per model')
end

function set_paths(S)
restoredefaultpath
for p = 1:size(S.path.code,1)
    if S.path.code{p,1}
        addpath(genpath(S.path.code{p,2}));
    else
        addpath(S.path.code{p,2});
    end
end
