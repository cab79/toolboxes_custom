function plot_eegstats_LMM_charts(S)
dbstop if error
%% raw data plots with cov on X axis: fitted data will not be linear trend lines due to other IV effects
% P1: input: fitted+residuals - scattered data points
% P1: fit: population trend line from fitted data
% P2: fit+random: individual subject trend lines.
% P3: inputfit_p: population trend lines from new OLS.
% P4: inputfit_s: subject trend lines from new OLS.
% P5: cov mean vs. slope from new OLS
%% partial regression plots: partial cov (after partialling out over IVs) on X axis. fitted data will be linear trend lines
% resid_nocov: input after partialing out other IV effects.
% fit: population trend lines from fitted data
% fit+random: individual subject trend lines, either from fitted data or coefficients
% residfit_p: population trend lines from new OLS.
% residfit_s: subject trend lines from new OLS.
% P5: partial cov mean vs. slope from fit+random
%% partial residual plots: cov on X axis. fitted data will not be linear trend lines
% resid_nocov: input after partialing out other IV effects
% fit: population trend lines from fitted data 
% fit+random: individual subject trend lines.
% inputfit_p: population trend lines from new OLS.
% inputfit_s: subject trend lines from new OLS.
% P5: cov mean vs. slope from fit+random
%% CODE
% close all
% set paths and file names 
S.path.code = {
   1, 'G:\Q_backup\MATLAB\toolboxes_external\cbrewer'
   1, 'G:\Q_backup\MATLAB\toolboxes_external\spm12'
    };

% set paths
set_paths(S)

if ~isempty(S.model.clus)
    % load cluster data
    load(fullfile(S.path.stats_load,'D.mat'),'D');
    pth=S.path.stats_load;
elseif ~isempty(S.model.component)
    % load component data
    load(S.path.stats_load,'D');
    [pth] = fileparts(S.path.stats_load);
end
cd(pth)

if ~isempty(S.model.def)
    D.model(S.model.index).def = S.model.def{:};
end

disp('loading data table')
Ddtab = load(S.path.dtab_inputs);
dtab = Ddtab.D.prep.dtab;

% find con
i = S.model.index;
try D.model(S.model.index)
catch
    error(['this model number does not exist'])
end
if ~any(strcmp({D.model(S.model.index).con(:).term},S.model.contrast_term))
    error([S.model.contrast_term ' is not a constrast term in this model'])
end
c = find(strcmp({D.model(S.model.index).con(:).term},S.model.contrast_term));



% ensure variables are in the right order - useful for extracting
% coefficients 
vnames=dtab.Properties.VariableNames;
vars = [S.IV_plot {S.cov_plot} {S.effect}]; % correct order
vars = fliplr(unique(fliplr(vars),'stable'));
for iv = 1:length(vars)
    if ~isempty(vars{iv}) && ~isempty(find(strcmp(vnames,vars{iv})))
        ivi(iv) = find(strcmp(vnames,vars{iv}));
    end
end
ivi(ivi==0)=[];
new_vnames = vnames;
new_vnames(ivi) = vnames(sort(ivi));
dtab = dtab(:,new_vnames);

% clusters from topo or components from CCA?
if ~isempty(S.model.clus)
    S.model.clus(S.model.clus>length(D.model(i).con(c).clus))=[];
    splot = S.model.clus;
elseif ~isempty(S.model.component)
    S.model.component(S.model.component>size(Ddtab.D.prep.PCA.CCA{1},2))=[];
    splot = S.model.component;
end

% subplots per cluster
for cl = splot
        
%     load_exist=questdlg('load existing models?');
    if ~isempty(S.model.clus)
        switch S.data_type
            case 'median'
                input=D.model(i).con(c).clus(cl).input_median'; % better to use mean?
            case 'mean'
                input=D.model(i).con(c).clus(cl).input_mean'; % better to use mean?
            case 'eig'
                input=D.model(i).con(c).clus(cl).eig(1,:)'; % better to use mean?
        end
        dtab.data=input;
    elseif ~isempty(S.model.component)
        dtab.data=Ddtab.D.prep.Y(cl).dtab.data;
    end
    
    
    
    lme_fname = fullfile(pth,['LME_model' num2str(S.model.index) '_cluster' num2str(cl) '.mat']);
    fit_models=1;
%     if exist(lme_fname,'file') && strcmp(load_exist,'Yes')
%         load(lme_fname);
%         fit_models=0;
%     end
    if fit_models
        
        % estimate LMM of control model to obtain residuals
        if isfield(S.model,'def_control') && ~isempty(S.model.def_control)
            disp(['estimating control LMM, cluster ' num2str(cl)])
            lme_control=fitlme(dtab,S.model.def_control{:},'FitMethod',S.fitmethod,'DummyVarCoding',S.coding,'CheckHessian', true);
            dtab.data = residuals(lme_control,'Conditional',true,'ResidualType','Raw');     
        end
        
        % re-estimate LMM to obtain fixed and random effect coefficients
        % and fitted responses
        disp(['re-estimating LMM, cluster ' num2str(cl)])
        lme=fitlme(dtab,D.model(i).def,'FitMethod',S.fitmethod,'DummyVarCoding',S.coding,'CheckHessian', true);
        
        disp(lme)
        anova(lme)

        % LMM without S.effects to obtain residuals
%         def=regexprep(D.model(i).def, '\s+', ''); % remove spaces
%         newdef = regexprep(def, ['[+*]\w*:?' S.effect '(:\w*)?'], ''); % remove main effects and interactions involving S.effect
        newdef=regexprep(S.model.def_reduced{:}, '\s+', '');
        disp(['estimating ' newdef ', cluster ' num2str(cl)])
        lme_nocovY=fitlme(dtab,newdef,'FitMethod',S.fitmethod,'DummyVarCoding',S.coding,'CheckHessian', true);

        disp(lme_nocovY)
        
        % LMM of cov vs. other IVs
        if strcmp(S.cov_plot,S.effect)
            newdef2 = regexprep(newdef, 'data', S.effect);
            disp(['estimating LMM ' newdef2 ', cluster ' num2str(cl)])
            lme_nocovX=fitlme(dtab,newdef2,'FitMethod',S.fitmethod,'DummyVarCoding',S.coding,'CheckHessian', true);
        else
            lme_nocovX=[];
        end
        
        % calculate R2 for model and for fixed effects for lme model only
        rp=covarianceParameters(lme);
        R2.fix.lme.total = sum(diag(lme.CoefficientCovariance),'all')/(sum(diag(lme.CoefficientCovariance),'all')+lme.MSE);
        R2.fixran.lme.total = (sum(diag(lme.CoefficientCovariance),'all')+sum(diag(rp{1}),'all'))/(sum(diag(lme.CoefficientCovariance),'all')+lme.MSE+sum(diag(rp{1}),'all'));
        for f = 1:length(lme.CoefficientCovariance)
            R2.fix.lme.fixed(f) = lme.CoefficientCovariance(f,f)/(sum(diag(lme.CoefficientCovariance),'all')+lme.MSE);
            R2.fixran.lme.fixed(f) = lme.CoefficientCovariance(f,f)/(sum(diag(lme.CoefficientCovariance),'all')+lme.MSE+sum(diag(rp{1}),'all'));
        end
        for r = 1:length(rp{1})
            R2.fix.lme.random(r) = rp{1}(r,r)/(sum(diag(lme.CoefficientCovariance),'all')+lme.MSE);
            R2.fixran.lme.random(r) = rp{1}(r,r)/(sum(diag(lme.CoefficientCovariance),'all')+lme.MSE+sum(diag(rp{1}),'all'));
        end
        
        % save
        save(lme_fname,'lme','lme_nocovX','lme_nocovY','R2')
    end

    resid = residuals(lme,'Conditional',true,'ResidualType','Raw');
    fit_conditional = fitted(lme,'Conditional',true);
    fit_marginal = fitted(lme,'Conditional',false);

    resid_nocovY = residuals(lme_nocovY,'Conditional',true,'ResidualType','Raw');
    fit_conditional_nocovY = fitted(lme_nocovY,'Conditional',true);
    fit_marginal_nocovY = fitted(lme_nocovY,'Conditional',false);

    if ~isempty(lme_nocovX)
        resid_nocovX = residuals(lme_nocovX,'Conditional',true,'ResidualType','Raw');
        fit_conditional_nocovX = fitted(lme_nocovX,'Conditional',true);
        fit_marginal_nocovX = fitted(lme_nocovX,'Conditional',false);
    else
        resid_nocovX = [];
        fit_conditional_nocovX = [];
        fit_marginal_nocovX = [];
    end

    % common variables
    disp(['Extracting data, cluster ' num2str(cl)])
    rdm = dtab.(S.random_name{:});
    if ~iscell(rdm) % temporary fix - can remove later
        for rd = 1:size(rdm,1)
            rdnew{rd,1} = rdm(rd,:);
        end
        rdm=rdnew;
    end
    [rdm_levels,rdm_idx] = unique(rdm,'stable');
    try
        cov = dtab.(S.cov_plot);
    catch
        cov=[];
    end
    try
        group = dtab.group(rdm_idx);
    catch
        cov=group;
    end

    % NEW: fitted values for plotting NOT YET IMPLEMENTED
    try % does not work unless there are categorical variables
        [fitM_val,fitC_val] = fitted_values(lme,fit_marginal,fit_conditional);
    end
    
    % IV_plot cells and labels FOR CATEGORICAL ERROR AND BOX PLOTS ONLY
    if isfield(S,'IV_plot_levels') && ~isempty(S.IV_plot_levels) && ~strcmp(S.IVplot_type,'line')
        levels=S.IV_plot_levels;
    else
        levels=[];
    end
    [cells_IVplot, celltypes_IVplot, IV_levels, label] = IV_cells_labels(dtab,S.IV_plot,levels);
    if any(ismember(S.IV_plot,S.effect))
        IV_noeffect = S.IV_plot; IV_noeffect(find(ismember(S.IV_plot,S.effect)))=[]; 
        levels_noeffect=levels; 
        levels_effect=levels; 
        if ~isempty(levels)
            levels_noeffect(find(ismember(S.IV_plot,S.effect)))=[];
            levels_effect(find(~ismember(S.IV_plot,S.effect)))=[];
        end
        [cells_IVplot_noeffect, celltypes_IVplot_noeffect, IV_levels_noeffect, label_noeffect] = IV_cells_labels(dtab,IV_noeffect,levels_noeffect);
        [cells_IVplot_effect, celltypes_IVplot_effect, IV_levels_effect, label_effect] = IV_cells_labels(dtab,{S.effect},levels_effect);
    end
    
    % coefficients for the effect of interest
    IV_plot = S.IV_plot; %intersect(S.IV_plot,lme.CoefficientNames);
    [coeff,~,Rcoeff] = FReffects(lme,S.effect,rdm_levels,IV_plot);
    %coeff_table = table(rdm_levels,group,coeff,Rcoeff);
    coeff_fname = fullfile(pth,['Coeff_model' num2str(S.model.index) '_cluster' num2str(cl) '_' S.effect '.mat']);
    save(coeff_fname,'coeff','Rcoeff')
   

    plotdata(1).title = {'observed EEG data','all effects','all effects'};
    plotdata(1).Y = dtab.data;
    plotdata(1).X_name = 'cov';
    plotdata(1).X = cov;
    plotdata(1).fit_marginal = fit_marginal;
    plotdata(1).fit_conditional = fit_conditional;
    plotdata(1).slope_type = 'ols'; % ols or fit
    plotdata(1).cellind = {cells_IVplot,cells_IVplot,cells_IVplot};
    plotdata(1).celltypes = {celltypes_IVplot,celltypes_IVplot,celltypes_IVplot};
    plotdata(1).celllabel = {label,label,label};
    plotdata(2).title = {['residual effect of ' S.effect],['without ' S.effect],['without ' S.effect]};
    plotdata(2).Y = resid_nocovY;
    plotdata(2).X_name = 'cov';
    plotdata(2).X = cov;
    plotdata(2).fit_marginal = fit_marginal_nocovY;
    plotdata(2).fit_conditional = fit_conditional_nocovY;
    plotdata(2).slope_type = 'fit'; % ols or fit
    plotdata(2).cellind = {cells_IVplot,cells_IVplot,cells_IVplot};
    plotdata(2).celltypes = {celltypes_IVplot,celltypes_IVplot,celltypes_IVplot};
    plotdata(2).celllabel = {label,label,label};
    plotdata(3).title = {['residual effect of ' S.effect],['without ' S.effect],['without ' S.effect]};
    plotdata(3).Y = resid_nocovY;
    plotdata(3).X_name = {'IV residuals:',['effect of ' S.effect]};
    plotdata(3).X = resid_nocovX;
    plotdata(3).fit_marginal = fit_marginal_nocovY;
    plotdata(3).fit_conditional = fit_conditional_nocovY;
    plotdata(3).slope_type = 'fit'; % ols or fit
    if any(ismember(S.IV_plot,S.effect))
        plotdata(3).cellind = {cells_IVplot_effect,cells_IVplot_noeffect,cells_IVplot_noeffect};
        plotdata(3).celltypes = {celltypes_IVplot_effect,celltypes_IVplot_noeffect,celltypes_IVplot_noeffect};
        plotdata(3).celllabel = {label_effect,label_noeffect,label_noeffect};
    else
        plotdata(3).cellind = {};
        plotdata(3).celltypes = {};
        plotdata(3).celllabel = {};
    end
    
    %% COVARIATE PLOTS
    if ~isempty(S.cov_plot)
        
        % initiate plotting 
%         col = cbrewer('seq', 'YlGnBu', length(IV_levels)+2, 'pchip'); col = col(3:end,:);
        if isfield(S,'IV_plot_col') && ~isempty(S.IV_plot_col)
            col = S.IV_plot_col{1};
        else
            col = cbrewer('qual', 'Set1', length(IV_levels), 'pchip');
        end
        
        figure('Name',['covariate: cluster ' num2str(cl)]);
        pi(1) = length(plotdata); pi(2) = 6; % 5x3 grid

        % to obtain mean estimate of the covariate, identify cells
        if ~isempty(S.cov_mean)
            cells_meanest = get_cells(dtab, dtab.Properties.VariableNames, S.cov_mean);
        else
            cells_meanest = ones(length(cov),1);
        end

        for pd = 1:length(plotdata)
            if isempty(plotdata(pd).X) || isempty(plotdata(pd).Y)
                continue
            end
            xSlice = [min(plotdata(pd).X) : range(plotdata(pd).X)/100 : max(plotdata(pd).X)];

            % fit based on LME coefficients for population
            if pd>1 
                IV_plot = S.IV_plot; %intersect(S.IV_plot,lme_nocovY.CoefficientNames);
                [coeff_plot,Fcoeff_plot] = FReffects(lme_nocovY,S.cov_plot,rdm_levels,IV_plot);
            else
                IV_plot = S.IV_plot; %intersect(S.IV_plot,lme.CoefficientNames);
                [coeff_plot,Fcoeff_plot] = FReffects(lme,S.cov_plot,rdm_levels,IV_plot);
            end
            if ~isempty(Fcoeff_plot)
                cf=0;
                e4iv = [];
                for g = 1:height(Fcoeff_plot)
                    for e=1:size(Fcoeff_plot.intercept,2)
                        cf=cf+1;
                        e4iv(cf) = e;
                        cfit{cf} = polyval([Fcoeff_plot.slope(g,e),Fcoeff_plot.intercept(g,e)], xSlice);
                    end
                end
            else
                cfit = {};
            end

            for iv=1:length(IV_levels) % e.g. groups
                if iscell(IV_levels)
                    IV_idx = strcmp(cells_IVplot,IV_levels{iv});
                else
                    IV_idx = cells_IVplot==IV_levels(iv);
                end
                if ~IV_idx
                    IV_idx = ~isnan(plotdata(pd).X);
                end
                [~,sort_idx] = sort(plotdata(pd).X(IV_idx));

                % linear fit coefficients for population
                [IVfitcoef(iv,:)] = robustfit(plotdata(pd).X(IV_idx),plotdata(pd).Y(IV_idx),S.regress);
                IVfit(iv,:) = polyval(IVfitcoef(iv,[2,1]), xSlice);

                % calculate slopes per random variable (e.g. subject)
                ri=0;
                sp=(pd-1)*pi(2);
                X={};Rsort_idx={};cRfit=[];
                for r=1:length(rdm_levels) % subjects
                    % get data
                    idx = IV_idx & strcmp(rdm,rdm_levels{r});
                    if ~any(idx)
                        continue
                    end
                    ri=ri+1;
                    X{ri}=plotdata(pd).X(idx);
                    Y{ri}=plotdata(pd).Y(idx);
%                         yFitM{ri}=plotdata(pd).fit_marginal(idx);
                    yFitC{ri}=plotdata(pd).fit_conditional(idx);
                    [~,Rsort_idx{ri}] = sort(X{ri});
% 
%                     if ~isempty(Rcoeff_plot)
%                         % LME fit random efficients in correct order
%                         Rcoeff_ordered(ri,:) = Rcoeff_plot(strcmp(rdm_levels{r}, Rnames),:);
%                         % fit based on LME coefficients for population
%                         coeff(ri,:) = Fcoeff_plot + Rcoeff_ordered(ri,:);
%                     else
%                         coeff(ri,:) = Fcoeff_plot;
%                     end
                    if ~isempty(coeff_plot)
                        cRfit(ri,:) = polyval([coeff_plot.slope(r,e4iv(iv)),coeff_plot.intercept(r,e4iv(iv))], xSlice);
                    end

                    % linear fit coefficients and y values per trial
                    [Rfitcoef(ri,:)] = robustfit(X{ri},Y{ri},S.regress);
                    Rfit(ri,:) = polyval(Rfitcoef(ri,[2,1]), xSlice);

                    % produces confidence intervals for each subjects' fit
    %                 [coef(ri,:),stats] = polyFitFit(x,Y,1);
    %                 [yFitFit(ri,:),yCI(ri,:)] = polyconf(coef(ri,[2,1]),xSlice,stats, 'predopt','curve');


                    % obtain mean estimate of the covariate
                    cu=unique(cells_meanest(idx));
                    cMX=[];
                    for ci = 1:numel(cu)
                        cMX(ci,1) = mean(X{ri}(cu==cu(ci)));
                    end
                    cMean(ri,1) = mean(cMX);

                end



                disp(['Plotting data, cluster ' num2str(cl) ', data type ' num2str(pd) ', IV level ' num2str(iv)])

                % scatter plot of observations: population/marginal fit
                sp=sp+1; subplot(pi(1),pi(2),sp); title(['obs & popn fit']); hold on
                xlabel(plotdata(pd).X_name)
%                 ylabel()
                scatter(plotdata(pd).X(IV_idx),plotdata(pd).Y(IV_idx),10,col(iv,:)/2);
                % trend line: sort by cov
                find_idx = find(IV_idx);
                p(iv)=plot(plotdata(pd).X(find_idx(sort_idx)),plotdata(pd).fit_marginal(find_idx(sort_idx)),'LineWidth',3,'Color',col(iv,:));

                % scatter plot of observations: population/marginal fit
                sp=sp+1; subplot(pi(1),pi(2),sp); title(['popn and coeff fit']); hold on
                xlabel(plotdata(pd).X_name)
%                 ylabel()
                for cf = 1:length(cfit)
                    plot(xSlice,cfit{cf},'k--'); 
                end
                find_idx = find(IV_idx);
                plot(plotdata(pd).X(find_idx(sort_idx)),plotdata(pd).fit_marginal(find_idx(sort_idx)),'LineWidth',2,'Color',col(iv,:));

                % fitted response per random variable (e.g. subject)
                sp=sp+1; subplot(pi(1),pi(2),sp); title(['fit by ' S.random_name{:}]); hold on
                xlabel(plotdata(pd).X_name)
%                 ylabel()
                for ri = 1:length(X)
                    plot(X{ri}(Rsort_idx{ri}),yFitC{ri}(Rsort_idx{ri}),'Color',col(iv,:))
                end

                % scatter plot of observations with OLS fit
                sp=sp+1; subplot(pi(1),pi(2),sp); title(['obs & OLS fit']); hold on
                xlabel(plotdata(pd).X_name)
%                 ylabel()
                scatter(plotdata(pd).X(IV_idx),plotdata(pd).Y(IV_idx),10,col(iv,:)/2);
                % trend line
                plot(xSlice,IVfit(iv,:),'LineWidth',3,'Color',col(iv,:))

                % fitted response per random variable (e.g. subject) with OLS fit
                sp=sp+1; ax(1)=subplot(pi(1),pi(2),sp); title(['OLS fit by ' S.random_name{:}]); hold on
                xlabel(plotdata(pd).X_name)
%                 ylabel()
                for ri = 1:length(X)
                    plot(xSlice,Rfit(ri,:),'Color',col(iv,:))
                end

                % plot for each group
%                 upper = yMean+dsp;
%                 lower = yMean-dsp;
%                 sp=sp+1; subplot(pi(1),pi(2),sp); title(['mean+' S.dispersion]); hold on
%                 fill([xSlice, fliplr(xSlice)], [(upper), fliplr((lower))], ...
%                 col(iv,:), 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
%                 plot(xSlice,yMean,col(iv,:)); 
%                 xlabel('cov')
%                 ylabel('fit')

                % fit lines from coefficients
                sp=sp+1; ax(2)=subplot(pi(1),pi(2),sp); title(['coefficient lines']); hold on
                xlabel(plotdata(pd).X_name)
%                 ylabel()
                if ~isempty(cRfit)
                    for ri = 1:length(X)
                        plot(xSlice,cRfit(ri,:),'Color',col(iv,:))
                    end
                end
                if ~isempty(cfit)
                    plot(xSlice,cfit{iv},'k--','LineWidth',2); 
                end
                y1=get(ax(1),'YLim');
                y2=get(ax(2),'YLim');
                linkaxes([ax(1) ax(2)])
                ylim = y2*1.5;
                set(ax(1),'YLim',ylim);

%                 % cov mean vs. slope
%                 sp=sp+1; subplot(pi(1),pi(2),sp); title(['cov mean vs. ' plotdata(pd).slope_type ' slope']); hold on
%                 if strcmp(plotdata(pd).slope_type,'fit')
%                     coef = repmat(Fcoeff,length(X),1) + Rcoeff_ordered;
%                     scatter(cMean,coef(:,2),10,col(iv,:)); 
%                 elseif strcmp(plotdata(pd).slope_type,'ols')
%                     scatter(cMean,Rfitcoef(:,2),10,col(iv,:)); 
%                 end
%                 xlabel('cov mean')
%                 ylabel('slope')

            end
            if ~isempty(label) && pd==1
                legend(p,label)
            end
        end

        hold off
    end

    %% IV plots
    if ~isempty(celltypes_IVplot)
        figure('Name',['IVs: cluster ' num2str(cl)]);
        pi(1) = length(plotdata); pi(2) = 3; 

        % 1. raw data
        % 2. popn fit (fixed)
        % 3. ind fit with error bars (fixed + random)

        for pd = 1:length(plotdata)
            sp=(pd-1)*pi(2);

            switch S.IVplot_type 
                case 'box'
                    col = cbrewer('seq', 'YlGnBu', length(IV_levels)+2, 'pchip'); col = col(3:end,:);
                    sp=sp+1; subplot(pi(1),pi(2),sp); title(plotdata(pd).title{1}); ylabel(['amplitude']); hold on
                    boxplot(plotdata(pd).Y,plotdata(pd).cellind{1},'labels',plotdata(pd).celllabel{1},'LabelOrientation','inline','colors',col,'plotstyle','compact');
                    sp=sp+1; subplot(pi(1),pi(2),sp); title(plotdata(pd).title{2}); ylabel(['marginal: fixed']); hold on
                    boxplot(plotdata(pd).fit_marginal,plotdata(pd).cellind{2},'labels',plotdata(pd).celllabel{2},'LabelOrientation','inline','colors',col,'plotstyle','compact');
                    sp=sp+1; subplot(pi(1),pi(2),sp); title(plotdata(pd).title{3}); ylabel(['conditional: fixed+random']); hold on
                    boxplot(plotdata(pd).fit_conditional,plotdata(pd).cellind{3},'labels',plotdata(pd).celllabel{3},'LabelOrientation','inline','colors',col,'plotstyle','compact');
                case {'error','error_grouped'}
                    if ~isempty(plotdata(pd).celllabel) && ~isempty(plotdata(pd).celllabel{1})
                        sp=sp+1; ax(1)=subplot(pi(1),pi(2),sp); title(plotdata(pd).title{1}); ylabel(['amplitude, mean(CI)']); hold on
                        if strcmp(S.IVplot_type,'error_grouped') && pd<3
                            col = cbrewer('seq', 'YlGnBu', length(IV_levels_noeffect), 'pchip'); col = col(1:end,:);
                            errorplot_grouped(plotdata(3).cellind{3},plotdata(pd).Y,S.dispersion,plotdata(3).celllabel{3},col,plotdata(3).celllabel{1},plotdata(3).cellind{1})
                        elseif strcmp(S.IVplot_type,'error_grouped') && pd==3 % use different colors
                            col = cbrewer('qual', 'Set1', length(IV_levels), 'pchip');
                            errorplot(plotdata(pd).cellind{1},plotdata(pd).Y,S.dispersion,plotdata(pd).celllabel{1},col)
                        else
                            col = cbrewer('seq', 'YlGnBu', length(IV_levels), 'pchip'); col = col(1:end,:);
                            errorplot(plotdata(pd).cellind{1},plotdata(pd).Y,S.dispersion,plotdata(pd).celllabel{1},col)
                        end
                    end
                    if ~isempty(plotdata(pd).celllabel) && ~isempty(plotdata(pd).celllabel{2})
                        sp=sp+1; ax(2)=subplot(pi(1),pi(2),sp); title(plotdata(pd).title{2}); ylabel(['marginal: fixed']); hold on
                        col = cbrewer('seq', 'YlGnBu', length(IV_levels)+2, 'pchip'); col = col(3:end,:);
                        errorplot(plotdata(pd).cellind{2},plotdata(pd).fit_marginal,'',plotdata(pd).celllabel{2},col)
                    end
                    if ~isempty(plotdata(pd).celllabel) && ~isempty(plotdata(pd).celllabel{3})
                        sp=sp+1; ax(3)=subplot(pi(1),pi(2),sp); title(plotdata(pd).title{3}); ylabel(['conditional: fixed+random']); hold on
                        col = cbrewer('seq', 'YlGnBu', length(IV_levels)+2, 'pchip'); col = col(3:end,:);
                        errorplot(plotdata(pd).cellind{3},plotdata(pd).fit_conditional,'',plotdata(pd).celllabel{3},col)
                    end
                    [ax(:).ActivePositionProperty] = deal('position');
                
                case {'line'}
                    if ~isempty(plotdata(pd).celllabel) && ~isempty(plotdata(pd).celllabel{1})
                        sp=sp+1; ax(1)=subplot(pi(1),pi(2),sp); title(plotdata(pd).title{1}); ylabel(['amplitude, mean(CI)']); hold on
                        if isfield(S,'IV_plot_col') && ~isempty(S.IV_plot_col)
                            col = S.IV_plot_col;
                        else
                            col = cbrewer('seq', 'YlGnBu', length(IV_levels_noeffect), 'pchip'); col = col(1:end,:);
                        end
                        lineplot(plotdata(pd).celltypes{1},dtab,plotdata(pd).Y,S.dispersion,col,S.IVplot_axes,S.IV_plot_levels)
                    end
                    if ~isempty(plotdata(pd).celllabel) && ~isempty(plotdata(pd).celllabel{2})
                        sp=sp+1; ax(2)=subplot(pi(1),pi(2),sp); title(plotdata(pd).title{2}); ylabel(['marginal: fixed']); hold on
                        if isfield(S,'IV_plot_col') && ~isempty(S.IV_plot_col)
                            col = S.IV_plot_col;
                        else
                            col = cbrewer('seq', 'YlGnBu', length(IV_levels)+2, 'pchip'); col = col(3:end,:);
                        end
                        lineplot(plotdata(pd).celltypes{2},dtab,plotdata(pd).fit_marginal,'',col,S.IVplot_axes,S.IV_plot_levels)
                    
                    end
                    if ~isempty(plotdata(pd).celllabel) && ~isempty(plotdata(pd).celllabel{3})
                        sp=sp+1; ax(3)=subplot(pi(1),pi(2),sp); title(plotdata(pd).title{3}); ylabel(['conditional: fixed+random']); hold on
                        if isfield(S,'IV_plot_col') && ~isempty(S.IV_plot_col)
                            col = S.IV_plot_col;
                        else
                            col = cbrewer('seq', 'YlGnBu', length(IV_levels)+2, 'pchip'); col = col(3:end,:);
                        end
                        lineplot(plotdata(pd).celltypes{3},dtab,plotdata(pd).fit_conditional,'',col,S.IVplot_axes,S.IV_plot_levels)
                    end
                    [ax(:).ActivePositionProperty] = deal('position');
            end
        end
    end
end


   

function set_paths(S)
for p = 1:size(S.path.code,1)
    if S.path.code{p,1}
        addpath(genpath(S.path.code{p,2}));
    else
        addpath(S.path.code{p,2});
    end
end

function [cells_IVplot, celltypes_IVplot, IV_levels, label] = IV_cells_labels(design,IVs,levels)
if ~isempty(IVs)
    [cells_IVplot, celltypes_IVplot] = get_cells(design, design.Properties.VariableNames, IVs);
    if ~isempty(levels)
        newcelltypes_IVplot = table;
        for lev=1:length(levels)
            if iscategorical(celltypes_IVplot.(IVs{lev})) && isprotected(celltypes_IVplot.(IVs{lev}))
                celltypes_IVplot.(IVs{lev}) = cellstr(celltypes_IVplot.(IVs{lev}));
            end
            oldlev=unique(celltypes_IVplot.(IVs{lev}),'stable');
            for ol = 1:length(oldlev)
                if iscell(oldlev)
                    newcelltypes_IVplot.(IVs{lev})(strcmp(celltypes_IVplot.(IVs{lev}),oldlev{ol}))=levels{lev}(ol);
                else
                    newcelltypes_IVplot.(IVs{lev})(celltypes_IVplot.(IVs{lev})==oldlev(ol))=levels{lev}(ol);
                end
            end
            celltypes_IVplot.(IVs{lev}) = cellstr(newcelltypes_IVplot.(IVs{lev}));
        end
    end
else
    cells_IVplot = [];
    celltypes_IVplot = [];
end
if ~isempty(celltypes_IVplot)
    IV = cells_IVplot;
    IV_levels = unique(IV,'stable'); 
else
    IV_levels = {'IV'}; 
    IV=repmat(IV_levels,height(design),1);
end

% cell labels
if ~isempty(celltypes_IVplot)
    nCells = size(celltypes_IVplot,1); 
    for ce = 1:nCells
        label{ce}='';
        for iv = 1:length(IVs)
            temp=celltypes_IVplot.(IVs{iv})(ce);
            if ~isempty(levels)
                label{ce} = [label{ce} ' ' temp{:}];
            else
                if ~iscell(temp)
                    temp=num2str(double(temp));
                    label{ce} = [label{ce} ' ' IVs{iv} '_' temp];
                else
                    temp=temp{:};
                    label{ce} = [label{ce} ' ' IVs{iv} ':' temp];
                end
            end
        end
    end
else
    label = {};
end

function [cells,celltypes] = get_cells(design,headernames,IVname)
for iv = 1:length(IVname)
    header_col(iv) = find(ismember(headernames,IVname(iv)));
end
temp=design(:,header_col);
[celltypes,~,cells] = unique(temp,'rows','stable');


function errorplot(cells_IVplot,plotdata,dispersion,label,col)
nCells=length(label);
for ce = 1:nCells
    idx=cells_IVplot==ce;
    meanY(ce)=mean(plotdata(idx));
    ySD(ce)=std(plotdata(idx));
    ySE(ce) = ySD(ce)/sqrt(sum(idx));
    tscore = tinv(0.025,sum(idx)-1);
    yCI(ce) = tscore*ySE(ce);
end

% dispersion
switch dispersion
    case 'CI'
        dsp = yCI;
    case 'SE'
        dsp = ySE;
    case 'SD'
        dsp = ySD;
    case ''
        dsp = [];
end

Xcat=categorical(label);
hold on
b=bar(Xcat,meanY);
b.FaceColor = 'flat';
for k = 1:length(meanY)
    b.CData(k,:) = col(k,:);
end
if any(dsp~=0)
    errorbar(Xcat,meanY,dsp,'.k')
end
hold off


function errorplot_grouped(cells_IVplot,plotdata,dispersion,label,col,grouping_label,grouping_cells)

Xcat=categorical(grouping_label);
meanY=[];
dsp=[];
for g=1:length(grouping_label)
    nCells=length(label);
    for ce = 1:nCells
        idx=cells_IVplot==ce & grouping_cells==g;
        temp_meanY(ce)=mean(plotdata(idx));
        ySD(ce)=std(plotdata(idx));
        ySE(ce) = ySD(ce)/sqrt(sum(idx));
        tscore = tinv(0.025,sum(idx)-1);
        yCI(ce) = tscore*ySE(ce);
    end

    % dispersion
    switch dispersion
        case 'CI'
            temp_dsp = yCI;
        case 'SE'
            temp_dsp = ySE;
        case 'SD'
            temp_dsp = ySD;
        case ''
            temp_dsp = [];
    end
    meanY=[meanY;temp_meanY];
    dsp=[dsp;temp_dsp];
end

hold on
clear b
b=bar(meanY);
pause(0.1)
if any(dsp~=0,'all')
    for k1 = 1:length(b)
        ctr(k1,:) = bsxfun(@plus, b(1).XData, [b(k1).XOffset]');
        ydt(k1,:) = b(k1).YData;
    end
    errorbar(ctr', ydt', dsp, '.k')
%     errorbar(Xcat,meanY,dsp,'.')
end
set(gca, 'XTick',1:length(Xcat),'XTickLabel', Xcat);
for b1=1:length(b)
    b(b1).FaceColor = 'flat';
    for k = 1:size(meanY,2)
        b(b1).CData(k,:) = col(b1,:);
    end
end
legend(label)
hold off


function lineplot(design,fulldesign,plotdata,dispersion,col,IVplot_axes,labels)

IV_plot = design.Properties.VariableNames;
fulldesign = fulldesign(:,ismember(fulldesign.Properties.VariableNames,IV_plot));
IVplot_axes = IVplot_axes(end-length(IV_plot)+1:end);
plot_axes = unique(IVplot_axes,'stable');
labels = labels(end-length(plot_axes)+1:end);

% subplot definitions
if ismember('subplot',plot_axes)
    splot_name = IV_plot{ismember(IVplot_axes,'subplot')};
    splot_levels = unique(design.(splot_name));
    splot_labels = labels{ismember(plot_axes,'subplot')};
    splots = 1:length(splot_levels);
else
    splot_name = '';
    splot_levels = '';
    splots = 1;
end

% lines
if ismember('line',plot_axes)
    line_name = IV_plot(ismember(IVplot_axes,'line'));
    line_levels = unique(design(:,ismember(design.Properties.VariableNames,line_name)));
    line_labels = labels{ismember(plot_axes,'line')};
else
    line_name = '';
    line_levels = {};
    line_labels = {};
end

% X axis
Xcat_name = IV_plot{ismember(IVplot_axes,'X')};
Xcat_levels = unique(design.(Xcat_name));
Xcat_labels = labels{ismember(plot_axes,'X')};

for g=splots
    meanY=[];
    dsp=[];
    for x = 1:length(Xcat_levels)
        
        % select design for this group/Xcat
        if ismember('subplot',plot_axes)
            design_sub = design(ismember(design.(splot_name),splot_levels(g)),:);
        else
            design_sub = design;
        end
        design_sub = design_sub(ismember(design_sub.(Xcat_name),Xcat_levels(x)),:);
        design_sub = unique(design_sub);
        
        % get mean and dispersion values from this group and X value
        for ce = 1:height(design_sub)
            
            % get column indices for 'line' factors
            idx = find(ismember(fulldesign,design_sub(ce,:),'rows'));

            % select data for this group/Xcat
            temp_meanY(1,ce)=mean(plotdata(idx));
            ySD(1,ce)=std(plotdata(idx));
            ySE(1,ce) = ySD(1,ce)/sqrt(length(idx));
            tscore = tinv(0.025,sum(idx)-1);
            yCI(1,ce) = tscore*ySE(ce);
        end

        % dispersion
        switch dispersion
            case 'CI'
                temp_dsp = yCI;
            case 'SE'
                temp_dsp = ySE;
            case 'SD'
                temp_dsp = ySD;
            case ''
                temp_dsp = [];
        end
        meanY=[meanY;temp_meanY];
        dsp=[dsp;temp_dsp];
    end
    hold all
    for i=1:size(meanY,2)
        xpos{g} = [(g-1)+0.25 g-0.25];
        line(xpos{g},meanY(:,i)','LineWidth',1,'LineStyle','--','Color',col(i,:));
        if isempty(dsp)
            scat(i)=scatter(xpos{g},meanY(:,i)',100,col(i,:),'filled','s');
        else
            scat(i)=errorbar(xpos{g},meanY(:,i)',dsp(:,i)','--s','MarkerSize',1,'LineWidth',1,'Color',col(i,:),'MarkerEdgeColor',col(i,:),'MarkerFaceColor',col(i,:));
        end
    end
end
yl=ylim;
ylim([yl(1)-diff(yl)*0.1,yl(2)+diff(yl)*0.1])
xlim([0 g])
xticks(horzcat(xpos{:}))
xticklabels(repmat(Xcat_labels,1,g));
if ismember('subplot',plot_axes)
    for g = splots
        text(mean(xpos{g}),yl(2),splot_labels{g},'HorizontalAlignment','center');
    end
end
% set(gca,'fontsize',fontsize)

% legend
if ismember('line',plot_axes)
    leg_order = [1:length(line_labels)];
    legend(scat(leg_order),line_labels(leg_order))
end

function [fitM_val,fitC_val] = fitted_values(lme,fitM,fitC)
% EXTRACT FITTED VALUES PER CELL/SUBJECT FOR PLOTTING

% define IVs
IVs = lme.PredictorNames;
catV=[];
contV=[];
for v = 1:length(IVs)
    dat=lme.Variables.(IVs{v});
    if strcmp(IVs{v},'ID')
        idV = v;
    elseif iscategorical(dat) || iscell(dat)
        catV=[catV v];
    else
        contV=[contV v];
    end
end
% obtain marginal fitted values per cell
[cells, fitM_val] = IV_cells_labels(lme.Variables,IVs(catV),[]);
[u,~,irows] = unique(cells);
for i = u'
    Ufitted(i,1)=mean(fitM(irows==i));
end
fitM_val.fitted=Ufitted;
% obtain conditional fitted values per cell/subject
[cells, fitC_val] = IV_cells_labels(lme.Variables,IVs([idV catV]),[]);
[u,~,irows] = unique(cells);
for i = u'
    Cfitted(i,1)=mean(fitC(irows==i));
end
fitC_val.fitted=Cfitted;

function [coeff,Fcoeff,Rcoeff] = FReffects(lme,effectstr,rdm_levels,IVs)
% EXTRACT COEFFICIENTS

effects=strsplit(effectstr,':');
% get effects
[F,Fn]=fixedEffects(lme);
Fn=Fn.Name;
[R,Rn]=randomEffects(lme);
Rnames = Rn.Level(strcmp(Rn.Name,'(Intercept)'));
[~,Ridx] = ismember(rdm_levels, Rnames);
Rn=Rn.Name;
% is group a factor?
if any(strcmp(IVs,'group'))
    grprange=lme.VariableInfo('group','Range');
    grps=grprange{:,:};
    grps=grps{:};
    gcell=[{''} grps(2:end)];
else
    gcell={''};
    grps={''};
end
% build within-subject terms
wIVs = IVs; 
wIVs(strcmp(wIVs,'group'))=[];
effect_is_IV = find(ismember(wIVs,effects));
if ~isempty(effect_is_IV)
    wIVs(effect_is_IV)=[];
    effectstr = [effectstr '_1'];
end
if ~isempty(wIVs)
    [~, celltypes] = get_cells(lme.Variables, lme.Variables.Properties.VariableNames, wIVs);
    terms={};term={};
    for c = 1:height(celltypes)
        for w = 1:length(wIVs)
            if double(celltypes.(wIVs{w})(c))>1
                terms{c,w} = [wIVs{w} '_' num2str(double(celltypes.(wIVs{w})(c))-1) ':'];
            else
                terms{c,w} = '';
            end
        end
        term{c,1}=horzcat(terms{c,:});
    end
    term = regexprep(term,'^:','');
    term = regexprep(term,':$','');
else
    term={''};
end

% initialise tables
Rcoeff=table;
Fcoeff=table;
% for each group
for g = 1:length(gcell)

    % get terms for this group
    if g==1
        eff_exp = '^((?!group).)*$'; % must not contain group term
    else
        eff_exp = ['.*group_' gcell{g} '.*']; % must contain group term
    end
    eff_fun = @(s)~cellfun('isempty',regexp(Fn,eff_exp));
    eff_out = cellfun(eff_fun,effects,'UniformOutput',false);
    eff_idx = all(horzcat(eff_out{:}),2);
    Fg=F(eff_idx);
    Fng=Fn(eff_idx);
    
    % remove group term to make later processing simpler
    Fng= regexprep(Fng,['(:?)group_' gcell{g} '(:?)'],'');
    
    % current version: random and fixed effects must correspond (i.e. full random slopes and intercepts model)
    if g==1
        
        if ~isequal(Fng,unique(Rn,'stable'))
%             error('random and fixed effect terms do not correspond')
            
        end
    end
    
    % get term indices and intercepts
    if g==1
        Fng{strcmp(Fng,'(Intercept)')}='';
    end
    for t = 1:size(term,1)
        term_idx(t) = find(strcmp(Fng,term(t)));
    end
    
    % identify the slopes of interest
%     eff_exp = '.*'; % expr = no colon before the first effect term
%     for eff = 1:length(effects)
%         eff_exp = [eff_exp effects{eff} '.*']; % expr = no colon after last effect term
%     end
%     eff_fun = @(s)~cellfun('isempty',regexp(Fng,eff_exp));
%     eff_out = cellfun(eff_fun,effects,'UniformOutput',false);
%     eff_idx = find(all(horzcat(eff_out{:}),2));
%     eff_terms = Fng(eff_idx);
    
    % build effect terms)
    eff_terms={};
    slopes=1;
    for t = 1:size(term,1)
        eff_terms{t,1} = [term{t} ':' effectstr];
        eff_terms{t,1} = regexprep(eff_terms{t,1},'^:','');
        if any(strcmp(Fng,eff_terms(t)))
            slp_idx(t) = find(strcmp(Fng,eff_terms(t)));
        else
            slopes = 0;
        end
    end
   
    % fixed effects
    Fx_int = [Fg(term_idx(1)), Fg(term_idx(2:end))'+Fg(term_idx(1))];
    Fcoeff.intercept(g,:)=Fx_int;
    if slopes
        Fx_slp = [Fg(slp_idx(1)), Fg(slp_idx(2:end))'+Fg(slp_idx(1))];
        Fcoeff.slope(g,:)=Fx_slp;
    else
        Fcoeff.slope(g,:)=zeros(size(Fcoeff.intercept(g,:)));
    end
    
    % if g>1 add to first group's (reference) data
    if g>1
%         Fcoeff{g,2:end} = Fcoeff{g,2:end} + Fcoeff{1,2:end}; % don't include intercept
        Fcoeff.intercept(g,:)=Fcoeff.intercept(g,:)+Fcoeff.intercept(1,:);
        Fcoeff.slope(g,:)=Fcoeff.slope(g,:)+Fcoeff.slope(1,:);
    end
end
% random effects
Rn(strcmp(Rn,'(Intercept)'))={''};
for t = 1:length(term)
    if t==1
        Rcoeff.intercept(:,t)=R(strcmp(Rn,term{t}));
        if slopes
            Rcoeff.slope(:,t)=R(strcmp(Rn,eff_terms{t}));
        end
    else
        if slopes
            Rcoeff.intercept(:,t)=R(strcmp(Rn,term{t})) + R(strcmp(Rn,term{1}));
            Rcoeff.slope(:,t)=R(strcmp(Rn,eff_terms{t})) + R(strcmp(Rn,eff_terms{1}));
        end
    end
end
if ~slopes
    Rcoeff.slope=zeros(size(Rcoeff.intercept));
end
Rcoeff=Rcoeff(Ridx,:);
% create coeff table
coeff=table;
coeff.rdm_levels=rdm_levels;
for r = 1:height(Rcoeff)
    coeff.group(r) = unique(lme.Variables.group(strcmp(lme.Variables.ID,rdm_levels(r))));
end
Rcoeff = [coeff,Rcoeff];
coeff=Rcoeff;
% add in fixed effects
if slopes
    for r = 1:height(coeff)
        gi = 1;
        if any(strcmp(coeff.group(r),gcell))
            gi=find(strcmp(coeff.group(r),gcell));
        end
        coeff{r,3:end} = coeff{r,3:end} + Fcoeff{gi,:};
    end
end
gtab=table; gtab.group=grps';
Fcoeff=[gtab Fcoeff];

