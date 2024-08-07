function varargout=HGF_group_model_comparison(S,varargin)
% OUTPUT: varargout is {posterior,out,gposterior,gout}
addpath('E:\Q_backup\MATLAB\toolboxes_external\cbrewer'); % (https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab)

% Set up display
scrsz = get(0,'screenSize');
outerpos = [0.05*scrsz(3),0.05*scrsz(4),0.95*scrsz(3),0.95*scrsz(4)];

if nargin>1 && ~isempty(varargin{1})
    LME = varargin{1};
    LME_input = 1;
else
    LME_input = 0; % load from HGF structures
end

if nargin>2 && ~isempty(varargin{2})
    LME_se = varargin{2};
    LME_se_input = 1;
else
    LME_se_input = 0; % load from HGF structures
end

if ~isfield(S,'family_on')
    S.family_on=0;
end

if ~iscell(S.fname_ext)
    S.fname_ext = {S.fname_ext};
end

if ~LME_input
    fitted_path = fullfile(S.path.hgf);
    cd(fitted_path)
end

% group info
if isfield(S,'grplist')
    grplist = S.grplist;
elseif isfield(S,'designmat')
    grplist = S.designmat(2:end,strcmp(S.designmat(1,:),'groups'));
elseif isfield(S,'designtab')
    grplist = S.designtab.groups;
end
grpuni = unique(grplist,'stable');

% for each model
mr=[];
rc=cell(1,1);

% testing perceptual models
mi=0;
if S.family_on
    pm_family{1,length(S.perc_models)} = [];
    rm_family{1,length(S.resp_models)} = [];
end

% axis labels
if isfield(S,'xlabel')
    yl = S.ylabel;
    xl = S.xlabel;
else
    yl = 'perceptual models';
    xl = 'response models';
end

% axis ticks
if isfield(S,'xticklabels')
    yt = S.yticklabels;
    xt = S.xticklabels;
else
    yt = S.perc_models;
    xt = S.resp_models;
end

% convert to strings
if isnumeric(S.resp_models)
    S.resp_models = strsplit(num2str(S.resp_models));
end
if isnumeric(S.perc_models)
    S.perc_models = strsplit(num2str(S.perc_models));
end

for pm = 1:length(S.perc_models)
%     load(fullfile(fitted_path,['CORE_fittedparameters_percmodel' num2str(S.perc_models(m)) '_respmodel' num2str(S.resp_model) '_fractrain' num2str(S.frac_train) '_' S.sname '.mat']));

    % testing response models
    for rm = 1:length(S.resp_models)
        
        % same model, different runs
        for om = 1:length(S.fname_ext)
            mi=mi+1;   

            if S.family_on
                % factor matrix and families
                fm(mi,1) = pm;
                fm(mi,2) = rm;
                pm_family{pm} = [pm_family{pm} mi];
                rm_family{rm} = [rm_family{rm} mi];
            end

            if ~LME_input
                try
                    fname=[S.fname_pref '_percmodel' S.perc_models{pm} '_respmodel' S.resp_models{rm} S.fname_ext{om}];
                    ls = load(fullfile(fitted_path,fname));
                catch
                    fname=[S.fname_pref '_pm' S.perc_models{pm} '_rm' S.resp_models{rm} S.fname_ext{om}];
                    ls = load(fullfile(fitted_path,fname));
                end
                D_fit=ls.D_fit;

                % select subjects (only those within designtab)
                D_fit(~ismember({D_fit(:).subname}, S.designtab.subjects)) = [];

                % mean evidence for each fitted model over repetitions
                for d = 1:length(D_fit)
                    try
                        LME(mi,d) = D_fit(d).HGF.fit.optim.LME; % rows: fitted models; columns: simulated models
                    catch
                        mr=[mr d];
                        try
                            disp(['Missing LME from file: ' fname])
                        end
                    end
                    if isfield(D_fit(d).HGF.fit,'eGBM')
                        D_fit(d).HGF.fit = rmfield(D_fit(d).HGF.fit,'eGBM');
                    end
                    rs(mi,d)=D_fit(d).HGF.fit;
                    rc{mi,d}=D_fit(d).HGF.fit;

                    % find number of subjects with unique priors (due to inabilty to optimise at the commonly use prior)
                    try
                        priors(mi,d,:) = D_fit(d).HGF.fit.c_prc.PL.ommu(2:3);
                    end
                end
            end
        end
        
    end
end

if ~LME_input
    % remove bad subjects
    mr = unique(mr);
    LME(:,mr) = [];
    if LME_se_input
        LME_se(:,mr) = [];
    end
    rs(:,mr) = [];
    rc(:,mr) = [];
    grplist(mr) = [];
    %priors(:,mr,:) = [];
end
    
for mi = 1:size(LME,1)
    
    % separate into groups
    for g = 1:length(grpuni)
        LMEgrp{1,g}(mi,:) = LME(mi,strcmp(grplist,grpuni{g}));
        if LME_se_input
            LMEgrp_se{1,g}(mi,:) = LME_se(mi,strcmp(grplist,grpuni{g}));
        end
    end

    if ~LME_input
        % diagnostics: parameter correlations
%         tapas_fit_plotCorr_group(rs(mi,:));
%         title(['model ' num2str(mi) ', simple average'])
        try 
            rc_out=tapas_bayesian_parameter_average_CAB(1,rc(mi,:));
%             tapas_fit_plotCorr(rc_out);
%             title(['model ' num2str(mi) ', bayesian average'])
        catch
            rc_out = [];
        end
    else
        rc_out = [];
    end
end
varargout={LMEgrp,rc_out};

% plot mean LME; separate by groups
if (rm>1 || pm>1) && om==1 
    
    % pm on x axis
%     pLME = reshape(sum(cat(2,LMEgrp{:}),2),rm,pm);
%     clims = [min(pLME(:)), max(pLME(:))];
%     figure
%     ax1=subplot('Position',[0.1 0.2 0.3 0.6])
%     imagesc(pLME,clims)
%     ax1.XAxisLocation = 'top';
%     xlabel('perceptual models')
%     ylabel('response models')
%     title(['Log-model evidence'])
%     ax2=subplot('Position',[0.42 0.2 0.1 0.6])
%     imagesc(mean(pLME,2),clims)
%     set(gca,'XTick',1,'XTickLabel','mean')
%     set(gca,'YTick',[])
% %     ax2.YAxisLocation = 'right';
%     ax2.XAxisLocation = 'top';
%     ax3=subplot('Position',[0.1 0.08 0.3 0.1])
%     imagesc(mean(pLME,1),clims)
%     set(gca,'YTick',1,'YTickLabel','mean')
%     set(gca,'XTick',[])
%     colorbar('Position',[0.54 0.08 0.03 0.72])
    
    pLME = reshape(sum(cat(2,LMEgrp{:}),2),rm,pm)';
    clims = [min(pLME(:)), max(pLME(:))];
    image_plot(pLME,S,clims,'Log-model evidence','mean',{},xl,yl,xt,yt,outerpos)
    
    % plot mean LME; separate by groups
%     for g = 1:length(grpuni)
%         % mean over subjects
%         pLME = reshape(sum(LMEgrp{1,g},2),rm,pm);
%         figure
%         imagesc(pLME)
%         colorbar
%         xlabel('perc models')
%         ylabel('resp models')
%         title(['sum LME, group ' num2str(g)])
%     end

    if LME_se_input
        % Assuming pLME is a 2D matrix of PSIS-LOO values and pLME_se is the corresponding standard error matrix

        % LMEgrp and LMEgrp_se are cell arrays where each cell contains a 2D matrix of values for different subjects
        
        % Summing LME and combining SE across subjects and reshaping into vectors
        pLME_vec = sum(cat(2, LMEgrp{:}), 2);
        pLME_se_vec = sqrt(sum(cat(2, LMEgrp_se{:}).^2, 2));
        
        % Number of models (total number of elements)
        num_models = length(pLME_vec);
        
        % Function to compute the significance of the difference
        is_significant = @(diff, se_diff) abs(diff) > 2 * se_diff;
        
        % Initialize matrices to store pairwise comparison results
        pairwise_diffs = zeros(num_models, num_models);
        pairwise_se_diffs = zeros(num_models, num_models);
        significant_wins = zeros(num_models, 1);
        
        % Pairwise comparison of models across all combinations
        for i = 1:num_models
            for j = 1:num_models
                if i ~= j
                    % Compute the difference in PSIS-LOO and its standard error
                    diff_loo = pLME_vec(i) - pLME_vec(j);
                    se_diff_loo = sqrt(pLME_se_vec(i)^2 + pLME_se_vec(j)^2);
        
                    % Store the differences and standard errors
                    pairwise_diffs(i, j) = diff_loo;
                    pairwise_se_diffs(i, j) = se_diff_loo;
        
                    % Check if the difference is significant
                    if is_significant(diff_loo, se_diff_loo)
                        % If model i has a significantly greater (negative) PSIS-LOO than model j, record a win
                        if diff_loo > 0
                            significant_wins(i) = significant_wins(i) + 1;
                        end
                    end
                end
            end
        end
        
        % Reshape significant_wins back to the original matrix shape
        significant_wins_matrix = reshape(significant_wins, rm, pm)';
        
        clims = [min(significant_wins_matrix(:)), max([1;significant_wins_matrix(:)])];
        image_plot(significant_wins_matrix,S,clims,'N sig. wins','mean',{},xl,yl,xt,yt,outerpos)
        
    end

end

% compare models
% if length(LMEgrp)==1
%     [posterior,out] = VBA_groupBMC(LME);
%     varargout(1:2) = {posterior,out};
% elseif length(LMEgrp)>1
%     [gposterior,gout] = VBA_groupBMC_btwGroups_CAB(LMEgrp)
%     varargout(3:4) = {gposterior,gout};
% end

options.DisplayWin = 0;
if S.family_on
    % compare perc model families
    options.families = pm_family;
    if length(LMEgrp)==1
        [~,pm_out] = VBA_groupBMC_cab(LME,options,0);
        pm_out = {pm_out};
    elseif length(LMEgrp)>1
        options.grpNames = grpuni;
        [~,~,pm_out] = VBA_groupBMC_btwGroups_CAB(LMEgrp,options,0);
    end

    % compare resp model families
    options.families = rm_family;
    if length(LMEgrp)==1
        [~,rm_out] = VBA_groupBMC_cab(LME,options,0);
        rm_out = {rm_out};
    elseif length(LMEgrp)>1
        [~,~,rm_out] = VBA_groupBMC_btwGroups_CAB(LMEgrp,options,0);
    end
    varargout=[varargout {pm_out,rm_out}];
    
    % plot Efs - model frequencies
    ef = pm_out{1}.Ef';
    ef_pm = pm_out{1}.families.Ef';
    ef_rm = rm_out{1}.families.Ef';
    ef = reshape(ef,rm,pm)';
    clims = [min(ef(:)), max(ef(:))];
    image_plot(ef,S,clims,'Estimated model frequencies','family',{ef_pm,ef_rm},xl,yl,xt,yt,outerpos)
    % figure
    % ax1=subplot('Position',[0.1 0.2 0.3 0.6])
    % imagesc(ef,clims)
    % ax1.XAxisLocation = 'top';
    % xlabel('perceptual models')
    % ylabel('response models')
    % title(['Estimated model frequencies'])
    % ax2=subplot('Position',[0.42 0.2 0.1 0.6])
    % imagesc(ef_rm',clims)
    % set(gca,'XTick',1,'XTickLabel','family')
    % set(gca,'YTick',[])
    % %     ax2.YAxisLocation = 'right';
    % ax2.XAxisLocation = 'top';
    % ax3=subplot('Position',[0.1 0.08 0.3 0.1])
    % imagesc(ef_pm,clims)
    % set(gca,'YTick',1,'YTickLabel','family')
    % set(gca,'XTick',[])
    % colorbar('Position',[0.54 0.08 0.03 0.72])

    % plot EPs
    ep = pm_out{1}.ep;
    ep_pm = pm_out{1}.families.ep;
    ep_rm = rm_out{1}.families.ep;
    ep = reshape(ep,rm,pm)';
    clims = [0 1];
    image_plot(ep,S,clims,'Exceedence probability (EP)','family',{ep_pm,ep_rm},xl,yl,xt,yt,outerpos)
    % figure
    % ax1=subplot('Position',[0.1 0.2 0.3 0.6])
    % imagesc(ep,clims)
    % ax1.XAxisLocation = 'top';
    % xlabel('perceptual models')
    % ylabel('response models')
    % title(['Exceedence probability (EP)'])
    % ax2=subplot('Position',[0.42 0.2 0.1 0.6])
    % imagesc(ep_rm',clims)
    % set(gca,'XTick',1,'XTickLabel','family')
    % set(gca,'YTick',[])
    % %     ax2.YAxisLocation = 'right';
    % ax2.XAxisLocation = 'top';
    % ax3=subplot('Position',[0.1 0.08 0.3 0.1])
    % imagesc(ep_pm,clims)
    % set(gca,'YTick',1,'YTickLabel','family')
    % set(gca,'XTick',[])
    % colorbar('Position',[0.54 0.08 0.03 0.72])

    % plot PEPs
    pep = pm_out{1}.pep;
    pep_pm = pm_out{1}.families.pep;
    pep_rm = rm_out{1}.families.pep;
    pep = reshape(pep,rm,pm)';
    clims = [0 1];
    image_plot(pep,S,clims,'Protected exceedence probability (PEP)','family',{pep_pm,pep_rm},xl,yl,xt,yt,outerpos)
    % figure
    % ax1=subplot('Position',[0.1 0.2 0.3 0.6])
    % imagesc(pep,clims)
    % ax1.XAxisLocation = 'top';
    % xlabel('perceptual models')
    % ylabel('response models')
    % title(['Protected exceedence probability (PEP)'])
    % ax2=subplot('Position',[0.42 0.2 0.1 0.6])
    % imagesc(pep_rm',clims)
    % set(gca,'XTick',1,'XTickLabel','family')
    % set(gca,'YTick',[])
    % %     ax2.YAxisLocation = 'right';
    % ax2.XAxisLocation = 'top';
    % ax3=subplot('Position',[0.1 0.08 0.3 0.1])
    % imagesc(pep_pm,clims)
    % set(gca,'YTick',1,'YTickLabel','family')
    % set(gca,'XTick',[])
    % colorbar('Position',[0.54 0.08 0.03 0.72])

    % plot num of unique priors
    % figure
    % for mi = 1:size(priors,3)
    %     for pi=1:2 % om priors
    %         subplot(size(priors,3))
    %         [~,~,ui]=length(unique(priors(mi,:,pi)));
    %     end
    % end
    
else
    options = {};
    % compare models
    if length(LMEgrp)==1
        [~,out] = VBA_groupBMC_cab(LME,options,S.pep_flag);
        out={out};
    elseif length(LMEgrp)>1
        options.grpNames = grpuni;
        [~,~,out] = VBA_groupBMC_btwGroups_CAB(LMEgrp,options,S.pep_flag);
    end
    if length(S.perc_models)>1
        pm_out=out;
        rm_out={};
    else
        rm_out=out;
        pm_out={};
    end
    varargout=[varargout {pm_out,rm_out}];
end

function image_plot(data,S,clims,title_text,XYlabel,family,xl,yl,xt,yt,outerpos)
if isempty(family)
    family{1} = mean(data,2)';
    family{2} = mean(data,1);
end
figure('OuterPosition', outerpos,'Name',title_text);
ax1=subplot('Position',[0.1 0.2 0.5 0.5]); % [left bottom width height]
imagesc(data,clims)
set(gca,'XTick',1:length(S.resp_models),'YTick',1:length(S.perc_models))
if isnumeric(xt)
    xt=strsplit(num2str(xt));
elseif iscell(xt{1})
    xt=horzcat(xt{:});
end
if isnumeric(yt)
    yt=strsplit(num2str(yt));
elseif iscell(yt{1})
    yt=horzcat(yt{:});
end
set(gca,'XTickLabel',xt,'YTickLabel',yt,'FontSize',10)
ax1.XAxisLocation = 'top';
ax1.YAxisLocation = 'left';
ylabel(yl)
xlabel(xl)
t=title(title_text,'position',[0 -0.8]);
set(t,'horizontalAlignment','left');
if size(data,1)>1 && size(data,2)>1
    % right family
    ax2=subplot('Position',[0.5+0.12 0.2 0.1 0.5]) % [left bottom width height] 
    imagesc(family{1}',clims)
    set(gca,'XTick',1,'XTickLabel',XYlabel,'FontSize',10)
    set(gca,'YTick',[])
    ax2.XAxisLocation = 'top';
    % lower family
    ax3=subplot('Position',[0.1 0.08 0.5 0.1])
    imagesc(family{2},clims)
    set(gca,'YTick',1,'YTickLabel',XYlabel,'FontSize',10)
    set(gca,'XTick',[])
    colorbar('Position',[0.5+0.24 0.08 0.03 0.5+0.12]) % [left bottom width height]
else
    colorbar('Position',[0.5+0.12 0.2 0.03 0.5+0.12]) % [left bottom width height]
end

colormap(flipud(cbrewer('seq', 'Reds', 100, 'pchip')))