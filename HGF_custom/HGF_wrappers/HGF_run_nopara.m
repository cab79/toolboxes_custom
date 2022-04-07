function D=HGF_run(D,S,sim_on)

% sim_on = 0: fit model to responses
% sim_on = 1: simulate responses
% sim_on = 2: fit to simulated responses

dbstop if error

if ~exist(S.path.hgf,'dir')
    mkdir(S.path.hgf)
end
if ~isfield(S.HGF,'plottraj')
    S.HGF.plottraj=0;
end

sdate = datestr(now,30);

    
%% HGF
% prc: perceptual; obs:observation; opt:optimisation
prc_model = S.prc_config;
pmodel = prc_model; %strrep(prc_model,'_config',''); % CAB updated - was causing problems for CORE sim_predict_EEG
obs_model = S.obs_config;
opt_algo = 'tapas_quasinewton_optim_config';

if ~S.bayes_opt
    S.use_y_col = find([any(strcmp(S.resp_modelspec.responses,'Ch')), any(strcmp(S.resp_modelspec.responses,'RT')),any(strcmp(S.resp_modelspec.responses,'EEG')) || any(strcmp(S.resp_modelspec.responses,'HGFvar')) || any(strcmp(S.resp_modelspec.responses,'CCA_FApred'))]);
end

% when simulating there may be multiple y created
if isfield(S.HGF,'selectrep')
    yn = S.HGF.selectrep;
else
    yn = 1;
end

if ~isfield(S,'parallel')
    S.parallel = 0;
end
if S.parallel
    checkp = gcp('nocreate');
    if isempty(checkp)
        myPool = parpool;
    end
    parforArg = Inf;
else
    parforArg = 0;
end

% If S contains D_perc, this is used to update the perceptual model
% parameters and fixed them (variance zero) per subject. Used for
% estimating response model parameters only for certain applications.
% Can easily be updated so that the priors are not fixed.
if isfield(S,'D_perc')
    r=struct;
    for d = 1:length(D)
        r(d).c_prc = eval([S.prc_config '(S,0)']); % runs the prc config function
        % prc: re-define within-subject priors
        r(d).c_prc.priormus = S.D_perc(d).HGF.fit.p_prc.p;
        r(d).c_prc.priorsas = S.D_perc(d).HGF.fit.c_prc.priorsas;
        % fixed variances to 0 for alpha and omega/theta
        r(d).c_prc.priorsas(S.fix_param_indices) = 0;
    end
end

% get parameter struct by running bayesopt, otherwise parallel fails
% fit = tapas_fitModel_CAB([], D(1).HGF(1).u, prc_model, S.bayesopt_config, opt_algo,S, 0);
% sim=fit;
fit=cell(length(D),1);sim=cell(length(D),1);
% parfor (d = 1:length(D),parforArg)
for d = 1:length(D)
    
    disp(['testing perc model ' num2str(S.perc_model) ', resp model ' num2str(S.resp_model) ', subject ' num2str(d) '/' num2str(length(D))])

    if ~isfield(S,'nstim')
        nst=length(D(d).HGF(1).u);
    else
        if isempty(S.nstim)
            nst=length(D(d).HGF(1).u);
        else
            nst=S.nstim;
        end
    end
    
    %cycle through until VarApprox is valid
    fin=0;
    failed=0;
    while fin==0
        
%         if isfield(S,'D_perc')
%             c_prc = r(d).c_prc;
%             if failed 
%                 c_prc.priormus(S.fix_param_indices(end-1)) = c_prc.priormus(S.fix_param_indices(end-1))-failed;
%             end
%         else
            c_prc = [];
%         end
        
        if sim_on==1
            if isfield(S,'sim') && ~isempty(S.sim)
                sim_param = S.sim;
                sim_obs_param=[];
            elseif isfield(D.HGF,'fit') 
                sim_param = D(d).HGF(yn).fit.p_prc;
                sim_obs_param = D(d).HGF(yn).fit.p_obs;
            end
            sim{d} = tapas_simModel_CAB(D(d).HGF(1).u(1:nst,:), pmodel, sim_param, obs_model, sim_obs_param, S, failed);
        elseif S.bayes_opt
            fit{d} = tapas_fitModel_CAB([], D(d).HGF(1).u(1:nst,:), prc_model, obs_model, opt_algo,S, failed); %BAYES OPTIMAL
        elseif strfind(S.prc_config,'tapas')
            fit{d} = tapas_fitModel(D(d).HGF(yn).y(1:nst,:), D(d).HGF(1).u(1:nst,:), prc_model, obs_model, opt_algo,S, failed);
            tapas_hgf_binary_plotTraj(fit{d});
        else
            if sim_on==2 && isfield(D(d).HGF(yn),'sim')
                y = D(d).HGF(yn).sim.y;
            elseif isfield(D(d).HGF(yn),'y')
                y = D(d).HGF(yn).y;
            end
            out = tapas_fitModel_CAB(y(1:nst,:), D(d).HGF(1).u(1:nst,:), prc_model, obs_model, opt_algo, S, failed, c_prc);
            if isfield(out,'traj')
                fit{d} = out;
            else
                failed = failed+1;
                continue
            end
        end

        if ~sim_on
            priormodels = fieldnames(S.perc_modelspec.priormodels);
            for pm = 1:length(priormodels)
                if isfield(fit{d}.traj,priormodels{pm})
                    ntrials = min(length(fit{d}.traj.(priormodels{pm}).mu),500);
                    dmu = diff(fit{d}.traj.(priormodels{pm}).mu(1:ntrials,2:end));
                    dpi = diff(1./(fit{d}.traj.(priormodels{pm}).sa(1:ntrials,2:end)));
                    rmdmu = repmat(sqrt(mean(dmu.^2)),length(dmu),1);
                    rmdpi = repmat(sqrt(mean(dpi.^2)),length(dpi),1);
                    
                    jumpTol = 16;
                    if any(abs(dmu(:)) > jumpTol*rmdmu(:)) || any(abs(dpi(:)) > jumpTol*rmdpi(:))
                        failed = failed+1;
                        disp('HGF_run: Variational approximation invalid. Parameters are in a region where model assumptions are violated.');
                        disp('Use plot for diagnosis: see within function'); % plot(abs(dpi(:))); hold on; plot(rmdpi(:),'r'); hold on; plot(jumpTol*rmdpi(:),'g')
                    else
                        fin=1;
                    end
                else
                    failed = failed+1;
                    fin=0;
                    disp('HGF_run: Variational approximation invalid. Parameters are in a region where model assumptions are violated.');
                    disp('Use plot for diagnosis: see within function'); % plot(abs(dpi(:))); hold on; plot(rmdpi(:),'r'); hold on; plot(jumpTol*rmdpi(:),'g')
                end
            end
        else
            fin=1;
        end
    end

    if S.parallel==0
        % assign parameters to D
        if sim_on==1
            D(d).HGF(yn).sim = sim{d};
        else
            D(d).HGF(yn).fit = fit{d};
        end
    end

    % PLOTS
    if sim_on==1
        plotstruct=sim{d};
        %sname = [S(1).select.subjects{1, d} '_' sdate '_sim.mat'];
        %save(fullfile(S.path.hgf,sname),'sim');
        %varargout={sim{d}};
    else
        plotstruct=fit{d};
        %sname = [S(1).select.subjects{1, d} '_' sdate '_fit.mat'];
        %save(fullfile(S.path.hgf,sname),'fit');
        %varargout={fit{d}};
    end
    if S.HGF.plottraj
        if isfield(plotstruct.traj,'AL')
            tapas_hgf_binary_plotTraj_CAB(plotstruct,'AL',true)
        end
        if isfield(plotstruct.traj,'PR')
            tapas_hgf_binary_plotTraj_CAB(plotstruct,'PR',false)
        end
        if isfield(plotstruct.traj,'PL')
            tapas_hgf_binary_plotTraj_CAB(plotstruct,'PL',true)
        end
    end
    if isfield(S.HGF,'plotdiag') && S.HGF.plotdiag
        tapas_fit_plotResidualDiagnostics(fit{d});
        tapas_fit_plotCorr(fit{d});
    end

    if S.parallel ==0 && sim_on==0 && S.bayes_opt==0
        save(fullfile(S.path.hgf,['temp_percmodel' num2str(S.perc_model) '_respmodel' num2str(S.resp_model) '_fractrain' num2str(S.frac_train) '_' S.sname '.mat']), 'D', 'S');
    end
end



% assign parameters to D
for d = 1:length(D)
    if sim_on==1
        D(d).HGF(yn).sim = sim{d};
    else
        D(d).HGF(yn).fit = fit{d};
    end
end

disp('FINISHED');
% LME
%LMEs can be used to calculate Bayes factors by exponentiating the difference in LME
%between two models applied to the same dataset. For example, an LME difference of 3
%implies a Bayes factor of about 20.
%For a fixed-effects analysis with several datasets (e.g., from different subjects), add up
%the LMEs for the different datasets and compare the LME sums.