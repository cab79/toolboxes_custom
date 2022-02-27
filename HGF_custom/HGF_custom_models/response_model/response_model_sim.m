function y = response_model_sim(r, infStates, pvec)
% Calculates the log-probability of log-reaction times y (in units of log-ms) according to the
% linear log-RT model developed with Louise Marshall and Sven Bestmann
% http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1002575
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2014-2016 Christoph Mathys, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.

% CAB: Number of levels
try
    l = r.c_prc.(r.c_obs.model).n_priorlevels+1;
catch
    l=1;
end

% Initialize returned log-probabilities, predictions,
% and residuals as NaNs so that NaN is returned for all
% irregualar trials
n = size(r.u,1);
y = nan(n,1);

% Create param struct
nme=r.c_obs.pnames;
nme_gen=r.c_obs.pnames_gen;
idx=r.c_obs.priormusi;

% Transform parameters to their native space
%pvec = logrt_softmax_binary_transp(r, pvec);

% Initialize random number generator
rng('shuffle');

% real or abs?
abs_switch = r.c_obs.abs;

%% SOFTMAX
if any(strcmp(r.c_obs.responses, 'Ch'))
    
    % Create param struct
    type='soft';
    for pn=1:length(nme)
        if strcmp(nme{pn,1}(1:length(type)),type)
            eval([nme_gen{pn} ' = pvec(idx{pn})'';']);
        end
    end
    
%     % Predictions or posteriors?
%     pop = 1; % Default: predictions
%     if r.c_obs.predorpost == 2
%         pop = 3; % Alternative: posteriors
%     end
% 
%     % Softmax parameter
%     be = p(8);

    % Softmax: Weed irregular trials out from inferred states, responses, and inputs
    x = r.traj.(r.c_obs.model).(r.c_obs.rep)(:,1);

    % Calculate log-probabilities for non-irregular trials
    % If input matrix has only one column, assume the weight (reward value)
    % of both options is equal to 1
    %if size(x,2) == 1
        if ~any(x<0) && ~any(x>1)
            % Apply the logistic sigmoid to the inferred states
            prob = tapas_sgm(be.*(2.*x-1),1);
        else
            error('tapas:hgf:SoftMaxBinary:InfStatesIncompatible', 'infStates incompatible with tapas_softmax_binary observation model.')
        end
    %end
    % If input matrix has three columns, the second contains the weights of
    % outcome 0 and the third contains the weights of outcome 1
%     if size(x,2) == 2 
%         % Apply the logistic sigmoid to the inferred states
%         prob = tapas_sgm(be.*(x(:,1)-x(:,2)),1);
%     end
    
    % Simulate
    y(:,1) = binornd(1, prob);
end
if any(strcmp(r.c_obs.responses, 'RT'))

    %% RT
    % Create param struct
    type='reg';
    for pn=1:length(nme)
        if strcmp(nme{pn,1}(1:length(type)),type)
            eval([nme_gen{pn} ' = pvec(idx{pn})'';']);
        end
    end


    u = r.u(:,1);

    if strcmp(r.c_obs.model,'like')
        % surprise-only model

        % Extract trajectories of interest
        mu1hat = r.traj.(r.c_obs.model).xchat(:,1);
        
        % Surprise: informational
        % ~~~~~~~~
        m1hreg = mu1hat;
        %m1hreg(r.irr) = [];
        poo = m1hreg.^u.*(1-m1hreg).^(1-u); % probability of observed outcome
        surp = -log2(poo);

%         % Bernoulli variance (aka irreducible uncertainty, risk) 
%         % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%         xcsahat = r.traj.(r.c_obs.model).xcsahat(:,1);
%         
%         bernvhat = xcsahat;
%         %bernvhat(r.irr) = [];
    
        % Calculate predicted log-response
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        logresp = be0 +be1.*surp;

    else
        % Extract trajectories of interest from infStates
        mu1 = r.traj.(r.c_obs.model).mu(:,1);
        mu1hat = r.traj.(r.c_obs.model).muhat(:,1);
        sa1hat = r.traj.(r.c_obs.model).sahat(:,1);
        sa1 = r.traj.(r.c_obs.model).sa(:,1);
        dau = r.traj.(r.c_obs.model).dau;
        ep1 = r.traj.(r.c_obs.model).epsi(:,1);
        da1 = r.traj.(r.c_obs.model).da(:,1);
        ep2 = r.traj.(r.c_obs.model).epsi(:,2);
        if l>1
            mu2    = r.traj.(r.c_obs.model).mu(:,2);
            sa2    = r.traj.(r.c_obs.model).sa(:,2);
            mu2hat = r.traj.(r.c_obs.model).muhat(:,2);
            sa2hat = r.traj.(r.c_obs.model).sahat(:,2);
        end
        if l>2
            mu3 = r.traj.(r.c_obs.model).mu(:,3);
            da2 = r.traj.(r.c_obs.model).da(:,2);
            ep3 = r.traj.(r.c_obs.model).epsi(:,3);
            da3 = r.traj.(r.c_obs.model).da(:,3);
        end
        
        
        % prediction error
        % ~~~~~~~~
        if abs_switch
            daureg = abs(dau);
            %daureg(r.irr) = [];
            ep1reg = abs(ep1);
            %ep1reg(r.irr) = [];
            da1reg = abs(da1);
            %da1reg(r.irr) = [];
            ep2reg = abs(ep2);
            %ep2reg(r.irr) = [];
            if l>2
                da2reg = abs(da2);
                %da2reg(r.irr) = [];
                ep3reg = abs(ep3);
                %ep3reg(r.irr) = [];
                da3reg = abs(da3);
                %da3reg(r.irr) = [];
            end
        else
            daureg = dau;
            %daureg(r.irr) = [];
            ep1reg = ep1;
            %ep1reg(r.irr) = [];
            da1reg = da1;
            %da1reg(r.irr) = [];
            ep2reg = ep2;
            %ep2reg(r.irr) = [];
            if l>2
                da2reg = da2;
                %da2reg(r.irr) = [];
                ep3reg = ep3;
                %ep3reg(r.irr) = [];
                da3reg = da3;
                %da3reg(r.irr) = [];
            end
        end
        
        % Posterior expectation
        % ~~~~~~~~
        m1reg = mu1;
        %m1reg(r.irr) = [];
        if l>1
            m2reg = mu2;
            %m2reg(r.irr) = [];
        end
        if l>2
            m3reg = mu3;
            %m3reg(r.irr) = [];
        end
    
        % Surprise: informational
        % ~~~~~~~~
        m1hreg = mu1hat;
        %m1hreg(r.irr) = [];
        poo = m1hreg.^u.*(1-m1hreg).^(1-u); % probability of observed outcome
        surp = -log2(poo);
    
        % Bernoulli variance (aka irreducible uncertainty, risk) 
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        bernv = sa1;
        %bernv(r.irr) = [];
        
        bernvhat = sa1hat;
        %bernvhat(r.irr) = [];
        
        mu1hreg = mu1hat;
        %mu1hreg(r.irr) = [];
    
        if l>1 % CAB
            % Inferential variance (aka informational or estimation uncertainty, ambiguity)
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %inferv = tapas_sgm(mu2, 1).*(1 -tapas_sgm(mu2, 1)).*sa2; % transform down to 1st level
            sigmoid_mu2 = 1./(1+exp(-mu2)); % transform down to 1st level
            inferv = sigmoid_mu2.*(1 -sigmoid_mu2).*sa2; 
            %inferv = sa2; 
            %inferv(r.irr) = [];
            
            mu2hreg = mu2hat;
            %mu2hreg(r.irr) = [];
            sa2hreg = sa2hat;
            %sa2hreg(r.irr) = [];
        end
    
        if l>2 % CAB
            % Phasic volatility (aka environmental or unexpected uncertainty)
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            pv = sigmoid_mu2.*(1-sigmoid_mu2).*exp(mu3); % transform down to 1st level
            %pv = exp(mu3); 
            %pv(r.irr) = [];
            
        end
        
    
        % Calculate predicted log-response
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if l>2
            logresp = be0 +be1.*surp +be2.*bernv +be3.*inferv +be4.*pv +be5.*daureg +be6.*ep1reg +be7.*da1reg +be8.*ep2reg +be9.*da2reg +be10.*ep3reg +be11.*da3reg +be12.*m1reg +be13.*m2reg +be14.*m3reg +be15.*sa2hreg +be16.*mu2hreg;
        elseif l>1
            logresp = be0 +be1.*surp +be2.*bernv +be3.*inferv +be5.*daureg +be6.*ep1reg +be7.*da1reg +be8.*ep2reg +be12.*m1reg +be13.*m2reg +be15.*sa2hreg +be16.*mu2hreg;
        else
            logresp = be0 +be1.*surp +be2.*bernv +be5.*daureg +be6.*ep1reg +be7.*da1reg +be8.*ep2reg +be12.*m1reg;
        end
    end

    
    % Simulate
    %y(:,2) = logresp +sqrt(ze)*randn(n, 1); % response time plus Gaussian noise
    y(:,2) = logresp; % response time without Gaussian noise
end


%% Regression
if any(strcmp(r.c_obs.responses, 'HGFvar'))

    % Create param struct
    type='reg'; % regression
    for pn=1:length(nme)
        if strcmp(nme{pn,1}(1:length(type)),type)
            eval([nme_gen{pn} ' = pvec(idx{pn})'';']);
        end
    end

%     u = r.u(:,1);
    
    ynames = r.c_obs.ynames;
    for yi=1:length(ynames)
        switch ynames{yi}
            case 'HGF_PL_epsi_1'
                ep1 = r.traj.(r.c_obs.model).epsi(:,1);
                logresp = be06 +be6.*ep1;
            case 'HGF_PL_epsi_2'
                ep2 = r.traj.(r.c_obs.model).epsi(:,2);
                logresp = be08 +be8.*ep2;
            case 'HGF_PL_epsi_3'
                ep3 = r.traj.(r.c_obs.model).epsi(:,3);
                logresp = be010 +be10.*ep3;
            case 'HGF_PL_epsi-abs_1'
                ep1 = abs(r.traj.(r.c_obs.model).epsi(:,1));
                logresp = be06 +be6.*ep1;
            case 'HGF_PL_epsi-abs_2'
                ep2 = abs(r.traj.(r.c_obs.model).epsi(:,2));
                logresp = be08 +be8.*ep2;
            case 'HGF_PL_epsi-abs_3'
                ep3 = abs(r.traj.(r.c_obs.model).epsi(:,3));
                logresp = be010 +be10.*ep3;
            case 'HGF_PL_mu_1'
                mu1 = r.traj.(r.c_obs.model).mu(:,1);
                logresp = be012 +be12.*mu1;
            case 'HGF_PL_mu_2'
                mu2 = r.traj.(r.c_obs.model).mu(:,2);
                logresp = be013 +be13.*mu2;
            case 'HGF_PL_mu_3'
                mu3 = r.traj.(r.c_obs.model).mu(:,3);
                logresp = be014 +be14.*mu3;
            case 'HGF_PL_sa_1'
                sa1 = r.traj.(r.c_obs.model).sa(:,1);
                logresp = be02 +be2.*sa1;
            case 'HGF_PL_inferv_1'
                mu2 = r.traj.(r.c_obs.model).mu(:,2);
                sa2 = r.traj.(r.c_obs.model).sa(:,2);
                sigmoid_mu2 = 1./(1+exp(-mu2)); % transform down to 1st level
                inferv = sigmoid_mu2.*(1 -sigmoid_mu2).*sa2; 
                logresp = be03 +be3.*inferv;
            case 'HGF_PL_pv_1'
                mu2 = r.traj.(r.c_obs.model).mu(:,2);
                mu3 = r.traj.(r.c_obs.model).mu(:,3);
                sigmoid_mu2 = 1./(1+exp(-mu2)); % transform down to 1st level
                pv = sigmoid_mu2.*(1-sigmoid_mu2).*exp(mu3); 
                logresp = be04 +be4.*pv;
        end
        y(:,2+yi) = logresp; % response time without Gaussian noise
    end
end

%% Regression with CCA components
if any(strcmp(r.c_obs.responses, 'CCA_FApred'))

    % Create param struct
    type='reg'; % regression
    for pn=1:length(nme)
        if strcmp(nme{pn,1}(1:length(type)),type)
            eval([nme_gen{pn} ' = pvec(idx{pn})'';']);
        end
    end

    u = r.u(:,1);
    
    % Extract trajectories of interest from infStates
    mu1 = r.traj.(r.c_obs.model).mu(:,1);
    mu1hat = r.traj.(r.c_obs.model).muhat(:,1);
    sa1hat = r.traj.(r.c_obs.model).sahat(:,1);
    sa1 = r.traj.(r.c_obs.model).sa(:,1);
    dau = r.traj.(r.c_obs.model).dau;
    ep1 = r.traj.(r.c_obs.model).epsi(:,1);
    da1 = r.traj.(r.c_obs.model).da(:,1);
    ep2 = r.traj.(r.c_obs.model).epsi(:,2);
    if l>1
        mu2    = r.traj.(r.c_obs.model).mu(:,2);
        sa2    = r.traj.(r.c_obs.model).sa(:,2);
        mu2hat = r.traj.(r.c_obs.model).muhat(:,2);
        sa2hat = r.traj.(r.c_obs.model).sahat(:,2);
    end
    if l>2
        mu3 = r.traj.(r.c_obs.model).mu(:,3);
        sa3 = r.traj.(r.c_obs.model).sa(:,3);
        da2 = r.traj.(r.c_obs.model).da(:,2);
        ep3 = r.traj.(r.c_obs.model).epsi(:,3);
        da3 = r.traj.(r.c_obs.model).da(:,3);
    end
    
    
%     % prediction error
%     % ~~~~~~~~
%     if abs_switch
%         daureg = abs(dau);
%         %daureg(r.irr) = [];
%         ep1reg = abs(ep1);
%         %ep1reg(r.irr) = [];
%         da1reg = abs(da1);
%         %da1reg(r.irr) = [];
%         ep2reg = abs(ep2);
%         %ep2reg(r.irr) = [];
%         if l>2
%             da2reg = abs(da2);
%             %da2reg(r.irr) = [];
%             ep3reg = abs(ep3);
%             %ep3reg(r.irr) = [];
%             da3reg = abs(da3);
%             %da3reg(r.irr) = [];
%         end
%     else
        daureg = dau;
        %daureg(r.irr) = [];
        ep1reg = ep1;
        %ep1reg(r.irr) = [];
        da1reg = da1;
        %da1reg(r.irr) = [];
        ep2reg = ep2;
        %ep2reg(r.irr) = [];
        if l>2
            da2reg = da2;
            %da2reg(r.irr) = [];
            ep3reg = ep3;
            %ep3reg(r.irr) = [];
            da3reg = da3;
            %da3reg(r.irr) = [];
        end
%     end
    
    % Posterior expectation
    % ~~~~~~~~
    m1reg = mu1;
    %m1reg(r.irr) = [];
    sa1reg = sa1;
    %sa1reg(r.irr) = [];
    if l>1
        m2reg = mu2;
        %m2reg(r.irr) = [];
        sa2reg = sa2;
        %sa2reg(r.irr) = [];
    end
    if l>2
        m3reg = mu3;
        %m3reg(r.irr) = [];
        sa3reg = sa3;
        %sa3reg(r.irr) = [];
    end

    % Surprise: informational
    % ~~~~~~~~
    %m1hreg = mu1hat;
    %m1hreg(r.irr) = [];
    poo = m1reg.^u.*(1-m1reg).^(1-u); % probability of observed outcome
    surp = -log2(poo);

    % Bernoulli variance (aka irreducible uncertainty, risk) 
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bernv = sa1;
    %bernv(r.irr) = [];
    
    bernvhat = sa1hat;
    %bernvhat(r.irr) = [];
    
    mu1hreg = mu1hat;
    %mu1hreg(r.irr) = [];

    if l>1 % CAB
        % Inferential variance (aka informational or estimation uncertainty, ambiguity)
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %inferv = tapas_sgm(mu2, 1).*(1 -tapas_sgm(mu2, 1)).*sa2; % transform down to 1st level
        sigmoid_mu2 = 1./(1+exp(-mu2)); % transform down to 1st level
        inferv = sigmoid_mu2.*(1 -sigmoid_mu2).*sa2; 
        %inferv = sa2; 
        %inferv(r.irr) = [];
        
        mu2hreg = mu2hat;
        %mu2hreg(r.irr) = [];
        sa2hreg = sa2hat;
        %sa2hreg(r.irr) = [];
    end

    if l>2 % CAB
        % Phasic volatility (aka environmental or unexpected uncertainty)
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pv = sigmoid_mu2.*(1-sigmoid_mu2).*exp(mu3); % transform down to 1st level
        %pv = exp(mu3); 
        %pv(r.irr) = [];
        
    end
    
    ynames = r.c_obs.ynames;
    weights = r.c_obs.weights;
    factors = r.c_obs.factors;
    fac = unique([factors{:}]);
    for f=fac
       yvar = ynames;%(facnum==f);
       yval=0;
       for yi = 1:length(yvar)
           wn = ismember(ynames,yvar{yi});
           switch yvar{yi}
                case 'HGF_PL_dau_1'
                    var = daureg;
                case 'HGF_PL_da_1'
                    var = da1reg;
                case 'HGF_PL_da_2'
                    var = da2reg;
                    if l<2; weights(wn,f)=0; end
                case 'HGF_PL_mu_1'
                    var = m1reg;
                case 'HGF_PL_mu_2'
                    var = m2reg;
                    if l<2; weights(wn,f)=0; end
                case 'HGF_PL_mu_3'
                    var = m3reg;
                    if l<3; weights(wn,f)=0; end
                case 'HGF_PL_sa_1'
                    var = sa1reg;
                case 'HGF_PL_sa_2'
                    var = sa2reg;
                    if l<2; weights(wn,f)=0; end
                case 'HGF_PL_sa_3'
                    var = sa3reg;
                    if l<3; weights(wn,f)=0; end
           end
           yval = yval + var*weights(wn,f);
        end
        eval(['f' num2str(f) ' = yval;']);
    end
        
    for yi=1:length(r.c_obs.params_cell)

        % assign intercept
        eval(['be0 = be0' num2str(yi) ';']);
        
        % assign slopes, one per factor (x4)
        param_nums = r.c_obs.params_cell{yi};
        eval(['be1 = be' num2str(param_nums(1)) ';']);
        eval(['be2 = be' num2str(param_nums(2)) ';']);
        eval(['be3 = be' num2str(param_nums(3)) ';']);
        eval(['be4 = be' num2str(param_nums(4)) ';']);
        
        % assign zeta
        eval(['ze = ze' num2str(yi) ';']);
        
        % Calculate predicted log-response
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if l>2
            logresp = be0 +be1.*f1 +be2.*f2 +be3.*f3 +be4.*f4;
        elseif l>1
            logresp = be0 +be1.*f1 +be2.*f2 +be3.*f3 +be4.*f4;
        else
            logresp = be0 +be1.*f1 +be2.*f2 +be3.*f3;
        end
        
        y(:,2+yi) = logresp; % response time without Gaussian noise
    end


end



return;
