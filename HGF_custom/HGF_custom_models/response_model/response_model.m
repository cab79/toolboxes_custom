function [logp, yhat, res] = response_model(r, infStates, pvec)
% Calculates the log-probability of log-reaction times y (in units of log-ms) according to the
% linear log-RT model developed with Louise Marshall and Sven Bestmann
% http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1002575

try
    l = r.c_prc.(r.c_obs.model).n_priorlevels+1;
catch
    l=1;
end

% Initialize returned log-probabilities, predictions,
% and residuals as NaNs so that NaN is returned for all
% irregualar trials
n = size(r.y,1);
reg = ~ismember(1:n,r.irr);
logp = NaN(n,1);
yhat = NaN(n,1);
res  = NaN(n,1);

% Transform parameters to their native space
pvec = response_model_transp(r, pvec);

% Create param struct
nme=r.c_obs.pnames;
nme_gen=r.c_obs.pnames_gen;
idx=r.c_obs.priormusi;

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
%     % Softmax parameter
%     be = exp(ptrans(8));

    ycol=1;

    % Softmax: Weed irregular trials out from inferred states, responses, and inputs
    x = r.traj.(r.c_obs.model).(r.c_obs.rep)(:,1);
    x(r.irr) = [];
    ys = r.y(:,ycol);
    ys(r.irr) = [];

    if size(r.u,2) == 3
        r0 = r.u(:,2);
        r0(r.irr) = [];
        r1 = r.u(:,3);
        r1(r.irr) = [];
    end

    % Calculate log-probabilities for non-irregular trials
    % If input matrix has only one column, assume the weight (reward value)
    % of both options is equal to 1
    if size(r.u,2) == 2 % CAB changed to 2 from 1, as second column now contains conditions
        % Probability of observed choice
        probc = 1./(1+exp(-be.*(2.*x-1).*(2.*ys-1)));
    end
    % If input matrix has three columns, the second contains the weights of
    % outcome 0 and the third contains the weights of outcome 1
    if size(r.u,2) == 3 % CAB changed to 3 from 2, as second column now contains conditions
        % Probability of observed choice
        probc = 1./(1+exp(-be.*(r1.*x-r0.*(1-x)).*(2.*ys-1)));
    end
    logp_so(reg) = log(probc);
    yh = ys.*probc +(1-ys).*(1-probc);
    yhat(reg) = yh;
    res(reg) = (ys -yh)./sqrt(yh.*(1 -yh));
end


%% Regression
if any(strcmp(r.c_obs.responses, 'RT')) || any(strcmp(r.c_obs.responses, 'EEG'))

    % Create param struct
    type='reg'; % regression
    for pn=1:length(nme)
        if strcmp(nme{pn,1}(1:length(type)),type)
            eval([nme_gen{pn} ' = pvec(idx{pn})'';']);
        end
    end
    
    % switch the data input column
    if any(strcmp(r.c_obs.responses, 'RT'))
        ycol=2;
    elseif any(strcmp(r.c_obs.responses, 'EEG'))
        ycol=3:size(r.y,2);
    end

    % Weed irregular trials out from responses and inputs
    yr = r.y(:,ycol); % RTs are in column 2, EEG in column 3 or greater
    yr(r.irr,:) = [];
    u = r.u(:,1);
    u(r.irr) = [];
    
    % check inputs are logged
%     if any(yr>2)
%         error('inputs must be logged and not contain any extreme outliers')
%     end

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
        daureg(r.irr) = [];
        ep1reg = abs(ep1);
        ep1reg(r.irr) = [];
        da1reg = abs(da1);
        da1reg(r.irr) = [];
        ep2reg = abs(ep2);
        ep2reg(r.irr) = [];
        if l>2
            da2reg = abs(da2);
            da2reg(r.irr) = [];
            ep3reg = abs(ep3);
            ep3reg(r.irr) = [];
            da3reg = abs(da3);
            da3reg(r.irr) = [];
        end
    else
        daureg = dau;
        daureg(r.irr) = [];
        ep1reg = ep1;
        ep1reg(r.irr) = [];
        da1reg = da1;
        da1reg(r.irr) = [];
        ep2reg = ep2;
        ep2reg(r.irr) = [];
        if l>2
            da2reg = da2;
            da2reg(r.irr) = [];
            ep3reg = ep3;
            ep3reg(r.irr) = [];
            da3reg = da3;
            da3reg(r.irr) = [];
        end
    end
    
    % Posterior expectation
    % ~~~~~~~~
    m1reg = mu1;
    m1reg(r.irr) = [];
    if l>1
        m2reg = mu2;
        m2reg(r.irr) = [];
    end
    if l>2
        m3reg = mu3;
        m3reg(r.irr) = [];
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
    bernv(r.irr) = [];
    
    bernvhat = sa1hat;
    bernvhat(r.irr) = [];
    
    mu1hreg = mu1hat;
    mu1hreg(r.irr) = [];

    if l>1 % CAB
        % Inferential variance (aka informational or estimation uncertainty, ambiguity)
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %inferv = tapas_sgm(mu2, 1).*(1 -tapas_sgm(mu2, 1)).*sa2; % transform down to 1st level
        sigmoid_mu2 = 1./(1+exp(-mu2)); % transform down to 1st level
        inferv = sigmoid_mu2.*(1 -sigmoid_mu2).*sa2; 
        %inferv = sa2; 
        inferv(r.irr) = [];
        
        mu2hreg = mu2hat;
        mu2hreg(r.irr) = [];
        sa2hreg = sa2hat;
        sa2hreg(r.irr) = [];
    end

    if l>2 % CAB
        % Phasic volatility (aka environmental or unexpected uncertainty)
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pv = sigmoid_mu2.*(1-sigmoid_mu2).*exp(mu3); % transform down to 1st level
        %pv = exp(mu3); 
        pv(r.irr) = [];
        
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

    % Calculate log-probabilities for non-irregular trials
    % Note: 8*atan(1) == 2*pi (this is used to guard against
    % errors resulting from having used pi as a variable).
    logp_reg(reg) = -1/2.*log(8*atan(1).*ze) -(yr-logresp).^2./(2.*ze);
    yhat(reg) = logresp;
    res(reg) = yr-logresp;
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
    
    % switch the data input column
    ycol = 3:size(r.y,2);
    if ~isempty(r.c_obs.ynum)
        ycol = ycol(r.c_obs.ynum);
    end

    % Weed irregular trials out from responses and inputs
    yr = r.y(:,ycol); % RTs are in column 2, EEG in column 3 or greater
%     yr = r.y(:,2); % HACK: NEED TO CHANGE BACK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    yr(r.irr,:) = [];
    u = r.u(:,1);
    u(r.irr) = [];
    
    logp_reg(reg,1) = 0;
    ynames = r.c_obs.ynames;
    for y=1:length(ynames)
        switch ynames{y}
            case 'HGF_PL_epsi_1'
                ep1 = r.traj.(r.c_obs.model).epsi(:,1);
                ep1(r.irr) = [];
                logresp = be06 +be6.*ep1;
            case 'HGF_PL_epsi_2'
                ep2 = r.traj.(r.c_obs.model).epsi(:,2);
                ep2(r.irr) = [];
                logresp = be08 +be8.*ep2;
            case 'HGF_PL_epsi_3'
                ep3 = r.traj.(r.c_obs.model).epsi(:,3);
                ep3(r.irr) = [];
                logresp = be010 +be10.*ep3;
            case 'HGF_PL_epsi-abs_1'
                ep1 = abs(r.traj.(r.c_obs.model).epsi(:,1));
                ep1(r.irr) = [];
                logresp = be06 +be6.*ep1;
            case 'HGF_PL_epsi-abs_2'
                ep2 = abs(r.traj.(r.c_obs.model).epsi(:,2));
                ep2(r.irr) = [];
                logresp = be08 +be8.*ep2;
            case 'HGF_PL_epsi-abs_3'
                ep3 = abs(r.traj.(r.c_obs.model).epsi(:,3));
                ep3(r.irr) = [];
                logresp = be010 +be10.*ep3;
            case 'HGF_PL_mu_1'
                mu1 = r.traj.(r.c_obs.model).mu(:,1);
                mu1(r.irr) = [];
                logresp = be012 +be12.*mu1;
            case 'HGF_PL_mu_2'
                mu2 = r.traj.(r.c_obs.model).mu(:,2);
                mu2(r.irr) = [];
                logresp = be013 +be13.*mu2;
            case 'HGF_PL_mu_3'
                mu3 = r.traj.(r.c_obs.model).mu(:,3);
                mu3(r.irr) = [];
                logresp = be014 +be14.*mu3;
            case 'HGF_PL_sa_1'
                sa1 = r.traj.(r.c_obs.model).sa(:,1);
                sa1(r.irr) = [];
                logresp = be02 +be2.*sa1;
            case 'HGF_PL_sa_2'
                sa2 = r.traj.(r.c_obs.model).sa(:,2);
                sa2(r.irr) = [];
                logresp = be02 +be2.*sa2;
            case 'HGF_PL_sa_13'
                sa3 = r.traj.(r.c_obs.model).sa(:,3);
                sa3(r.irr) = [];
                logresp = be02 +be2.*sa3;
            case 'HGF_PL_inferv_1'
                mu2 = r.traj.(r.c_obs.model).mu(:,2);
                mu2(r.irr) = [];
                sa2 = r.traj.(r.c_obs.model).sa(:,2);
                sa2(r.irr) = [];
                sigmoid_mu2 = 1./(1+exp(-mu2)); % transform down to 1st level
                inferv = sigmoid_mu2.*(1 -sigmoid_mu2).*sa2; 
                logresp = be03 +be3.*inferv;
            case 'HGF_PL_pv_1'
                mu2 = r.traj.(r.c_obs.model).mu(:,2);
                mu2(r.irr) = [];
                mu3 = r.traj.(r.c_obs.model).mu(:,3);
                mu3(r.irr) = [];
                sigmoid_mu2 = 1./(1+exp(-mu2)); % transform down to 1st level
                pv = sigmoid_mu2.*(1-sigmoid_mu2).*exp(mu3); 
                logresp = be04 +be4.*pv;
        end
        logp_reg(reg) = logp_reg(reg) + -1/2.*log(8*atan(1).*ze) -(yr(:,y)-logresp).^2./(2.*ze);
        
    end


end


%% Regression with CCA components
if any(strcmp(r.c_obs.responses, 'CCAcomp'))

    % Create param struct
    type='reg'; % regression
    for pn=1:length(nme)
        if strcmp(nme{pn,1}(1:length(type)),type)
            eval([nme_gen{pn} ' = pvec(idx{pn})'';']);
        end
    end
    
    % switch the data input column
    ycol = 3:size(r.y,2);
    if isfield(r.c_obs,'ynum') && ~isempty(r.c_obs.ynum)
        ycol = ycol(r.c_obs.ynum);
    end

    % Weed irregular trials out from responses and inputs
    yr = r.y(:,ycol); % RTs are in column 2, EEG in column 3 or greater
    yr(r.irr,:) = [];
    u = r.u(:,1);
    u(r.irr) = [];
    
    
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
    
    
    % prediction error
    % ~~~~~~~~
    if abs_switch
        daureg = abs(dau);
        daureg(r.irr) = [];
        ep1reg = abs(ep1);
        ep1reg(r.irr) = [];
        da1reg = abs(da1);
        da1reg(r.irr) = [];
        ep2reg = abs(ep2);
        ep2reg(r.irr) = [];
        if l>2
            da2reg = abs(da2);
            da2reg(r.irr) = [];
            ep3reg = abs(ep3);
            ep3reg(r.irr) = [];
            da3reg = abs(da3);
            da3reg(r.irr) = [];
        end
    else
        daureg = dau;
        daureg(r.irr) = [];
        ep1reg = ep1;
        ep1reg(r.irr) = [];
        da1reg = da1;
        da1reg(r.irr) = [];
        ep2reg = ep2;
        ep2reg(r.irr) = [];
        if l>2
            da2reg = da2;
            da2reg(r.irr) = [];
            ep3reg = ep3;
            ep3reg(r.irr) = [];
            da3reg = da3;
            da3reg(r.irr) = [];
        end
    end
    
    % Posterior expectation
    % ~~~~~~~~
    m1reg = mu1;
    m1reg(r.irr) = [];
    sa1reg = sa1;
    sa1reg(r.irr) = [];
    if l>1
        m2reg = mu2;
        m2reg(r.irr) = [];
        sa2reg = sa2;
        sa2reg(r.irr) = [];
    end
    if l>2
        m3reg = mu3;
        m3reg(r.irr) = [];
        sa3reg = sa3;
        sa3reg(r.irr) = [];
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
    bernv(r.irr) = [];
    
    bernvhat = sa1hat;
    bernvhat(r.irr) = [];
    
    mu1hreg = mu1hat;
    mu1hreg(r.irr) = [];

    if l>1 % CAB
        % Inferential variance (aka informational or estimation uncertainty, ambiguity)
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %inferv = tapas_sgm(mu2, 1).*(1 -tapas_sgm(mu2, 1)).*sa2; % transform down to 1st level
        sigmoid_mu2 = 1./(1+exp(-mu2)); % transform down to 1st level
        inferv = sigmoid_mu2.*(1 -sigmoid_mu2).*sa2; 
        %inferv = sa2; 
        inferv(r.irr) = [];
        
        mu2hreg = mu2hat;
        mu2hreg(r.irr) = [];
        sa2hreg = sa2hat;
        sa2hreg(r.irr) = [];
    end

    if l>2 % CAB
        % Phasic volatility (aka environmental or unexpected uncertainty)
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pv = sigmoid_mu2.*(1-sigmoid_mu2).*exp(mu3); % transform down to 1st level
        %pv = exp(mu3); 
        pv(r.irr) = [];
        
    end
    
    logp_reg(reg,1) = 0;
    for y=1:size(yr,2)
        param_nums = r.c_obs.params_cell{y};

        % assign intercept
        eval(['be0 = be0' num2str(y) ';']);
        
        % assign other params to zero if not needed
        for pn = 1:16
            if ~any(param_nums==pn)
                eval(['be' num2str(pn) ' = 0;']);
            end
        end
        
        logresp = be0 +be1.*surp +be2.*bernv +be3.*inferv +be4.*pv +be5.*daureg +be6.*ep1reg +be7.*da1reg +be8.*ep2reg +be9.*da2reg +be10.*ep3reg +be11.*da3reg +be12.*m1reg +be13.*m2reg +be14.*m3reg +be15.*sa1reg +be16.*sa2reg +be17.*sa3reg;

        logp_reg(reg) = logp_reg(reg) + -1/2.*log(8*atan(1).*ze) -(yr(:,y)-logresp).^2./(2.*ze);
        
    end


end

if ~exist('logp_so','var')
    logp_so=0;
end
if ~exist('logp_reg','var')
    logp_reg=0;
end

% Joint probability:
% The logarithm of the probability of multiple joint probabilities simplifies 
% to the sum of the logarithms of the individual probabilities
logp = logp_so + logp_reg;

return;
