function [traj,infStates] = GBM(r, pvec, varargin)
% GENERAL BINARY (INPUT) MODEL

% NOTES
% 1. in and N defined separately
% 2. representation Mu has no variance in SDT, there is just noise (alpha)
% that influences it's mean value.
% 3. must be a two-level model because the 2nd level sets the prior for the
% first, even if the prior is static

% PLANNED UPDATES
% 1. Allow separate evaluation of alpha for targets and non-targets - done
% 2. Estimate prior expectation (static)
% 3. Prior and alpha variable by block type
% 4. Prior and alpha variable by trial (Kalman)
% 5. Hierarchical priors

if ~isfield(r,'eGBM')
    r.eGBM=0;
end

% Transform paramaters back to their native space if needed
if ~isempty(varargin) 
    if strcmp(varargin{1},'trans')
        pvec = GBM_transp(r, pvec);
    end
end

% Add dummy "zeroth" trial
u = [zeros(1,size(r.u,2)); r.u(:,:)];

% Number of trials (including prior)
n = size(u,1);

% Assume that if u has more than one column, the last contains t
try
    if r.c_prc.irregular_intervals
        if size(u,2) > 1
            t = [0; r.u(:,end)];
        else
            error('tapas:hgf:InputSingleColumn', 'Input matrix must contain more than one column if irregular_intervals is set to true.');
        end
    else
        t = ones(n,1);
    end
catch
    if size(u,2) > 1
        t = [0; r.u(:,end)];
    else
        t = ones(n,1);
    end
end

%% INITIALISE

% Create param struct
nme=r.c_prc.pnames;
nme_gen=r.c_prc.pnames_gen;
idx=r.c_prc.priormusi;
type='like';
for pn=1:length(nme)
    if strcmp(nme{pn,1}(1:length(type)),type)
        eval([nme_gen{pn} ' = pvec(idx{pn})'';']);
    end
end

%levels
for m=1:r.c_prc.nModels
    type = r.c_prc.type{m};
    l(m) = r.c_prc.(type).n_priorlevels+1;
end
maxlev=max(l);

% Representations
nmod=r.c_prc.nModels;
mu0 = NaN(n,maxlev,nmod);
mu = NaN(n,maxlev,nmod);
pi = NaN(n,maxlev,nmod);

% Other quantities
muhat = NaN(n,maxlev,nmod);
pihat = NaN(n,maxlev,nmod);
v     = NaN(n,maxlev,nmod);
w     = NaN(n,maxlev,nmod);
da    = NaN(n,maxlev,nmod);
dau   = NaN(n,1,nmod);
g  = NaN(n,maxlev,nmod); % Kalman gain (optional)

al = NaN(n,1);
% joint trajectories (if needed)
vj_mu = NaN(n,1);
xc = NaN(n,1);
xchat = NaN(n,1);
xcpihat = NaN(n,1);

% parameters
th = NaN(1,1,nmod);
    
for m=1:r.c_prc.nModels
    type = r.c_prc.type{m};
    
    % Create param struct
    for pn=1:length(nme)
        if strcmp(nme{pn,1}(1:length(type)),type)
            lenpvec=length(pvec(idx{pn}));
            thispvec=nan(lenpvec,1);
            thispvec(1:lenpvec)=pvec(idx{pn})';
            eval([nme_gen{pn} '(1:lenpvec,1,m) = thispvec;']);
        end
    end

    if exist('om','var')
        if l(m)==3
            th(1,1,m)   = exp(om(3,1,m));
            om(3,1,m)=NaN;
        elseif l(m)==2
            th(1,1,m)   = exp(om(2,1,m));
            om(2,1,m)=NaN;
        end
    end

    % Representation priors
    % Note: first entries of the other quantities remain
    % NaN because they are undefined and are thrown away
    % at the end; their presence simply leads to consistent
    % trial indices.
    mu(1,1,m) = 1./(1+exp(-mu_0(1,1,m)));
    pi(1,1,m) = Inf;
    if l(m)>1
        mu(1,2:end,m) = mu_0(2:l(m),1,m);
        pi(1,2:end,m) = 1./sa_0(2:l(m),1,m);
    end
    if exist('g_0','var') && size(g_0,3)>=m
        g(1,:,m)  = g_0(1,1,m); % Kalman gain (optional)
        expom(:,1,m) = exp(om(:,1,m));
    end
end
xc(1) = 0.5; % for joint models

%% Efficiacy
% Unpack from struct to double for efficiency
% type='like';
% pnames = fieldnames(M.(type).p);
% for pn=1:length(pnames)
%     eval([pnames{pn} ' = M.(type).p.(pnames{pn});']);
% end
% for m=1:r.c_prc.nModels
%     type = r.c_prc.type{m};
% 
%     % Unpack prior parameters
%     pnames = fieldnames(M.(type).p);
%     for pn=1:length(pnames)
%         eval([pnames{pn} '(:,:,m) = M.(type).p.(pnames{pn})'';']);
%     end
% 
%     %Unpack prior traj
% %     tnames = fieldnames(M.(type).tr);
% %     for tn=1:length(tnames)
% %         eval([tnames{tn} '(:,1:size(M.(type).tr.(tnames{tn}),2),m) = M.(type).tr.(tnames{tn});']);
% %     end
% end

% find out if certain variables exist or not
% no_al1=0;
% if ~exist('al1','var')
%     no_al1=1;
% end
no_rb=0;
if ~exist('rb','var')
    no_rb=1;
end
no_g=0;
if ~exist('g','var')
    no_g=1;
end

% get model structure
for m=1:r.c_prc.nModels
    type = r.c_prc.type{m};
    hierarchical(m)=0;
    static(m)=0;
    dynamic(m)=0;
    state(m)=0;
    if strcmp(r.c_prc.(type).priortype,'hierarchical')
        hierarchical(m)=1;
        if strcmp(r.c_prc.(type).priorupdate,'static')
            static(m)=1;
        elseif strcmp(r.c_prc.(type).priorupdate,'dynamic')
            dynamic(m)=1;
        end
    elseif strcmp(r.c_prc.(type).priortype,'state')
        state(m)=1;
    end
    AL(m)=0;
    PL(m)=0;
    if strcmp(r.c_prc.type{m},'AL')
        AL(m)=1;
    elseif strcmp(r.c_prc.type{m},'PL')
        PL(m)=1;
    end
end
inputfactors=r.c_prc.inputfactors;
n_inputfactors=length(inputfactors);
nModels=r.c_prc.nModels;

% ignore trials
no_ign = ones(1,n);
no_ign(r.ign) = 0;

%% UPDATES
for k=2:1:n
    
    for m=1:nModels
        if no_ign(k-1)
            
        
            %% Predictions (from previous trial or static parameters)
            if hierarchical(m)

                % 2nd level prediction
                muhat(k,2,m) = mu(k-1,2,m) +t(k) *rho(2,1,m);

                % Prediction from level 2 (which can be either static or dynamic)
                muhat(k,1,m) = 1./(1+exp(-muhat(k,2,m)));

            elseif state(m)
                % Prediction from prior state, e.g. Kalman filter
                muhat(k,1,m) =  mu(k-1,1,m);

            end

            % Precision of prediction
            pihat(k,1,m) = 1/(muhat(k,1,m)*(1 -muhat(k,1,m)));


            %% Updates
            if r.eGBM
                pi(k,1,m) = Inf; % eHGF
            end
            % Value prediction error, e.g. for Kalman filter, also known as the
            % "innovation"
            dau(k,1,m) = u(k,1) -muhat(k,1,m);

            al(k,1,m)=al0(1);
            if n_inputfactors >0
                % apply gain to alpha according to input factors
                % alpha multipliers are not variances, but have to be positive
                if_str = num2str(u(k,2));
                for nif = 1:n_inputfactors 
                    if strcmp(if_str(inputfactors(nif)),'1') % e.g. lower intensity stimulus where options are '1' or '2' - more uncertainty expected for lower intensity
                        al(k,1,m)=al(k,1,m) * al0(inputfactors(nif)+1);
                    end
                end
            end
            
            if hierarchical(m)
                % Likelihood functions: one for each
                % possible signal
%                 if n_inputcond >1
%                     und1 = exp(-(u(k) -eta1)^2/(2*al1(u(k,2))));
%                     und0 = exp(-(u(k) -eta0)^2/(2*al0(u(k,2))));
%                 else
%                     und1 = exp(-(u(k) -eta1)^2/(2*al1));
%                     und0 = exp(-(u(k) -eta0)^2/(2*al0));
%                 end
                if al0(1)==0 % certain
                    mu(k,1,m)=u(k,1);
                    mu0(k,1,m)=u(k,1);
                else % uncertain
                    und1 = exp(-(u(k,1) -eta1)^2/(2*al(k,1,m)));
                    und0 = exp(-(u(k,1) -eta0)^2/(2*al(k,1,m)));

                    if AL(m)
                        if u(k,3)==2
                            mu0(k,1,m) = muhat(k,1,m) *und1 /(muhat(k,1,m) *und1 +(1 -muhat(k,1,m)) *und0);
                            mu(k,1,m) = mu0(k,1,m);
                        elseif u(k,3)==1
                            mu0(k,1,m) = (1-muhat(k,1,m)) *und1 /(muhat(k,1,m) *und0 +(1 -muhat(k,1,m)) *und1);
                            mu(k,1,m) = 1-mu0(k,1,m);
                        end
                    else
                        mu(k,1,m) = muhat(k,1,m) *und1 /(muhat(k,1,m) *und1 +(1 -muhat(k,1,m)) *und0);
                        mu0(k,1,m) = mu(k,1,m);
                    end
                end


                %%
                % Representation prediction error
                da(k,1,m) = mu(k,1,m) -muhat(k,1,m);

                % second level predictions and precisions
                % Precision of prediction
                if l(m)>2 
                    % if there is a third level
                    pihat(k,2,m) = 1/(1/pi(k-1,2,m) +exp(ka(2,1,m) *mu(k-1,3,m) +om(2,1,m)));
                else
                    % otherwise, use theta since this is the highest
                    % level
                    pihat(k,2,m) = 1/(1/pi(k-1,2,m)  +t(k) *th(1,1,m));
                end
                % Updates
                pi(k,2,m) = pihat(k,2,m) +1/pihat(k,1,m);
                mu(k,2,m) = muhat(k,2,m) +1/pi(k,2,m) *da(k,1,m);

                % Volatility prediction error
                da(k,2,m) = (1/pi(k,2,m) +(mu(k,2,m) -muhat(k,2,m))^2) *pihat(k,2,m) -1;
                

                if l(m) > 3
                    % Pass through higher levels
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~
                    for j = 3:l(m)-1
                        % Prediction
                        muhat(k,j,m) = mu(k-1,j,m) +t(k) *rho(j,1,m);

                        % Precision of prediction
                        pihat(k,j,m) = 1/(1/pi(k-1,j,m) +t(k) *exp(ka(j,1,m) *mu(k-1,j+1,m) +om(j,1,m)));

                        % Weighting factor
                        v(k,j-1,m) = t(k) *exp(ka(j-1,1,m) *mu(k-1,j,m) +om(j-1,1,m));
                        w(k,j-1,m) = v(k,j-1,m) *pihat(k,j-1,m);

                        
                        if r.eGBM
                            % Mean update
                            mu(k,j,m) = muhat(k,j,m) +1/2 *1/pihat(k,j,m) *ka(j-1,1,m) *w(k,j-1,m) *da(k,j-1,m);
            
                            % Ingredients of precision update which depend on the mean
                            % update
                            vv = t(k) *exp(ka(j-1,1,m) *mu(k,j,1) +om(j-1),1,m);
                            pimhat = 1/(1/pi(k-1,j-1,m) +vv); 
                            ww = vv *pimhat;
                            rr = (vv -1/pi(k-1,j-1,m)) *pimhat;
                            dd = (1/pi(k,j-1,m) +(mu(k,j-1,m) -muhat(k,j-1,m))^2) *pimhat -1;
                            
                            % Precision update
                            pi(k,j,m) = pihat(k,j,m) +max(0, 1/2 *ka(j-1,1,m)^2 *ww*(ww +rr*dd));
                        else
                            % Updates
                            pi(k,j,m) = pihat(k,j,m) +1/2 *ka(j-1,1,m)^2 *w(k,j-1,m) *(w(k,j-1,m) +(2 *w(k,j-1,m) -1) *da(k,j-1,m));
    
                            if pi(k,j,m) <= 0
                                error('tapas:hgf:NegPostPrec', 'Negative posterior precision. Parameters are in a region where model assumptions are violated.');
                            end
    
                            mu(k,j,m) = muhat(k,j,m) +1/2 *1/pi(k,j,m) *ka(j-1,1,m) *w(k,j-1,m) *da(k,j-1,m);
                        
                        end

                        % Volatility prediction error
                        da(k,j,m) = (1/pi(k,j,m) +(mu(k,j,m) -muhat(k,j,m))^2) *pihat(k,j,m) -1;
                    end
                end
                
                if l(m)>2
                    % Last level
                    % ~~~~~~~~~~
                    % Prediction
                    muhat(k,l(m),m) = mu(k-1,l(m),m) +t(k) *rho(l(m),1,m);

                    % Precision of prediction
                    pihat(k,l(m),m) = 1/(1/pi(k-1,l(m),m) +t(k) *th(1,1,m));

                    % Weighting factor
                    v(k,l(m),m)   = t(k) *th(1,1,m);
                    v(k,l(m)-1,m) = t(k) *exp(ka(l(m)-1,1,m) *mu(k-1,l(m),m) +om(l(m)-1,1,m));
                    w(k,l(m)-1,m) = v(k,l(m)-1,m) *pihat(k,l(m)-1,m);
                    
                    if r.eGBM
                        % Mean update
                        mu(k,l(m),m) = muhat(k,l(m),m) +1/2 *1/pihat(k,l(m),m) *ka(l(m)-1,1,m) *w(k,l(m)-1,m) *da(k,l(m)-1,m);
                        
                        % Ingredients of the precision update which depend on the mean
                        % update
                        vv = t(k) *exp(ka(l(m)-1,1,m) *mu(k,l(m),m) +om(l(m)-1,1,m));
                        pimhat = 1/(1/pi(k-1,l(m)-1,m) +vv); 
                        ww = vv *pimhat;
                        rr = (vv -1/pi(k-1,l(m)-1,m)) *pimhat;
                        dd = (1/pi(k,l(m)-1,m) +(mu(k,l(m)-1,m) -muhat(k,l(m)-1,m))^2) *pimhat -1;
                                
                        pi(k,l(m),m) = pihat(k,l(m),m) +max(0, 1/2 *ka(l(m)-1,1,m)^2 *ww*(ww +rr*dd));
                    else
                        % Updates
                        pi(k,l(m),m) = pihat(k,l(m),m) +1/2 *ka(l(m)-1,1,m)^2 *w(k,l(m)-1,m) *(w(k,l(m)-1,m) +(2 *w(k,l(m)-1,m) -1) *da(k,l(m)-1,m));
    
                        if pi(k,l(m),m) <= 0
                            error('tapas:hgf:NegPostPrec', 'Negative posterior precision. Parameters are in a region where model assumptions are violated.');
                        end
    
                        mu(k,l(m),m) = muhat(k,l(m),m) +1/2 *1/pi(k,l(m),m) *ka(l(m)-1,1,m) *w(k,l(m)-1,m) *da(k,l(m)-1,m);
                    end

                    % Volatility prediction error
                    da(k,l(m),m) = (1/pi(k,l(m),m) +(mu(k,l(m),m) -muhat(k,l(m),m))^2) *pihat(k,l(m),m) -1;
                end
                
                % Implied posterior precision at first level
                sgmmu2 = 1./(1+exp(-mu(k,2,m)));
                pi(k,1,m) = pi(k,2,m)/(sgmmu2*(1-sgmmu2));


            elseif state(m) 
                if ~isinf(expom(1,1,m)) % otherwise, RP
                    %nModels==1 % Kalman - CANNOT APPLY TO HYBRID MODELS? Yes, can
                    % Gain update - optimal gain is calculated from ratio of input
                    % variance to representation variance
    
                    % Same gain function modified by two different alphas
                    if expom(1,1,m)==0
                        % no 1st level uncertainty - similar to RW learning
                        g(k,1,m) = al(k,1,m);

                        % Ensure numerical stability by avoiding extremes
                        muhat(k,1,m) = max(muhat(k,1,m), 0.001);
                        muhat(k,1,m) = min(muhat(k,1,m), 0.999);
        
                    else
                        % Kalman gain
                        g(k,1,m) = (g(k-1,1,m) +al(k,1,m)*expom(1,1,m))/(g(k-1,1,m) +al(k,1,m)*expom(1,1,m) +1);
                    end
                    % Hidden state mean update
                    mu(k,1,m) = muhat(k,1,m)+g(k,1,m)*dau(k,1,m);
                    pi(k,1,m) = (1-g(k,1,m)) *al(k,1,m)*expom(1,1,m); 

                    if AL(m)
                        % convert cue-outcome association back to outcome
                        if u(k,3)==2
                            mu0(k,1,m) = mu(k,1,m);
                        elseif u(k,3)==1
                            mu0(k,1,m) = 1-mu(k,1,m);
                        end
                    else
                        mu0(k,1,m) = mu(k,1,m);
                    end

                    % Alternative: separate gain functions for each stimulus type
              %      if r.c_prc.(type).one_alpha
              %          pi_u=al0(u(k,2));
              %          g(k,1) = (g(k-1,1) +pi_u*expom)/(g(k-1,1) +pi_u*expom +1);
              %          % Hidden state mean update
              %          mu(k,1,m) = muhat(k,1,m)+g(k,1)*dau(k);
              %          pi(k,1,m) = (1-g(k,1)) *pi_u*expom;
              %      else
              %          if u(k,1)==0
              %              pi_u=al0(u(k,2));
              %              g(k,1) = (g(k-1,1) +pi_u*expom)/(g(k-1,1) +pi_u*expom +1);
              %              g(k,2) = g(k-1,2);
              %              % Hidden state mean update
              %              mu(k,1,m) = muhat(k,1,m)+g(k,1)*dau(k);
              %              pi(k,1,m) = (1-g(k,1)) *pi_u*expom;
              %          elseif u(k,1)==1
              %              pi_u=al1(u(k,2));
              %              g(k,2) = (g(k-1,2) +pi_u*expom)/(g(k-1,2) +pi_u*expom +1);
              %              g(k,1) = g(k-1,1);
              %              % Hidden state mean update
              %              mu(k,1,m) = muhat(k,1,m)+g(k,2)*dau(k);
              %              pi(k,1,m) = (1-g(k,2)) *pi_u*expom;
              %          end
              %      end
    
                    % Representation prediction error
                    da(k,1,m) = mu(k,1,m) -muhat(k,1,m);
                else
                    % just make these zero for plotting, but they are not
                    % used for estimating the model (only muhat for joint models)
                    muhat(k,1,m) =  xc(k-1,1);
                    mu(k,1,m) = 0;
                    mu0(k,1,m) = 0;
                end

            end

            % RESPONSE BIAS
            if ~no_rb
                mu(k,1,m) = mu(k,1,m)+rb(1,1,m);
            end

        else
            mu(k,:,m) = mu(k-1,:,m); 
            pi(k,:,m) = pi(k-1,:,m);

            muhat(k,:,m) = muhat(k-1,:,m);
            pihat(k,:,m) = pihat(k-1,:,m);

            v(k,:,m)  = v(k-1,:,m);
            w(k,:,m)  = w(k-1,:,m);
            da(k,:,m) = da(k-1,:,m);
            dau(k,1,m) = dau(k-1,1,m);
            if no_g==0
                g(k,:,m)=g(k-1,:,m);
            end
            al(k,1,m)  = al(k-1,1,m);
        end
    end
    
    % Joint prediction if more than one model
    if nModels>1
        
        % set alpha - allow m to be the max model index, m (since it doesn't
        % matter which we use for estimation)
        al(k,1,m)=al0(1);
        if n_inputfactors >0
            % apply gain to alpha according to input factors
            % alpha multipliers are not variances, but have to be positive
            if_str = num2str(u(k,2));
            for nif = 1:n_inputfactors 
                if strcmp(if_str(inputfactors(nif)),'1')
                    al(k,1,m)=al(k,1,m) * al0(inputfactors(nif)+1);
                end
            end
        end

        if exist('xi','var')
            % convert xi’s into simplex weights
            expvec = exp(xi(1,1,:));
            phi    = expvec ./ sum(expvec);    
            
            % % 4) final branch precisions --------------------------------------------
            % for m = 1:nModels
            %     phi(1,1,m) = phiT * wphi(m);     % π_m  for branch m
            % end
        end
                
        % joint probability
        vj_phi = sum(phi(1,1,:));
%         vj_mu(k,1) = sum(phi(1,1,:).*mu0(k,1,:))/vj_phi;
        vj_mu(k,1) = sum(phi(1,1,:).*muhat(k,1,:))/vj_phi; % new - weight muhat for use by RT response model - surprise
        xchat(k,1) = vj_mu(k,1);

        % joint prediction
            % ((eta1-vj_mu(k,1))^2 - (eta0-vj_mu(k,1))^2) gives a value between
            % -1 (if vj_mu is closer to eta1) and 1 (if closer to eta0).
            % dividing by vj_phi makes log(rt) either very large positive (if phi
            % is very precise) or very large negative (if uncertain). Exp of
            % these results in a large positive or small positive number
            % respectively.
            % So, when vj_mu is close to 0 and precise, xchat approaches 0
        if 0
            % bimodal likelihood - calculate xc here, not xchat
            rt=exp(((eta1-vj_mu(k,1))^2 - (eta0-vj_mu(k,1))^2)/(vj_phi^-2));
            xc(k,1) = 1/(1+rt);
%             rt=exp(((eta1-vj_muhat(k,1))^2 - (eta0-vj_muhat(k,1))^2)/(vj_phi^-2));
%             xchat(k,1) = 1/(1+rt);
%         elseif 0
%             % bimodal likelihood - probably wrong
%             rt=exp(((eta1-vj_mu(k,1))^2 - (eta0-vj_mu(k,1))^2)/(vj_phi^-2));
%             xchat(k,1) = 1/(1+rt);  
        end
        
        % Precision of prediction - INCORRECT SINCE XCHAT IS EQUIVALENT TO
        % MU NOT MUHAT
        %xcpihat(k,1) = 1/(xchat(k,1)*(1 -xchat(k,1)));
        
        if 1
            % Likelihood functions: one for each
            % possible signal - INCORRECT SINCE WE ALREADY APPLIED THIS FOR
            % AL/PL MODELS
            und1 = exp(-(u(k) -eta1)^2/(2*al(k,1,m)));
            und0 = exp(-(u(k) -eta0)^2/(2*al(k,1,m)));
    
            % Update
            xc(k,1) = xchat(k,1) *und1 /(xchat(k,1) *und1 +(1 -xchat(k,1)) *und0);
        end
        
    end
    
end

%% COMPILE RESULTS

for m=1:nModels
    if l(m)>1
%         if r.eGBM
%             sgmmu2 = tapas_sgm(mu(:,2,m), 1);
%             dasgmmu2 = u(:,1) -sgmmu2;
%             lr1    = diff(sgmmu2)./dasgmmu2(2:n,1);
%             lr1(da(2:n,1,m)==0) = 0;
%         else
            % Implied learning rate at the first level
            sgmmu2 = 1./(1+exp(-mu(:,2,m)));
            lr1    = diff(sgmmu2)./da(2:n,1,m);
            lr1(da(2:n,1,m)==0) = 0;
%         end
    end
end

% Remove representation priors
mu0(1,:,:)  = [];
mu(1,:,:)  = [];
pi(1,:,:)  = [];
al(1,:,:)     = [];
 
% joint trajectories (if needed)
vj_mu(1) = [];
xc(1) = [];
xchat(1) = [];
xcpihat(1) = [];

% Remove other dummy initial values
muhat(1,:,:) = [];
pihat(1,:,:) = [];
v(1,:,:)     = [];
w(1,:,:)     = [];
da(1,:,:)    = [];
dau(1,:,:)     = [];
if no_g==0
    g(1,:,:)  = [];
end

% Create result data structure
traj = struct;

for m=1:nModels

%     % Unpack likelihood parameters
%     type='like';
%     pnames = fieldnames(M.(type).p);
%     for pn=1:length(pnames)
%         p.(pnames{pn}) = (pnames{pn});
%     end

    type = r.c_prc.type{m};
    
%     % Unpack prior parameters
%     pnames = fieldnames(M.(type).p);
%     for pn=1:length(pnames)
%         p.(pnames{pn}) = (pnames{pn});
%     end
% 
%     %Unpack prior traj
%     tnames = fieldnames(M.(type).tr);
%     for tn=1:length(tnames)
%         tr.(tnames{tn}) = (tnames{tn});
%     end


    % Check validity of trajectories
    if any(isnan(mu(:,:,m))) %|| any(isnan(pi(:)))
        error('tapas:hgf:VarApproxInvalid', 'Variational approximation invalid. Parameters are in a region where model assumptions are violated.');
    else
        % Check for implausible jumps in trajectories
        % CAB: only use first 500 trials - after that changes in precision become too small
        ntrials = min(length(mu(:,:,m)),500);
        dmu = diff(mu(1:ntrials,2:end,m));
        dpi = diff(pi(1:ntrials,2:end,m));
        rmdmu = repmat(sqrt(mean(dmu.^2)),length(dmu),1);
        rmdpi = repmat(sqrt(mean(dpi.^2)),length(dpi),1);

        jumpTol = 16;
        if any(abs(dmu(:)) > jumpTol*rmdmu(:)) || any(abs(dpi(:)) > jumpTol*rmdpi(:))
            %disp('hgf:VarApproxInvalid: GBM Variational approximation invalid. Parameters are in a region where model assumptions are violated.');
            %disp('Use plot for diagnosis: see within function'); 
            % figure; plot(abs(dpi(:))); hold on; plot(rmdpi(:),'r'); hold on; plot(jumpTol*rmdpi(:),'g')
            % plot(abs(dmu(:))); hold on; plot(rmdmu(:),'r'); hold on; plot(jumpTol*rmdmu(:),'g')
            infStates = [];
            return
        end
    end

    traj.like.vj_mu = vj_mu;
    traj.like.xc = xc;
    traj.like.xchat = xchat;
    traj.like.xcsahat = 1./xcpihat;

    traj.(type).mu0    = mu0(:,:,m);
    traj.(type).mu     = mu(:,:,m);
    traj.(type).sa     = 1./pi(:,:,m);

    traj.(type).muhat  = muhat(:,:,m);
    traj.(type).sahat  = 1./pihat(:,:,m);
    traj.(type).v      = v(:,:,m);
    traj.(type).w      = w(:,:,m);
    traj.(type).da     = da(:,:,m);
    traj.(type).dau    = dau(:,:,m);
    if no_g==0 && size(g,3)>=l(m)
        traj.(type).g     = g(:,:,m);
    end

    % Updates with respect to prediction
    traj.(type).ud = muhat(:,:,m) -mu(:,:,m);

    % Psi (precision weights on prediction errors)
    psi        = NaN(n-1,l(m)+1);
    psi(:,1)   = 1./(al(:,:,m).*pi(:,1,m)); % dot multiply only if al is a vector. More simply: psi1 = precision of input ./ 1st level precision
    if state(m); psi(:,2)   = 1./pi(:,1,m);end
    if l(m)>1; psi(:,2)   = 1./pi(:,2,m);end
    if l(m)>2; psi(:,3:l(m)) = pihat(:,2:l(m)-1,m)./pi(:,3:l(m),m);end
    traj.(type).psi   = psi;

    % Epsilons (precision-weighted prediction errors)
    epsi        = NaN(n-1,l(m));
    epsi(:,1)   = psi(:,1) .*dau(:,:,m);
    if state(m); epsi(:,2) = psi(:,2) .*da(:,1,m);end
    if l(m)>1; epsi(:,2:l(m)) = psi(:,2:l(m)) .*da(:,1:l(m)-1,m);end
    traj.(type).epsi   = epsi;

    % Full learning rate (full weights on prediction errors)
    wt        = NaN(n-1,l(m));
    if l(m)==1; lr1=psi(:,1);end
    wt(:,1)   = lr1;
    if l(m)>1; wt(:,2)   = psi(:,2); end
    if l(m)>2; wt(:,3:l(m)) = 1/2 *(v(:,2:l(m)-1,m) *diag(ka(2:l(m)-1,m))) .*psi(:,3:l(m)); end
    traj.(type).wt   = wt;

    % Create matrices for use by tapas observation models
     infStates = NaN(n-1,l(m),11);
%      infStates(:,:,1) = traj.(type).muhat;
%      infStates(:,:,2) = traj.(type).sahat;
%      infStates(:,:,3) = traj.(type).mu;
%      infStates(:,:,4) = traj.(type).sa;
%      infStates(:,:,5) = traj.(type).da;
%      infStates(:,:,6) = traj.(type).epsi;
%      infStates(:,1,7) = traj.(type).dau;
%     infStates(:,1,8) = traj.mu0;
%     infStates(:,1,9) = traj.vj_mu;
%     infStates(:,1,10) = traj.xc;
%     infStates(:,1,11) = traj.xchat;
end

return;
