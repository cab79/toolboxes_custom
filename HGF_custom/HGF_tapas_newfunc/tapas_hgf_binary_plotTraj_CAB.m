function tapas_hgf_binary_plotTraj_CAB(r,priormodel,plotsd)
% Plots the estimated or generated trajectories for the binary HGF perceptual model
% Usage example:  est = tapas_fitModel(responses, inputs); tapas_hgf_binary_plotTraj(est);
%
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2012-2013 Christoph Mathys, TNU, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.

% Optional plotting of standard deviations (true or false)
%plotsd = true;

% Optional plotting of responses (true or false)
ploty = true;
fontsize=16;

% plot colours
% input_colour = {'green', [0 0.6 0]};
input_colour = {'red', [1 0 0]};
% mu1_colour = {'red', [0.6 0 0]};
mu1_colour = {'purple', [0.6 0 0.6]};
mu2_colour = {'blue', [0 0 1]};
mu3_colour = {'green', [0 0.8 0]};

% Set up display
scrsz = get(0,'screenSize');
outerpos = [0.2*scrsz(3),0.2*scrsz(4),0.8*scrsz(3),0.8*scrsz(4)];
figure(...
    'OuterPosition', outerpos,...
    'Name', ['HGF trajectories ' priormodel]);

% Time axis
%if size(r.u,2) > 1
%    t = r.u(:,end)';
%else
    t = ones(1,size(r.u,1));
%end

ts = cumsum(t);
ts = [0, ts];

% Number of levels
%try
    l = r.c_prc.(priormodel).n_priorlevels+1;
%catch
%    l = (length(r.p_prc.p)+1)/5;
%end

mu_0=r.p_prc.([priormodel '_mu_0']);
mu = r.traj.(priormodel).mu;
muhat = r.traj.(priormodel).muhat;
wt = r.traj.(priormodel).wt;
if isfield(r.p_prc,[priormodel '_om'])
    om = r.p_prc.([priormodel '_om']);
else
    om = NaN(1,l);
end
if isfield(r.p_prc,[priormodel '_rho'])
    rho = r.p_prc.([priormodel '_rho']);
else
    rho = NaN(1,l);
end
if isfield(r.p_prc,[priormodel '_ka'])
    ka = r.p_prc.([priormodel '_ka']);
else
    ka = NaN(1,l);
end

if plotsd
    if isfield(r.traj.(priormodel),'g') && l==1
        sa = r.traj.(priormodel).sa;
        sa_0=sa(1);
    else
        sa_0=r.p_prc.([priormodel '_sa_0']);
        sa = r.traj.(priormodel).sa;
    end
end

ll=0;
if all(~isnan(r.traj.like.xc))
    ll = 1;
end

% Upper levels
for j = 1:(l+ll)-(1+ll)

    % Subplots
    subplot(l+ll,1,j);

    if j==(l+ll)-(1+ll)
        mcolour = mu2_colour;
    else
        mcolour =  mu3_colour;
    end

    if plotsd == true
        upperprior = mu_0(l-j+1) +sqrt(sa_0(l-j+1));
        lowerprior = mu_0(l-j+1) -sqrt(sa_0(l-j+1));
        upper = [upperprior; mu(:,l-j+1)+sqrt(sa(:,l-j+1))];
        lower = [lowerprior; mu(:,l-j+1)-sqrt(sa(:,l-j+1))];
    
        plot(0, upperprior, 'o','Color',mcolour{2}, 'LineWidth', 1);
        hold all;
        plot(0, lowerprior, 'o','Color',mcolour{2}, 'LineWidth', 1);
        fill([ts, fliplr(ts)], [(upper)', fliplr((lower)')], ...
             mcolour{2}, 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
    end
    plot(ts, [mu_0(l-j+1); mu(:,l-j+1)], 'Color', mcolour{2}, 'LineWidth', 2);
    hold all;
    plot(0, mu_0(l-j+1), 'o','Color',mcolour{2}, 'LineWidth', 2); % prior
    xlim([0 ts(end)]);
    title(['Posterior expectation of x_' num2str(l-j+1)], 'FontWeight', 'bold');
    ylabel(['\mu_', num2str(l-j+1)]);
    set(gca,'FontSize',fontsize)
end

% Input level
subplot(l+ll,1,l);

try
    plot(ts, [tapas_sgm(mu_0(2), 1); tapas_sgm(mu(:,2), 1)], 'Color', mu2_colour{2}, 'LineWidth', 2);
    hold all;
    plot(0, tapas_sgm(mu_0(2), 1), 'o','Color',mu2_colour{2}, 'LineWidth', 2); % prior
    plot(ts, [mu_0(1); mu(:,1)], 'Color', mu1_colour{2}, 'LineWidth', 1);
    plot(0, mu_0(1), 'o','Color',mu1_colour{2},  'LineWidth', 1); % prior
catch
    plot(ts, [muhat(:,1); nan], 'Color', mu2_colour{2}, 'LineWidth', 2);
    hold all;
    plot(0, muhat(1,1), 'o','Color',mu2_colour{2}, 'LineWidth', 2); % prior
    plot(ts, [mu_0(1); mu(:,1)], 'Color', mu1_colour{2}, 'LineWidth', 1);
    plot(0, mu_0(1), 'o','Color',mu1_colour{2},  'LineWidth', 1); % prior
%     plot(ts, [mu_0(1); mu(:,1)], 'Color', mu2_colour{2},  'LineWidth', 1);
%     hold all;
%     plot(0, mu_0(1), 'o','Color',mu2_colour{2},  'LineWidth', 1); % prior
end
plot(ts(2:end), r.u(:,1), '.', 'Color', input_colour{2}); % inputs
%plot(ts(2:end), wt(:,1), 'k') % implied learning rate 
if (ploty == true) && ~isempty(find(strcmp(fieldnames(r),'y'))) && ~isempty(r.y)
    if ~isempty(find(strcmp(fieldnames(r),'c_sim'))) && strcmp(r.c_sim.obs_model,'tapas_beta_obs')
        y = r.y(:,1);
    else
        y = r.y(:,1) -0.5; y = 1.16 *y; y = y +0.5; % stretch
        if ~isempty(find(strcmp(fieldnames(r),'irr')))
            y(r.irr) = NaN; % weed out irregular responses
            if sum(isnan(y))<0.5*length(y)
                plot(ts(r.irr),  1.08.*ones([1 length(r.irr)]), 'x', 'Color', [1 0.7 0], 'Markersize', 11, 'LineWidth', 2); % irregular responses
                plot(ts(r.irr), -0.08.*ones([1 length(r.irr)]), 'x', 'Color', [1 0.7 0], 'Markersize', 11, 'LineWidth', 2); % irregular responses
            end
        end
    end
    plot(ts(2:end), y, '.', 'Color', [1 0.7 0]); % responses
%     title(['Response y (orange), input u (green), learning rate (fine black), posterior expectation of input s(\mu_2) (blue), ', ...
    title(['Response y (orange), input u (' input_colour{1} '), posterior expectation of input s(\mu_2) (' mu2_colour{1} '), ', ...
           '\mu_1 (' mu1_colour{1} ')'], ...
      'FontWeight', 'bold');
           %', for ',...
           %'\rho=', num2str(rho(2:end)), ', '...
           %'\kappa=', num2str(ka(2:end)), ', '...
           %'\omega=', num2str(om(2:end))], ...
    ylabel('y, u, s(\mu_2), \mu_1');
    axis([0 ts(end) -0.15 1.15]);
else
%     title(['Response y (orange), input u (green), learning rate (fine black), posterior expectation of input s(\mu_2) (blue), ', ...
    title(['Response y (orange), input u (' input_colour{1} '), posterior expectation of input s(\mu_2) (' mu2_colour{1} '), ', ...
           '\mu_1 (' mu1_colour{1} ')'], ...
      'FontWeight', 'bold');
           %', for ',...
           %'\rho=', num2str(rho(2:end)), ', '...
           %'\kappa=', num2str(ka(2:end)), ', '...
           %'\omega=', num2str(om(2:end))], ...
    ylabel('u, s(\mu_2), \mu_1');
    axis([0 ts(end) -0.1 1.1]);
end
plot(ts(2:end), 0.5, 'k');
xlabel('Trial number');
set(gca,'FontSize',fontsize)
hold off;

% Joint models
if ll
    subplot(l+ll,1,l+ll);
    
    xc = r.traj.like.xc;
    vj_mu = r.traj.like.vj_mu;
    plot(ts, [NaN; xc(:,1)], 'm', 'LineWidth', 2);
    hold all;
    plot(ts, [NaN; vj_mu(:,1)], 'k--', 'LineWidth', 2);
    
    title(['posterior expectation of input xc (magenta) and the joint prior vj_mu (black) '], ...
      'FontWeight', 'bold');
    axis([0 ts(end) -0.1 1.1]);
    
    plot(ts(2:end), 0.5, 'k');
    xlabel('Trial number');
    hold off;
end

if 0
   figure
   hold all
   plot(bopars.traj.PR.mu0,'r')
   plot(bopars.traj.AL.mu0,'b')
   plot(bopars.traj.PL.mu0,'g')
   axis([0 length(bopars.traj.PR.mu0) -0.1 1.1]);
   %legend({'PR','AL','PL'})
   hold off
end