function [logp, yhat, res] = bayes_optimal_binary_CAB(r, infStates, ptrans)
% Calculates the log-probabilities of the inputs given the current prediction
%
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2012-2016 Christoph Mathys, TNU, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.

% Inputs
u = r.u(:,1);
n = size(u,1);
u(r.irr) = [];

% Initialize returned log-probabilities as NaNs so that NaN is
% returned for all irregualar trials
logp = NaN(n,1);
yhat = NaN(n,1);
res  = NaN(n,1);

% Weed irregular trials out from inputs and predictions
%

% Predictions
%x = infStates(:,1,1);
if strcmp(r.c_prc.response.rep,'xc')
    x = r.traj.(r.c_obs.model).xchat(:,1);
else
    x = r.traj.(r.c_obs.model).muhat(:,1);
end
x(r.irr) = [];

% Calculate log-probabilities for non-irregular trials
reg = ~ismember(1:n,r.irr);
logp(reg) = u.*log(x) + (1-u).*log(1-x);
yhat(reg) = x;
res(reg) = (u-x)./sqrt(x.*(1-x));

return;
