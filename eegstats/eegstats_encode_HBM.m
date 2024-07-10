function [varargout] = eegstats_encode_HBM(varargin)
% Runs Stan models on one or a number of data samples (s)

if isempty(varargin)
    % assume Condor
    load input.mat;
    Y = struct;
    for s = 1:size(tempY, 2)
        temp = array2table(tempY(:, s), 'VariableNames', {'data'});
        Y(s).dtab = horzcat(dtab, temp);
    end
else
    % assume S, Y inputs
    S = varargin{1};
    Y = varargin{2};
    c = varargin{3};
end

LM = length(S.encode.model);
MC = length(S.encode.model_compare);
disp(['number of models: ' num2str(LM)]);
disp(['number of model comparisons: ' num2str(MC)]);

M = struct;
failed_s = cell(1, LM);
for s = 1:length(Y)

    % for each model
    model = struct;
    for i = 1:LM
        % fit model
        % try
            if ismember('train', Y(s).dtab.Properties.VariableNames)
                train_data = Y(s).dtab(find(double(Y(s).dtab.train)), :);
            else
                train_data = Y(s).dtab;
            end

            if S.encode.zscore
                try
                    mustBeNumericOrLogical(train_data.data)
                    train_data.data = applyz(train_data.data);
                catch
                    disp(['EEG data not numeric or logical for LM ' num2str(LM) ', sample ' num2str(s)])
                end
            end

            % predictors
            if s==1
                [terms, groups, train_pred, categ] = pred_from_formula(S.encode.model{i},Y(s).dtab(find(double(Y(s).dtab.train)),:));
            end

            % Define Stan model code (this should be customized per model)
            if LM~=length(S.encode.stan_code)
                % assume there is a second code for single predictors
                if size(train_pred, 2)==1
                    stan_code = S.encode.stan_code{2};
                else
                    stan_code = S.encode.stan_code{1};
                end
            else
                % assume one code per model
                stan_code = S.encode.stan_code{i};
            end

            % Create temporary directory for Stan model files
            temp_output_dir = fullfile(S.encode.path.temp,['c' num2str(c) 's' num2str(s) 'i' num2str(i)]);
            %temp_output_dir = fullfile(S.encode.path.temp, 'StanTemp');
            mkdir(temp_output_dir);

            % Write Stan model code to temporary file
            model_filename = fullfile(temp_output_dir, 'model.stan');
            fid = fopen(model_filename, 'w');
            if fid == -1
                error('Cannot open file for writing: %s', model_filename);
            end
            fprintf(fid, '%s\n', stan_code{:});

            % Use a while loop to ensure the file is closed before proceeding
            status = -1;
            while status ~= 0
                status = fclose(fid);  % Attempt to close the file
                if status == 0
                    disp('File closed successfully.');
                else
                    disp('Error closing the file, retrying...');
                    pause(0.1); % Short pause before retrying
                end
            end

            % Convert ID to numeric
            [~, ~, subj_numeric] = unique(train_data.ID);

            % Prepare data structure for Stan
            data = struct('N', size(train_data, 1), ...
                          'J', numel(unique(subj_numeric)), ...
                          'K', size(train_pred, 2), ...
                          'subj', subj_numeric, ...
                          'X', train_pred, ...
                          'y', train_data.data);

            % Fit the Stan model
            model(i).stan_model = StanModel('file', model_filename);
            disp('fitting the Stan model')
            fit = model(i).stan_model.sampling('data', data, 'iter', S.encode.hbm.nsamples, 'chains', S.encode.hbm.chains, 'warmup', S.encode.hbm.nsamples/2, 'thin', S.encode.hbm.thin, 'working_dir', temp_output_dir);
            %fit = model(i).stan_model.sampling('data', data);

            % Wait for Stan sampling to complete
            fit.block();

            % Extract relevant data from Stan model
            model(i).fit = fit;
            posterior_samples = fit.extract('permuted', true);

            % Define output structure
            if s == 1
                M.model(i).samples(1).def = char(stan_code);
                M.model(i).samples(1).fixeddesign = train_pred; % Customize if applicable
                M.model(i).samples(1).CoefficientNames = terms; % Customize if applicable
                M.model(i).samples(1).failed_s = failed_s{i};
            end

            % Save parameters and other relevant outputs
            M.model(i).samples(s).group_betas = mean(posterior_samples.mu); % Group-level betas (fixed effects)
            M.model(i).samples(s).cred_intervals(:,1:size(train_pred, 2)) = quantile(posterior_samples.mu, [0.025, 0.975]); % 95% credible intervals
            M.model(i).samples(s).bayesian_p_values = mean(posterior_samples.mu > 0); % Bayesian p-values
            M.model(i).samples(s).individual_betas(:,1:size(train_pred, 2)) = squeeze(mean(posterior_samples.beta, 1)); % Individual betas (random effects)
            M.model(i).samples(s).cred_intervals_individual(:,:,1:size(train_pred, 2)) = quantile(posterior_samples.beta, [0.025, 0.975], 1); % 95% credible intervals for individual betas
            M.model(i).samples(s).logl = mean(posterior_samples.log_lik(:)); % Log-likelihood
            y_pred_mean = mean(posterior_samples.y_pred, 1)';
            M.model(i).samples(s).mse = mean((data.y - y_pred_mean).^2); % MSE
            M.model(i).samples(s).r2_ord = 1 - sum((data.y - y_pred_mean).^2) / sum((data.y - mean(data.y)).^2); % R-squared
            M.model(i).samples(s).coeffcov = cov(posterior_samples.mu); % Coefficient covariance

            % Calculate LOO and standard error
            [loo, loos, pk] = psisloo(posterior_samples.log_lik);
            M.model(i).samples(s).loo = loo;
            M.model(i).samples(s).loo_se = std(loos) * sqrt(length(loos));

            resid = mean(posterior_samples.residuals, 1);
            M.model(i).samples(s).ktest_normality = kstest(resid);
            % residuals
            if S.encode.save_residuals
                M.model(i).samples(s).resid=resid;
            end
            % fitted
            if S.encode.save_fitted
                M.model(i).samples(s).fitted=y_pred_mean;
            end
            % input
            if S.encode.save_input
                M.model(i).samples(s).input=data.y;
            end

            % Clean up temporary directory
            rmdir(temp_output_dir, 's');

        % catch
        %     failed_s{i} = [failed_s{i} s];
        % end
    end

    % Compare models using PSIS-LOO
    for i = 1:MC
        modi = S.encode.model_compare{i};
        M.model_comp(i).samples(s).pair = modi;

        try
            % log_lik1 = model(modi(1)).fit.extract('log_lik', 'permuted', true);
            % log_lik2 = model(modi(2)).fit.extract('log_lik', 'permuted', true);
            % 
            % loo1 = compute_loo(log_lik1);
            % loo2 = compute_loo(log_lik2);
            % 
            % M.model_comp(i).samples(s).loo1 = loo1;
            % M.model_comp(i).samples(s).loo2 = loo2;
            % M.model_comp(i).samples(s).delta_loo = loo1 - loo2;

            % Alternative comparison using LOO and standard error
            s1 = model(modi(1)).fit.extract('permuted', true);
            s2 = model(modi(2)).fit.extract('permuted', true);
            [loo1, loos1, pk1] = psisloo(s1.log_lik);
            [loo2, loos2, pk2] = psisloo(s2.log_lik);
            loodiff = loos1 - loos2;
            M.model_comp(i).samples(s).elpd_diff = sum(loodiff);
            M.model_comp(i).samples(s).elpd_diff_se = std(loodiff) * sqrt(length(loodiff));
        catch
            M.model_comp(i).samples(s).loo1 = NaN;
            M.model_comp(i).samples(s).loo2 = NaN;
            M.model_comp(i).samples(s).delta_loo = NaN;
            M.model_comp(i).samples(s).elpd_diff = NaN;
            M.model_comp(i).samples(s).elpd_diff_se = NaN;
        end
    end

    % % Move all files from temporary to final output directory
    % final_output_dir = S.encode.path.outputs;
    % if ~exist(final_output_dir, 'dir')
    %     mkdir(final_output_dir);
    % end
    % movefile(fullfile(temp_output_dir, '*'), final_output_dir);

end

% for s = 1:length(Y)
%     for i = 1:LM
%         % Clean up temporary directory
%         rmdir(temp_output_dir, 's');
%     end
% end

if isempty(varargin)
    % Condor output
    save('output.mat', 'M', 'S', 'chunk_info');
else
    varargout = {M};
end

end


function [terms, groups, pred, categ] = pred_from_formula(modeldef,dtab)

% parse the modeldef term into predictors
def=regexprep(modeldef, '\s+', ''); % remove spaces
terms = strsplit(def,'+'); % identify terms occuring after a '+'
terms(1)=[];

% % expand
maineffects = {};
interactions = {};
for tm = 1:length(terms) % find interactions
    terms{tm} = regexprep(terms{tm},'+','');
    if any(regexp(terms{tm}, '.[*].'))
        maineffects = strsplit(terms{tm},'*');
        interactions={};
        for i = 2:length(maineffects)
            allcombs = nchoosek(maineffects,i);
            for ac = 1:size(allcombs,1)
                interactions = [interactions {strjoin(allcombs(ac,:),':')}];
            end
        end
        terms{tm}='';
    end
end
terms = [terms maineffects interactions];
terms = unique(terms,'stable');

% if there are interactions, obtain predictors and dummy coding from LMM
if ~isempty(interactions)
    lmm=fitlme(dtab,modeldef); % MUST use reference coding for Bayesreg to recognise cat variables
    pred=designMatrix(lmm,'Fixed');
    terms=lmm.CoefficientNames;
    pnames=lmm.PredictorNames;
else % otherwise don't use LMM, to save time
    pnames = terms;
    for p = 1:length(terms)
        pred(:,p) = double(dtab.(terms{p}));
    end
end
% remove intercept if present
rm = strcmp(terms,'(Intercept)');
pred(:,rm)=[]; terms(rm)=[];
% get grouped and categorical predictors
for pn=1:length(pnames)
    groups{pn} = find(contains(terms,pnames{pn}));
end
for tm = 1:length(terms)
    categ(tm)=all(ismember(unique(pred(:,tm)),[0 1]));
end
categ=find(categ);

% subterms={};
% for tm = 1:length(terms) % find interactions
%     if any(regexp(terms{tm}, '.:.'))
%         subterms{tm} = strsplit(terms{tm},':');
%     else
%         subterms{tm} = terms{tm};
%     end
% end

% generate new interaction variables
% pred=[];
% categ=zeros(1,length(terms));
% for tm = 1:length(subterms) 
%     if iscell(subterms{tm}) % interactions
%         temp=[];
%         subcateg=[];
%         for stm = 1:length(subterms{tm})
%             var=dtab.(subterms{tm}{stm});
%             subcateg(stm)=iscategorical(var);
%             if subcateg(stm)
%                 temp(:,stm)=double(var)-1;
%             else
%                 temp(:,stm)=var;
%             end
%         end
%         pred(:,tm) = prod(temp,2);
%         if all(subcateg)
%             categ(tm)=1;
%         end
%     else
%         var=dtab.(subterms{tm});
%         categ(tm)=iscategorical(var);
%         if categ(tm)
%             pred(:,tm) = double(var)-1;
%         else
%             pred(:,tm) = var;
%         end
%     end
% end
end

function data = applyz(data)
    data_mean = nanmean(data);
    data_std = nanstd(data);
    if data_std > 0
        data = (data - data_mean) / data_std;
    end
end

function loo_val = compute_loo(log_lik)
    log_lik_matrix = log_lik;
    [~, N] = size(log_lik_matrix);
    loo_val = 0;
    for i = 1:N
        log_lik_i = log_lik_matrix(:, i);
        max_log_lik = max(log_lik_i);
        loo_val = loo_val - 2 * log(mean(exp(log_lik_i - max_log_lik))) - 2 * max_log_lik;
    end
end

function [loo,loos,pk] = psisloo(log_lik,varargin)
%PSISLOO Pareto smoothed importance sampling leave-one-out log predictive densities
%   
%  Description
%    [LOO,LOOS,KS] = PSISLOO(LOG_LIK) computes Pareto smoothed importance
%    sampling leave-one-out log predictive densities given posterior
%    samples of the log likelihood terms p(y_i|\theta^s) in LOG_LIK.
%    Returns a sum of the leave-one-out log predictive densities LOO,
%    individual leave-one-out log predictive density terms LOOS and an
%    estimate of Pareto tail indeces KS. The estimates are unreliable if 
%    tail index k>0.7 (see more in the references).
%
%    [LOO,LOOS,KS] = PSISLOO(LOG_LIK,Reff) passes optional
%    arguments for Pareto smoothed importance sampling.
%      Reff - relative MCMC efficiency N_eff/N
%
%  References
%    Aki Vehtari, Andrew Gelman and Jonah Gabry (2017). Practical
%    Bayesian model evaluation using leave-one-out cross-validation
%    and WAIC. Statistics and Computing, 27(5):1413–1432. 
%    doi:10.1007/s11222-016-9696-4. https://arxiv.org/abs/1507.04544
%
%    Aki Vehtari, Andrew Gelman and Jonah Gabry (2017). Pareto
%    smoothed importance sampling. https://arxiv.org/abs/arXiv:1507.02646v5
%
%  Copyright (c) 2015 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

% log raw weights from log_lik
lw=-log_lik;
% compute Pareto smoothed log weights given raw log weights
[lw,pk]=psislw(lw,varargin{:});
% compute
loos=sumlogs(log_lik+lw);
loo=sum(loos);
end

function [lw,kss] = psislw(lw,Reff)
%PSIS Pareto smoothed importance sampling
%
%  Description
%    [LW,K] = PSISLW(LW,Reff) returns log weights LW
%    and Pareto tail indeces K, given log weights and optional arguments:
%      Reff - relative MCMC efficiency N_eff/N
%
%  Reference
%    Aki Vehtari, Daniel Simpson, Andrew Gelman, Yuling Yao, and Jonah
%    Gabry (2024). Pareto smoothed importance sampling. Journal of Machine
%    Learning Research, accepted for publication.
%    https://arxiv.org/abs/arXiv:1507.02646
%
% Copyright (c) 2015-2017 Aki Vehtari, Tuomas Sivula

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.
if size(lw,1)<=1
    error('psislw: more than one log-weight needed');
end
if nargin<2
    Reff=1;
end

for i1=1:size(lw,2)
    % Loop over sets of log weights
    x=lw(:,i1);
    % improve numerical accuracy
    x=x-max(x);
    % Divide log weights into body and right tail
    n=numel(x);
    xs=sort(x);
    xcutoff=xs(end-ceil(min(0.2*n,3*sqrt(n/Reff))));
    if xcutoff<log(realmin)
        % need to stay above realmin
        xcutoff=-700;
    end
    x1=x(x<=xcutoff);
    x2=x(x>xcutoff);
    n2=numel(x2);
    if n2<=4
        % not enough tail samples for gpdfitnew
        qx=x;
        k=Inf;
    else
        % fit generalized Pareto distribution to the right tail samples
        [k,sigma]=gpdfitnew(exp(x2)-exp(xcutoff));
    end
    if k<1/3 || isinf(k)
        % no smoothing if short tail or GPD fit failed
        qx=x;
    else
        x1=x(x<=xcutoff);
        x2=x(x>xcutoff);
        n2=numel(x2);
        [~,x2si]=sort(x2);
        % compute ordered statistic for the fit
        qq=gpinv(([1:n2]-0.5)/n2,k,sigma)+exp(xcutoff);
        % remap back to the original order
        slq=zeros(n2,1);slq(x2si)=log(qq);
        % join body and GPD smoothed tail
        qx=x;qx(x<=xcutoff)=x1;qx(x>xcutoff)=slq;
        % truncate smoothed values to the largest raw weight 0
        lwtrunc=0;
        qx(qx>lwtrunc)=lwtrunc;
    end
    % renormalize weights
    lwx=bsxfun(@minus,qx,sumlogs(qx));
    % return log weights and tail index k
    lw(:,i1)=lwx;
    kss(1,i1)=k;
end
end

function x = gpinv(p,k,sigma)
% Octave compatibility by Tuomas Sivula
    x = NaN(size(p));
    if sigma <= 0
        return
    end
    ok = (p>0) & (p<1);
    if abs(k) < eps
        x(ok) = -log1p(-p(ok));
    else
    x(ok) = expm1(-k * log1p(-p(ok))) ./ k;
    end
    x = sigma*x;
    if ~all(ok)
        x(p==0) = 0;
        x(p==1 & k>=0) = Inf;
        x(p==1 & k<0) = -sigma/k;
    end
end

function [k,sigma,ks,w] = gpdfitnew(x)
%GPDFITNEW Estimate the paramaters for the Generalized Pareto Distribution
%   
%  Description
%    [K,SIGMA] = GPDFITNEW(X) returns empirical Bayes estimate for the
%    parameters k (ksi) and sigma of the two-parameter generalized Parato
%    distribution (GPD) given the data in X.
%    
%    [K,SIGMA,KS,W] = GPDFITNEW(X) returns also the marginal posterior
%    distribution of k in a format of quadrature points KS and
%    quadrature weights W.
%    
%  References
%    Jin Zhang & Michael A. Stephens (2009) A New and Efficient
%    Estimation Method for the Generalized Pareto Distribution,
%    Technometrics, 51:3, 316-325, DOI: 10.1198/tech.2009.08017
%
%    Aki Vehtari, Andrew Gelman and Jonah Gabry (2017). Pareto
%    smoothed importance sampling. https://arxiv.org/abs/1507.02646v5
%
%  Note
%    This function returns a negative of Zhang and Stephens's k,
%    because it is more common parameterisation
%
% Copyright (c) 2015-2017 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

% fit generalized Pareto distribution
n=numel(x);
prior=3;
x=sort(x);
m=30+floor(sqrt(n));
bs=1/x(n)+(1-sqrt(m./([1:m]'-.5)))./prior./x(floor(n./4+0.5));
% loop version matching Zhang and Stephens paper
% w=zeros(m,1);
% L=zeros(m,1);
% ks=zeros(m,1);
% for i=1:m
%     k=-mean(log1p(-bs(i).*x));
%     ks(i)=-k;
%     L(i)=n*(log(bs(i)/k)+k-1);
% end
% for i=1:m
%     w(i)=1/sum(exp(L-L(i)));
% end
% faster vectorized version
% we use negative of Zhang and Stephens's k, because it
% is more common parameterisation
ks=mean(log1p(bsxfun(@times,-bs,x')),2);
L=n*(log(bs./-ks)-ks-1);
w=1./sum(exp(bsxfun(@minus,L,L')))';
sigmas=-ks./bs;

% remove negligible weights
dii=w<eps*10;
w(dii)=[];w=w./sum(w);
bs(dii)=[];
ks(dii)=[];

% posterior mean for b
b=sum(bs.*w);
% estimate for k, note that we return a negative of Zhang and
% Stephens's k, because it is more common parameterisation
k=mean(log1p(-b*x));
ks=mean(log1p(-bs*x'),2);
% estimate for sigma
sigma=-k/b;
% weakly informative prior for k
a=10;
k=k*n/(n+a)+a*0.5/(n+a);
ks=ks*n/(n+a)+a*0.5/(n+a);
end

function y = sumlogs(x,dim)
%SUMLOGS Sum of vector where numbers are represented by their logarithms.
%
%  Description
%    Y=SUMLOGS(X) Computes Y=log(sum(exp(X))) in such a fashion that
%    it works even when elements have large magnitude.
%
%    Y=SUMLOGS(X,DIM) sums along the dimension DIM. 
%
%  Copyright (c) 2013 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

if nargin<2
  dim=find(size(x)>1,1);
end
maxx=max(x(:));
y=maxx+log(sum(exp(x-maxx),dim));
end