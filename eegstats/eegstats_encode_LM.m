function [varargout] = eegstats_encode_LM(varargin)
% Runs Linear Model (now using fitlm) on one or a number of data samples (s)

if isempty(varargin)
    % assume Condor
    load input.mat;
    Y=struct;
    for s = 1:size(tempY,2)
        temp = array2table(tempY(:,s),'VariableNames',{'data'});
        Y(s).dtab = horzcat(dtab,temp);
    end
else
    % assume S,Y inputs
    S=varargin{1};
    Y=varargin{2};
end

LM=length(S.encode.model);
MC=length(S.encode.model_compare);
disp(['number of models:' num2str(LM)]);
disp(['number of model comparisons:' num2str(MC)]);

M = struct;
failed_s = cell(1,LM);
for s = 1:length(Y)

    % for each model
    model=struct;
    for i = 1:LM

        % fit model
        try
            if ismember('train',Y(s).dtab.Properties.VariableNames)
                train_data = Y(s).dtab(find(double(Y(s).dtab.train)),:);
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
            lm = fitlm(train_data,S.encode.model{i}); % Replaced fitlme with fitlm
            last_s_worked=s;
        catch
            failed_s{i} = [failed_s{i} s];
            if isfield(Y(s).dtab,'train')
                train_data = Y(last_s_worked).dtab(find(double(Y(last_s_worked).dtab.train)),:);
            else
                train_data = Y(last_s_worked).dtab;
            end
            if S.encode.zscore
                try
                    mustBeNumericOrLogical(train_data.data)
                    train_data.data = applyz(train_data.data);
                catch
                    disp(['EEG data not numeric or logical for LM ' num2str(LM) ', sample ' num2str(s)])
                end
            end
            lm = fitlm(train_data,S.encode.model{i}); % Replaced fitlme with fitlm
        end
        model(i).lm = lm; % Store the linear model (lm)

        % Outputs common to all samples
        if s == 1
            M.model(i).samples(1).def = char(lm.Formula); % Formula from fitlm
            M.model(i).samples(1).pred = lm.PredictorNames;
            %M.model(i).samples(1).fixeddesign = lm.Formula.Terms; % Design matrix (no random effects)
            % Commented out random effects related outputs
            % M.model(i).samples(1).randomdesign = designMatrix(lmm,'Random');
            M.model(i).samples(1).CoefficientNames = lm.CoefficientNames;
            % M.model(i).samples(1).RandomNames = Rn;
            M.model(i).samples(1).failed_s = failed_s{i};
        end

        % Outputs for each sample
        M.model(i).samples(s).b = double(lm.Coefficients.Estimate); % Coefficients from fitlm
        M.model(i).samples(s).t = double(lm.Coefficients.tStat); % Coefficients from fitlm
        M.model(i).samples(s).se = double(lm.Coefficients.SE); % Coefficients from fitlm
        % M.model(i).samples(s).r = R; % Random effects - not applicable in fitlm
        % M.model(i).samples(s).mse = lmm.MSE; % MSE not directly available in fitlm
        M.model(i).samples(s).logl = lm.LogLikelihood; % LogLikelihood (if available in fitlm)
        M.model(i).samples(s).r2_ord = lm.Rsquared.Ordinary;
        M.model(i).samples(s).r2_adj = lm.Rsquared.Adjusted;
        M.model(i).samples(s).coeffcov = lm.CoefficientCovariance;
        % M.model(i).samples(s).psi = covarianceParameters(lmm); % Covariance parameters not applicable

        % Contrasts (fitlm doesn't support contrasts directly, so leaving it as is)
        if strcmp(S.encode.lmm.contrasts(1).Term,'anova') 
            anovastats = anova(lm,'component');
            if length(anovastats.Properties.RowNames)>1
                M.model(i).samples(s).Term = anovastats.Properties.RowNames(1:end-1);
                M.model(i).samples(s).DF = [anovastats.DF(1:end-1),anovastats.DF(end)*ones(length(anovastats.Properties.RowNames)-1,1)];
                M.model(i).samples(s).F = anovastats.F(1:end-1);
                M.model(i).samples(s).p = anovastats.pValue(1:end-1);
            end
        else
            for h = 1:length(S.encode.lmm.contrasts)
                [pValue,FStat,DF1,DF2] = coefTest(lm,S.encode.lmm.contrasts(h).H);
                M.model(i).samples(s).Term = S.encode.lmm.contrasts(h).Term;
                M.model(i).samples(s).DF(h,:) = [DF1,DF2];
                M.model(i).samples(s).F(h,1) = FStat;
                M.model(i).samples(s).p(h,1) = pValue;
            end
        end

        % Residuals
        resid = lm.Residuals.Raw; % Residuals from fitlm (standardized not directly available)
        try
            if length(resid) >= 50
                M.model(i).samples(s).hnorm = kstest(resid);
            else
                M.model(i).samples(s).hnorm = NaN;
            end
            if length(resid) < 50
                M.model(i).samples(s).swtest = swtest(resid, 0.05);
            else
                M.model(i).samples(s).swtest = NaN;
            end
            M.model(i).samples(s).skew = skewness(resid) / (std(resid) / sqrt(length(resid)));
            M.model(i).samples(s).kurt = kurtosis(resid) / (std(resid) / sqrt(length(resid)));
        catch
            M.model(i).samples(s).hnorm = NaN;
            M.model(i).samples(s).swtest = NaN;
            M.model(i).samples(s).skew = NaN;
            M.model(i).samples(s).kurt = NaN;
        end
        if S.encode.save_residuals
            M.model(i).samples(s).resid = resid;
        end

        % Fitted values
        if S.encode.save_fitted
            ftd = lm.Fitted;
            M.model(i).samples(s).fitted = ftd;
        end

        % Input
        if S.encode.save_input
            input = Y(s).dtab.data;
            M.model(i).samples(s).input = input;
        end
    end

    % Compare models (model comparison logic using fitlm)
    for i = 1:MC
        modi = S.encode.model_compare{i};
        M.model_comp(i).samples(s).pair = modi;
        try
            results = compare(model(modi(1)).lm,model(modi(2)).lm);
            M.model_comp(i).samples(s).pval = results.pValue(2);
            M.model_comp(i).samples(s).LRStat = results.LRStat(2);
            M.model_comp(i).samples(s).deltaDF = results.deltaDF;
        catch
           M.model_comp(i).samples(s).pval = nan;
           M.model_comp(i).samples(s).LRStat = nan;
           M.model_comp(i).samples(s).deltaDF = [nan nan];
        end
    end
end

if isempty(varargin)
    % Condor output
    save('output.mat','M','S','chunk_info');
else
    varargout = {M};
end

function data = applyz(data)

data_mean = nanmean(data);
data_std = nanstd(data);
if data_std > 0
    data = (data - data_mean) / data_std;
end


% Helper function for the Shapiro-Wilk test
function [H, pValue, W] = swtest(x, alpha)
%SWTEST Shapiro-Wilk parametric hypothesis test of composite normality.
%   [H, pValue, SWstatistic] = SWTEST(X, ALPHA) performs the
%   Shapiro-Wilk test to determine if the null hypothesis of
%   composite normality is a reasonable assumption regarding the
%   population distribution of a random sample X. The desired significance 
%   level, ALPHA, is an optional scalar input (default = 0.05).
%
%   The Shapiro-Wilk and Shapiro-Francia null hypothesis is: 
%   "X is normal with unspecified mean and variance."
%
%   This is an omnibus test, and is generally considered relatively
%   powerful against a variety of alternatives.
%   Shapiro-Wilk test is better than the Shapiro-Francia test for
%   Platykurtic sample. Conversely, Shapiro-Francia test is better than the
%   Shapiro-Wilk test for Leptokurtic samples.
%
%   When the series 'X' is Leptokurtic, SWTEST performs the Shapiro-Francia
%   test, else (series 'X' is Platykurtic) SWTEST performs the
%   Shapiro-Wilk test.
% 
%    [H, pValue, SWstatistic] = SWTEST(X, ALPHA)
%
% Inputs:
%   X - a vector of deviates from an unknown distribution. The observation
%     number must exceed 3 and less than 5000.
%
% Optional inputs:
%   ALPHA - The significance level for the test (default = 0.05).
%  
% Outputs:
%  SWstatistic - The test statistic (non normalized).
%
%   pValue - is the p-value, or the probability of observing the given
%     result by chance given that the null hypothesis is true. Small values
%     of pValue cast doubt on the validity of the null hypothesis.
%
%     H = 0 => Do not reject the null hypothesis at significance level ALPHA.
%     H = 1 => Reject the null hypothesis at significance level ALPHA.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Copyright (c) 17 March 2009 by Ahmed Ben Sa�da          %
%                 Department of Finance, IHEC Sousse - Tunisia           %
%                       Email: ahmedbensaida@yahoo.com                   %
%                    $ Revision 3.0 $ Date: 18 Juin 2014 $               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% References:
%
% - Royston P. "Remark AS R94", Applied Statistics (1995), Vol. 44,
%   No. 4, pp. 547-551.
%   AS R94 -- calculates Shapiro-Wilk normality test and P-value
%   for sample sizes 3 <= n <= 5000. Handles censored or uncensored data.
%   Corrects AS 181, which was found to be inaccurate for n > 50.
%   Subroutine can be found at: http://lib.stat.cmu.edu/apstat/R94
%
% - Royston P. "A pocket-calculator algorithm for the Shapiro-Francia test
%   for non-normality: An application to medicine", Statistics in Medecine
%   (1993a), Vol. 12, pp. 181-184.
%
% - Royston P. "A Toolkit for Testing Non-Normality in Complete and
%   Censored Samples", Journal of the Royal Statistical Society Series D
%   (1993b), Vol. 42, No. 1, pp. 37-43.
%
% - Royston P. "Approximating the Shapiro-Wilk W-test for non-normality",
%   Statistics and Computing (1992), Vol. 2, pp. 117-119.
%
% - Royston P. "An Extension of Shapiro and Wilk's W Test for Normality
%   to Large Samples", Journal of the Royal Statistical Society Series C
%   (1982a), Vol. 31, No. 2, pp. 115-124.
%
%
% Ensure the sample data is a VECTOR.
%
if numel(x) == length(x)
    x  =  x(:);               % Ensure a column vector.
else
    error(' Input sample ''X'' must be a vector.');
end
%
% Remove missing observations indicated by NaN's and check sample size.
%
x  =  x(~isnan(x));
if length(x) < 3
   error(' Sample vector ''X'' must have at least 3 valid observations.');
end
if length(x) > 5000
    warning('Shapiro-Wilk test might be inaccurate due to large sample size ( > 5000).');
end
%
% Ensure the significance level, ALPHA, is a 
% scalar, and set default if necessary.
%
if (nargin >= 2) && ~isempty(alpha)
   if ~isscalar(alpha)
      error(' Significance level ''Alpha'' must be a scalar.');
   end
   if (alpha <= 0 || alpha >= 1)
      error(' Significance level ''Alpha'' must be between 0 and 1.'); 
   end
else
   alpha  =  0.05;
end
% First, calculate the a's for weights as a function of the m's
% See Royston (1992, p. 117) and Royston (1993b, p. 38) for details
% in the approximation.
x       =   sort(x); % Sort the vector X in ascending order.
n       =   length(x);
mtilde  =   norminv(((1:n)' - 3/8) / (n + 1/4));
weights =   zeros(n,1); % Preallocate the weights.
if kurtosis(x) > 3
    
    % The Shapiro-Francia test is better for leptokurtic samples.
    
    weights =   1/sqrt(mtilde'*mtilde) * mtilde;
    %
    % The Shapiro-Francia statistic W' is calculated to avoid excessive
    % rounding errors for W' close to 1 (a potential problem in very
    % large samples).
    %
    W   =   (weights' * x)^2 / ((x - mean(x))' * (x - mean(x)));
    % Royston (1993a, p. 183):
    nu      =   log(n);
    u1      =   log(nu) - nu;
    u2      =   log(nu) + 2/nu;
    mu      =   -1.2725 + (1.0521 * u1);
    sigma   =   1.0308 - (0.26758 * u2);
    newSFstatistic  =   log(1 - W);
    %
    % Compute the normalized Shapiro-Francia statistic and its p-value.
    %
    NormalSFstatistic =   (newSFstatistic - mu) / sigma;
    
    % Computes the p-value, Royston (1993a, p. 183).
    pValue   =   1 - normcdf(NormalSFstatistic, 0, 1);
    
else
    
    % The Shapiro-Wilk test is better for platykurtic samples.
    c    =   1/sqrt(mtilde'*mtilde) * mtilde;
    u    =   1/sqrt(n);
    % Royston (1992, p. 117) and Royston (1993b, p. 38):
    PolyCoef_1   =   [-2.706056 , 4.434685 , -2.071190 , -0.147981 , 0.221157 , c(n)];
    PolyCoef_2   =   [-3.582633 , 5.682633 , -1.752461 , -0.293762 , 0.042981 , c(n-1)];
    % Royston (1992, p. 118) and Royston (1993b, p. 40, Table 1)
    PolyCoef_3   =   [-0.0006714 , 0.0250540 , -0.39978 , 0.54400];
    PolyCoef_4   =   [-0.0020322 , 0.0627670 , -0.77857 , 1.38220];
    PolyCoef_5   =   [0.00389150 , -0.083751 , -0.31082 , -1.5861];
    PolyCoef_6   =   [0.00303020 , -0.082676 , -0.48030];
    PolyCoef_7   =   [0.459 , -2.273];
    weights(n)   =   polyval(PolyCoef_1 , u);
    weights(1)   =   -weights(n);
    
    if n > 5
        weights(n-1) =   polyval(PolyCoef_2 , u);
        weights(2)   =   -weights(n-1);
    
        count  =   3;
        phi    =   (mtilde'*mtilde - 2 * mtilde(n)^2 - 2 * mtilde(n-1)^2) / ...
                (1 - 2 * weights(n)^2 - 2 * weights(n-1)^2);
    else
        count  =   2;
        phi    =   (mtilde'*mtilde - 2 * mtilde(n)^2) / ...
                (1 - 2 * weights(n)^2);
    end
        
    % Special attention when n = 3 (this is a special case).
    if n == 3
        % Royston (1992, p. 117)
        weights(1)  =   1/sqrt(2);
        weights(n)  =   -weights(1);
        phi = 1;
    end
    %
    % The vector 'WEIGHTS' obtained next corresponds to the same coefficients
    % listed by Shapiro-Wilk in their original test for small samples.
    %
    weights(count : n-count+1)  =  mtilde(count : n-count+1) / sqrt(phi);
    %
    % The Shapiro-Wilk statistic W is calculated to avoid excessive rounding
    % errors for W close to 1 (a potential problem in very large samples).
    %
    W   =   (weights' * x) ^2 / ((x - mean(x))' * (x - mean(x)));
    %
    % Calculate the normalized W and its significance level (exact for
    % n = 3). Royston (1992, p. 118) and Royston (1993b, p. 40, Table 1).
    %
    newn    =   log(n);
    if (n >= 4) && (n <= 11)
    
        mu      =   polyval(PolyCoef_3 , n);
        sigma   =   exp(polyval(PolyCoef_4 , n));    
        gam     =   polyval(PolyCoef_7 , n);
    
        newSWstatistic  =   -log(gam-log(1-W));
    
    elseif n > 11
    
        mu      =   polyval(PolyCoef_5 , newn);
        sigma   =   exp(polyval(PolyCoef_6 , newn));
    
        newSWstatistic  =   log(1 - W);
    
    elseif n == 3
        mu      =   0;
        sigma   =   1;
        newSWstatistic  =   0;
    end
    %
    % Compute the normalized Shapiro-Wilk statistic and its p-value.
    %
    NormalSWstatistic   =   (newSWstatistic - mu) / sigma;
    
    % NormalSWstatistic is referred to the upper tail of N(0,1),
    % Royston (1992, p. 119).
    pValue       =   1 - normcdf(NormalSWstatistic, 0, 1);
    
    % Special attention when n = 3 (this is a special case).
    if n == 3
        pValue  =   6/pi * (asin(sqrt(W)) - asin(sqrt(3/4)));
        % Royston (1982a, p. 121)
    end
    
end
%
% To maintain consistency with existing Statistics Toolbox hypothesis
% tests, returning 'H = 0' implies that we 'Do not reject the null 
% hypothesis at the significance level of alpha' and 'H = 1' implies 
% that we 'Reject the null hypothesis at significance level of alpha.'
%
H  = (alpha >= pValue);