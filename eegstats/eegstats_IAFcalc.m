function D=eegstats_IAFcalc(S,D)
% inputs: S is a settings structure (see example script)
% outputs: D is the prepared EEG data and predictors as a struct

if nargin<2
    D=struct;
end

% change to the input directory
eval(sprintf('%s', ['cd(''' S.prep.path.inputs.eeg ''')']));

%% find data using getfilelist.m function
GFL=S;
GFL.func = 'prep';
GFL.path.file = S.prep.path.inputs.eeg;
GFL.path.datfile = S.prep.path.inputs.subjects;
if isfield(S.prep.select,'prefixes') && ~isempty(S.prep.select.prefixes); GFL.(GFL.func).load.prefix = S.prep.select.prefixes; end
if isfield(S.prep.select,'suffixes') && ~isempty(S.prep.select.suffixes); GFL.(GFL.func).load.suffix = S.prep.select.suffixes; end
GFL = getfilelist(GFL);
if isfield(GFL.(GFL.func),'designtab')
    designmat=GFL.(GFL.func).designtab;
else
    designmat=cell2table(GFL.(GFL.func).designmat(2:end,:));
    designmat.Properties.VariableNames = GFL.(GFL.func).designmat(1,:);
end
filelist = designmat.file;
subjects = GFL.(GFL.func).select.subjects;
if isempty(filelist)
    error('No files found to import!\n');
end

% if regularised values exist already, load them
if S.prep.calc.eeg.TF.iaf.plot_reg
     loaded_dtab = readtable(fullfile(S.prep.path.outputs,S.prep.sname));
end

% run though all EEG files in a loop
dtab=table;
dtab.ID = subjects';
for d = 1:length(subjects)
    
    % FIND THE FILES FOR THIS SUBJECT
    subfiles = filelist(find(not(cellfun('isempty', strfind(filelist,subjects{d})))));
    
    if isempty(subfiles)
        continue
    end
    
    disp(['subject ' num2str(d) '/' num2str(length(subjects))])  

    % CYCLE THROUGH EACH FILE FOR THIS SUBJECT
    % add data to D structure which has a common format for any data type
    for f = 1:length(subfiles)

        if S.prep.calc.eeg.TF.iaf.calc_PSD

            filename = subfiles{f};
            loadsuffix = S.prep.select.suffixes{1};
            fprintf('\nImporting %s.\n\n', filename);
    
            switch S.prep.fname.ext{:}
                case 'mat'
                    temp = load(filename);
                case 'set'
                    temp = pop_loadset(filename);
            end
            fename = fieldnames(temp(f));
    
    
            % Specify channels to include in IAF calculation
            channels = find(ismember({temp.chanlocs.labels}, S.prep.calc.eeg.TF.iaf.channels));
        end

        % if regularised values exist already, load them
        if S.prep.calc.eeg.TF.iaf.plot_reg
            iaf_reg = loaded_dtab.iaf_reg(find(ismember(loaded_dtab.ID,subjects{d})));
        else 
            iaf_reg = [];
        end
        
        % Calculate IAF using specified channels
        if S.prep.calc.eeg.TF.iaf.calc_PSD
            % calculate PSD
            [iaf, iaf_var, median_psd, fr] = calculateIAF(temp.data, temp.srate, channels, S.prep.calc.eeg.TF.iaf.Nvar, S.prep.calc.eeg.TF.iaf.plot, iaf_reg);
            D(d).median_psd = median_psd;
            D(d).fr = fr;
        else
            % use existing PSD
            [iaf, iaf_var] = calculateIAF(D(d).median_psd, D(d).fr, [], S.prep.calc.eeg.TF.iaf.Nvar, S.prep.calc.eeg.TF.iaf.plot, iaf_reg);
        end

        dtab.iaf(d)=iaf;
        dtab.iaf_var(d)=iaf_var;

        
    end
   

end

% Example vectors
means = dtab.iaf'; % Vector of means
variances = dtab.iaf_var'; % Vector of variances (uncertainties)

% Calculate the group mean
group_mean = mean(means);

% Define a range of lambda values to test
lambda_values = logspace(-2, 2, 100); % From 0.01 to 100, adjust range and resolution as needed
best_lambda = lambda_values(1);
best_p_value = 0;
regularised_means_best = means;

% Function to calculate the Shapiro-Wilk p-value
shapiro_test = @(data) swtest(data);

% Initialize array to store regularised means for each lambda
all_regularised_means = zeros(length(means), length(lambda_values));

% Initialize array to store R-squared values for each lambda
r_squared_values = zeros(1, length(lambda_values));
p_values = zeros(1, length(lambda_values));

% Iterate through lambda values to find the best one
for j = 1:length(lambda_values)
    lambda = lambda_values(j);
    regularised_means = means;
    for i = 1:length(means)
        weight = 1 / variances(i); % Higher weight for lower variance
        regularised_means(i) = (means(i) * weight + lambda * group_mean) / (weight + lambda);
    end
    
    % Store the regularised means
    all_regularised_means(:, j) = regularised_means;
    
    % Compute the R-squared value
    r_squared_values(j) = corr(means', regularised_means')^2;
    
    % Perform the Shapiro-Wilk test for normality
    [h, p_value] = shapiro_test(regularised_means');
    p_values(j) = p_value;
    
    % Update the best lambda if this one is better
    if p_value > best_p_value
        best_lambda = lambda;
        best_p_value = p_value;
        regularised_means_best = regularised_means;
    end
end

dtab.iaf_reg=regularised_means_best';
if S.prep.output.save
    writetable(dtab,fullfile(S.prep.path.outputs,S.prep.sname))
end

% Display the results
disp('Best lambda:')
disp(best_lambda)
disp('P-value of Shapiro-Wilk test with best lambda:')
disp(best_p_value)
disp('Regularised means with best lambda:')
disp(regularised_means_best)

% Create the figure with scatter plot, boxplot, histograms, lambda vs R-squared, and variance vs mean difference
figure;

% Scatter plot of regularised_means_best against original means
subplot(3, 2, 1);
scatter(means, regularised_means_best, 'filled');
xlabel('Original Means');
ylabel('Regularised Means (Best \lambda)');
title(['Scatter Plot (R^2 = ', num2str(corr(means', regularised_means_best')^2), ')']);
grid on;

% Boxplot of R-squared values for each lambda
subplot(3, 2, 2);
boxplot(r_squared_values);
xlabel('\lambda');
ylabel('R^2 Value');
title('R^2 Values for Different \lambda');
grid on;

% Histogram of original means with normal distribution curve
subplot(3, 2, 3);
histogram(means, 'Normalization', 'pdf');
hold on;
x = linspace(min(means), max(means), 100);
y = normpdf(x, mean(means), std(means));
plot(x, y, 'r', 'LineWidth', 2);
hold off;
xlabel('Original Means');
ylabel('Probability Density');
title('Histogram of Original Means');
grid on;

% Histogram of best regularised means with normal distribution curve
subplot(3, 2, 4);
histogram(regularised_means_best, 'Normalization', 'pdf');
hold on;
x = linspace(min(regularised_means_best), max(regularised_means_best), 100);
y = normpdf(x, mean(regularised_means_best), std(regularised_means_best));
plot(x, y, 'r', 'LineWidth', 2);
hold off;
xlabel('Regularised Means (Best \lambda)');
ylabel('Probability Density');
title('Histogram of Best Regularised Means');
grid on;

% Scatter plot of lambda vs R-squared values
subplot(3, 2, 5);
scatter(lambda_values, r_squared_values, 'filled');
hold on;
scatter(best_lambda, corr(means', regularised_means_best')^2, 'filled', 'r');
hold off;
xlabel('\lambda');
ylabel('R^2 Value');
title('Scatter Plot of \lambda vs R^2');
set(gca, 'XScale', 'log'); % Log scale for lambda
grid on;

% Scatter plot of variances vs absolute differences between original and best regularised means
subplot(3, 2, 6);
mean_differences = abs(means - regularised_means_best);
scatter(variances, mean_differences, 'filled');
xlabel('Variances');
ylabel('Absolute Difference (|Original - Regularised|)');
title('Scatter Plot of Variances vs Absolute Mean Differences');
grid on;





end

function [iaf,iaf_var, median_psd, f] = calculateIAF(data, srate, channels, Nvar, plot_psd, iaf_reg)
    % Default to no plotting if not provided
    if nargin < 6
        iaf_reg = [];
    end
    if nargin < 5
        plot_psd = false;
    end

    % if data is 3D, we assume it is EEG chans x time x trials, otherwise
    % we use it is an existing PSD
    if ndims(data)==3

        % Subtract the average event-related response from each trial to obtain induced responses
        mean_response = mean(data, 3);
        induced_responses = data - mean_response;
    
        % Calculate power spectral density (PSD) for each specified channel
        num_channels = length(channels);
        psd_all = [];
    
        % Multitaper parameters
        params.Fs = srate;
        params.tapers = [3 5]; % Time-bandwidth product and number of tapers
        params.fpass = [5 15];
    
        % Loop through each specified channel and calculate PSD using multitaper method
        for ch = 1:num_channels
            channel_data = reshape(induced_responses(channels(ch), :, :), 1, []);
            [S, f] = mtspectrumc(channel_data, params);
            psd_all = [psd_all; S']; % Collect PSD for each channel
        end
        
        % Compute median PSD across specified channels
        median_psd = median(psd_all, 1);
    else
        median_psd = data;
        f = srate;
    end
    
    % Apply Gaussian smoothing to the PSD
    smooth_psd = smoothdata(median_psd, 'gaussian', max(5,round(length(median_psd)*0.01))); % Adjust the window size if needed


    % Define alpha range
    alpha_range = [6 14];
    alpha_indices = find(f >= alpha_range(1) & f <= alpha_range(2));

    % f = f(alpha_indices);
    % median_psd = median_psd(alpha_indices);
    % smooth_psd = smooth_psd(alpha_indices);
    
    
    % Find top N peak frequencies within alpha range for both unsmoothed PSD
    [~,sorted] = sort(median_psd(alpha_indices), 'descend');
    iaf_unsmoothed_top_N = f(alpha_indices(sorted(1:Nvar)));
    iaf_var = var(iaf_unsmoothed_top_N);
    
    % Find peak frequency within alpha range for smoothed PSD
    [~, max_idx_smoothed] = max(smooth_psd(alpha_indices));
    iaf_smoothed = f(alpha_indices(max_idx_smoothed));
    
    % Output the smoothed IAF
    iaf = iaf_smoothed;
    
    % Optional plotting of the PSD
    if plot_psd
        figure;
        plot(f, 10*log10(median_psd), 'b');
        hold on;
        plot(f, 10*log10(smooth_psd), 'r');
        plot(iaf_unsmoothed_top_N, 10*log10(median_psd(alpha_indices(sorted(1:Nvar)))), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
        plot(iaf_smoothed, 10*log10(smooth_psd(alpha_indices(max_idx_smoothed))), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
        xlabel('Frequency (Hz)');
        ylabel('Power Spectral Density (dB)');
        title(['PSD with IAF: ' num2str(iaf_smoothed) ', var: ' num2str(iaf_var)]);
        legend('Unsmoothed PSD', 'Smoothed PSD', 'Unsmoothed IAF peaks', 'IAF');
        xlim(alpha_range);
        if ~isempty(iaf_reg)
            plot(iaf_reg, 10*log10(smooth_psd(alpha_indices(max_idx_smoothed))), 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 8);
            legend('Unsmoothed PSD', 'Smoothed PSD', 'Unsmoothed IAF peaks', 'IAF','IAF reg');
        end
        hold off;
    end
end

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

end



