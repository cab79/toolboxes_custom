function eegstats_cluster_between_subjects(S, varargin)
    dbstop if error

    %% Load D structure if not provided in varargin
    if isempty(varargin)
        try
            load(fullfile(S.clusterstats.path.inputs, 'D.mat'), 'D')
        catch
            error('Supply a D structure or path/file of the data')
        end
    else
        D = varargin{1};
    end

    % Set paths
    S.path.code = {};
    set_paths(S);

    if length(D) > 1
        save_pref = 'subject_';
    else
        save_pref = 'group_';
    end

    % get model and contrast indices
    if isempty(S.clusterstats.model.index)
        S.clusterstats.model.index = 1:length(D.model);
    end
    if isempty(S.clusterstats.model.contrast) && isfield(D.model(1),'con')
        for i = S.clusterstats.model.index
            S.clusterstats.model.contrast{i} = 1:length(D.model(i).con);
        end
    end

    % Participant IDs
    if isfield(D, 'ID')
        ID = D.ID;
    else
        error('D structure must contain a cell array of subject IDs in D.ID');
    end

    % Results Table
    if ischar(S.clusterstats.path.pred) % If path is provided, load predictor data
        results = readtable(S.clusterstats.path.pred);
    elseif istable(S.clusterstats.path.pred) % If MATLAB table is provided
        results = S.clusterstats.path.pred;
    else
        error('S.clusterstats.path.pred must be either an Excel file path or a MATLAB table');
    end

    if ~isequal(results.ID, ID)
        error('IDs are in the incorrect order');
    end

    npred = width(results);

    % Ensure predictors to be used are specified
    if ~isfield(S.clusterstats, 'predictors') || isempty(S.clusterstats.predictors)
        error('Please specify which predictors to use in S.clusterstats.predictors');
    end

    predictors = S.clusterstats.predictors;
    covariates = []; % Default to no covariates

    % Check if covariates are specified
    if isfield(S.clusterstats, 'covariates') && ~isempty(S.clusterstats.covariates)
        covariates = S.clusterstats.covariates;
    end

    % Outer loop to process each model separately
    for i = S.clusterstats.model.index
        fprintf('Processing model %d...', i);

        % Initialize results, stats_output, and p_values for each model
        model_results = results;
        model_stats_output = [];
        model_p_values = [];

        % Summary Data Handling
        for c = S.clusterstats.model.contrast{i}
            nc = numel(D.model(i).con(c).vox);
            for ci = 1:nc
                % Input summaries (optional)
                if ismember('input', S.clusterstats.summary_data)
                    if any(strcmp(S.clusterstats.summary_types, 'mean'))
                        model_results.(['model' num2str(i) '_con' num2str(c) '_clus' num2str(ci) '_inputmean']) = ...
                            D.model(i).con(c).clus(ci).input_mean_mean';
                    end
                    if any(strcmp(S.clusterstats.summary_types, 'median'))
                        model_results.(['model' num2str(i) '_con' num2str(c) '_clus' num2str(ci) '_inputmedian']) = ...
                            D.model(i).con(c).clus(ci).input_median_median';
                    end
                    if any(strcmp(S.clusterstats.summary_types, 'eig'))
                        model_results.(['model' num2str(i) '_con' num2str(c) '_clus' num2str(ci) '_inputeig']) = ...
                            D.model(i).con(c).clus(ci).input_eig_mean';
                    end
                end

                % Coeff summaries (optional)
                if ismember('coeffs', S.clusterstats.summary_data)
                    for co = 1:length(D.model(i).coeff)
                        inp_ind = find(strcmp(S.clusterstats.summary_coeffs(:, 2), D.model(i).coeff(co).name)); % input index
                        if isempty(inp_ind); continue; end
                        con_name = S.clusterstats.summary_coeffs{inp_ind, 1};
                        c_idx = find(strcmp({D.model(i).con(:).term}, con_name));
                        if isempty(c_idx); continue; end
                        nc = numel(D.model(i).con(c_idx).vox);
                        for ci = 1:nc
                            model_results.([strrep(con_name, ':', '_') '_clus' num2str(ci)]) = D.model(i).con(c_idx).clus(ci).coeff_median;
                        end
                    end
                end
            end
        end

        % index of response variables
        resp_index = (npred+1):width(model_results);
        
        % Clean up variable names to be compatible with MATLAB formulas
        model_results.Properties.VariableNames = matlab.lang.makeValidName(model_results.Properties.VariableNames);


        %% Apply Include index to model_results for each model
        if isfield(S.clusterstats, 'Include_only') && S.clusterstats.Include_only
            include_idx = results.Include == 1;
            model_results = model_results(include_idx, :);
        end

        %% Standardize continuous variables in model_results
        continuous_vars = {};
        for v = [predictors, covariates]
            var_name = v{1};
            if ismember(var_name, model_results.Properties.VariableNames) && isnumeric(model_results.(var_name))
                continuous_vars = [continuous_vars, var_name];
            end
        end
        
        for v = 1:length(continuous_vars)
            var_name = continuous_vars{v};
            if ~contains(var_name, '_std') && ismember(var_name, model_results.Properties.VariableNames)
                model_results.([var_name '_std']) = (model_results.(var_name) - mean(model_results.(var_name))) / std(model_results.(var_name));
            end
        end

        %% Statistical Analysis using fitlm (Parametric) or Residuals for Non-Parametric
        for p = 1:length(predictors)
            predictor_name = predictors{p};
            predictor_vars = [covariates, predictors(p)];
            standardized_predictor_vars = predictor_vars;
            for k = 1:length(predictor_vars)
                if ismember(predictor_vars{k}, continuous_vars)
                    standardized_predictor_vars{k} = [predictor_vars{k}, '_std'];
                end
            end % Standardized predictors where applicable

            % Initialize columns for this predictor
            effect_size_column = [];
            p_value_column = [];

            for v = resp_index
                response_variable = model_results.Properties.VariableNames{v};

                % Standardize response variable if it is numeric
                if isnumeric(model_results.(response_variable)) && ~contains(response_variable, '_std')
                    response_data = model_results.(response_variable);
                    response_data_standardized = (response_data - mean(response_data)) / std(response_data);
                    model_results.([response_variable '_std']) = response_data_standardized;
                    response_variable = [response_variable '_std']; % Use standardized version
                end

                % Construct formula for fitlm
                if ~isempty(covariates)
                    formula_vars = standardized_predictor_vars;
                    % Ensure only standardize continuous covariates
                    % for k = 1:length(covariates)
                    %     if ~ismember(covariates{k}, continuous_vars)
                    %         formula_vars{k} = covariates{k}; % Keep categorical as is
                    %     end
                    % end
                    formula = sprintf('%s ~ %s', response_variable, strjoin(formula_vars, ' + '));
                else
                    formula = sprintf('%s ~ %s', response_variable, predictor_name); % Only the predictor if no covariates
                end
                lm = fitlm(model_results, formula);

                % Extract the effect size (standardized coefficient) for predictor of interest
                predictor_idx = contains(lm.Coefficients.Row, predictor_name);
                if any(predictor_idx)
                    coef = lm.Coefficients.Estimate(predictor_idx);
                    effect_size_column = [effect_size_column; coef];
                end
                p_value = lm.Coefficients.pValue(predictor_idx);
                p_value_column = [p_value_column; p_value];
            end
            model_stats_output = [model_stats_output, effect_size_column]; % Append effect sizes horizontally for each predictor
            model_p_values = [model_p_values, p_value_column]; % Append p-values horizontally for each predictor
        end

        %% Plotting for Each Model
        if ~isempty(model_stats_output)
            figure;
            imagesc(model_stats_output, [-1 1]); colorbar; colormap('jet');
            xticks(1:length(predictors));
            xticklabels(predictors);
            xtickangle(90);
            yticks(1:(width(model_results) - npred));
            yticklabels(model_results.Properties.VariableNames(npred+1:end));
            xlabel('Predictors');
            ylabel('Cluster Variables');
            title(['Predictor Statistics (Model ' num2str(i) ')']);
            hold on;

            % Highlight statistically significant cells (p < 0.05)
            [rows, cols] = find(model_p_values < 0.05);
            for k = 1:length(rows)
                rectangle('Position', [cols(k) - 0.5, rows(k) - 0.5, 1, 1], 'EdgeColor', 'black', 'LineWidth', 1.5);
            end

            %% New Scatter Plots for Statistically Significant Results
            if ~isempty(rows)
                figure; % Create a new figure for scatter plots
                num_significant = length(rows);
                num_cols = ceil(sqrt(num_significant));
                num_rows = ceil(num_significant / num_cols);

                for k = 1:num_significant
                    subplot(num_rows, num_cols, k);
                    predictor_name = predictors{cols(k)};
                    response_variable = model_results.Properties.VariableNames{npred + rows(k)};

                    % Extract data for scatter plot
                    predictor_data = model_results.(predictor_name);
                    response_data = model_results.(response_variable);

                    % Plot scatter depending on the type of predictor (categorical or continuous)
                    if iscategorical(predictor_data) || numel(unique(predictor_data)) < 3
                        % If predictor is binary/categorical
                        boxplot(response_data, predictor_data);
                        xlabel(predictor_name);
                        ylabel(response_variable);
                        title(['Boxplot: ', predictor_name, ' vs ', response_variable]);
                    else
                        % If predictor is continuous
                        scatter(predictor_data, response_data, 'filled');
                        xlabel(predictor_name);
                        ylabel(response_variable);
                        title(['Scatter: ', predictor_name, ' vs ', response_variable]);
                        % Add a linear fit line to visualize the relationship
                        hold on;
                        lm = fitlm(predictor_data, response_data);
                        plot(lm);
                        hold off;
                    end
                end

                % Adjust the layout for better visualisation
                sgtitle(['Scatter Plots for Significant Results (Model ' num2str(i) ')']);
            end
        end


        %% Save Results for Each Model
        if S.clusterstats.save_table
            model_IDresults = [ID, model_results];
            writetable(model_IDresults, fullfile(S.clusterstats.path.inputs, ['between_subjects_summary_model_' num2str(i) '.xls']));
        end
    end

end

function effect_size = compute_effect_size(group1, group2)
    % Computes Cohen's d for effect size between two groups
    pooled_sd = sqrt((std(group1)^2 + std(group2)^2) / 2);
    effect_size = (mean(group1) - mean(group2)) / pooled_sd;
end

function set_paths(S)
    for p = 1:length(S.path.code)
        if S.path.code{p, 1}
            addpath(genpath(S.path.code{p, 2}));
        else
            addpath(S.path.code{p, 2});
        end
    end
end
