function eegstats_encode_compile_samples(varargin) % local version
% Compiles samples from chunks and reshapes back to EEG data dimensions

% Format: eegstats_encode_compile_samples(S,C)
% C: structure of size c, containing Y (size s), chunk_info and samples (size s). s = number of samples per chunk.

%% find outputted data
if ~isempty(varargin)
    iscondor = 0;
    S = varargin{1};
    C = varargin{2};
    pth = S.encode.path.outputs;
    pthData = pth;
else
    % assume using condor
    iscondor = 1;
    [S, C] = load_condor_data;
    pth = S.temp.pth;
    pthData = S.temp.pthData;
end

% initiate main output, D
D = struct;

% data info
cidx = find(~cellfun(@isempty, {C(:).chunk_info}));
nD = C(cidx(1)).chunk_info.nD;
ddim = C(cidx(1)).chunk_info.dim; % EEG data dimensions for reshaping
if nD > 1
    save_pref = 'subject_';
else
    save_pref = 'group_';
end

% loop through chunks of samples

for d = 1:nD % subject

    D(d).ID = S.ID{d};

    M = struct; % temporary structure
    si = []; sinz = [];
    
    %% compile samples from chunks
    
    % assume same number of chunks per subject
    n_chunks_d = C(cidx(1)).chunk_info.n_chunks_d;
    
    % chunks per subject - THIS VERSION WORKS FOR BRRs, UNSURE IF WORKS FOR LMMs
    for nc = 1:n_chunks_d 
        
        % chunk index, c
        c = (d - 1) * n_chunks_d + nc;
        
        if c > length(C); continue; end
        
        disp(['compiling output from file ' num2str(c)])
        
        si = [si C(c).chunk_info.sample_index]; % samples
        if isfield(C(c).chunk_info, 'sample_index_nonzero')
            sinz = [sinz C(c).chunk_info.sample_index_nonzero]; % samples
        else
            sinz = si; % samples
        end
    
        if ~isfield(M, 'model')
            M = C(c).M;
        else
            if isfield(C(c).M, 'model')
                for i = 1:length(M.model)
                    M.model(i).samples = [M.model(i).samples, C(c).M.model(i).samples];
                end
            end
            if isfield(C(c).M, 'model_comp')
                for i = 1:length(M.model_comp)
                    M.model_comp(i).samples = [M.model_comp(i).samples, C(c).M.model_comp(i).samples];
                end
            end
        end
        if d == nD; C(c).M = []; end % save memory
    end

    %% reshape to EEG data dimensions
    % loop through chunks of samples
    LM = length(M.model);
    % for each model
    for i = 1:LM
        % items per model
        D(d).model(i).def = M.model(i).samples(1).def;
        D(d).model(i).fixeddesign = M.model(i).samples(1).fixeddesign;
        
        %OLD: D(d).model(i).s = reshape(vertcat([M.model(i).samples(:).mse]'), ddim(1), ddim(2));
        D(d).model(i).s = nan(ddim(1:end-1)');
        D(d).model(i).s(sinz) = vertcat([M.model(i).samples(:).mse]');
        
        % D(d).model(i).logl = reshape(vertcat([M.model(i).samples(:).logl]'), ddim(1), ddim(2));
        D(d).model(i).logl = nan(ddim(1:end-1)');
        D(d).model(i).logl(sinz) = vertcat([M.model(i).samples(:).logl]');
        
        % D(d).model(i).r2_ord = reshape(vertcat([M.model(i).samples(:).r2_ord]'), ddim(1), ddim(2));
        D(d).model(i).r2_ord = nan(ddim(1:end-1)');
        D(d).model(i).r2_ord(sinz) = vertcat([M.model(i).samples(:).r2_ord]');
        
        if ismember(S.encode.method, {'LMM','BRR'})
            % items per model/coefficient
            for ii = 1:length(M.model(i).samples(1).CoefficientNames)
                D(d).model(i).coeff(ii).name = M.model(i).samples(1).CoefficientNames{ii};
                
                D(d).model(i).coeff(ii).b = nan(ddim(1:end-1)');
                D(d).model(i).coeff(ii).b(sinz) = arrayfun(@(X) X.b(ii), M.model(i).samples);
            end
        end

        % HBM-specific fields
        if strcmp(S.encode.method, 'HBM')
            
            % % Group-level betas
            % D(d).model(i).group_betas = nan(ddim(1:end-1)');
            % D(d).model(i).group_betas(sinz) = vertcat([M.model(i).samples(:).group_betas]');

            % items per model/coefficient
            for ii = 1:length(M.model(i).samples(1).CoefficientNames)
                D(d).model(i).coeff(ii).name = M.model(i).samples(1).CoefficientNames{ii};
                
                % Group-level betas
                D(d).model(i).coeff(ii).group_betas = nan([ddim(1:end-1)']);
                D(d).model(i).coeff(ii).group_betas(sinz) = arrayfun(@(X) X.group_betas(ii), M.model(i).samples);
            
                % Credible intervals for group-level betas
                cred_intervals = arrayfun(@(X) X.cred_intervals(:, ii), M.model(i).samples, 'UniformOutput', false);
                cred_intervals = reshape([cred_intervals{:}], 2, []);
                D(d).model(i).coeff(ii).cred_intervals = squeeze(nan([2 ddim(1:end-1)']));
                D(d).model(i).coeff(ii).cred_intervals(:, sinz) = cred_intervals;
            
                % Bayesian p-values
                D(d).model(i).coeff(ii).bayesian_p_values = nan([ddim(1:end-1)']);
                D(d).model(i).coeff(ii).bayesian_p_values(sinz) = arrayfun(@(X) X.bayesian_p_values(ii), M.model(i).samples);
            
                % Individual betas
                individual_betas = arrayfun(@(X) X.individual_betas(:, ii), M.model(i).samples, 'UniformOutput', false);
                individual_betas = reshape([individual_betas{:}], [], length(M.model(i).samples));
                D(d).model(i).coeff(ii).individual_betas = squeeze(nan([ddim(1:end-1)' size(M.model(i).samples(1).individual_betas, 1)]));
                D(d).model(i).coeff(ii).individual_betas(sinz, :) = individual_betas';
            
                % Credible intervals for individual betas
                cred_intervals_individual = arrayfun(@(X) X.cred_intervals_individual(:, :, ii), M.model(i).samples, 'UniformOutput', false);
                cred_intervals_individual = reshape([cred_intervals_individual{:}], 2, [], length(M.model(i).samples));
                D(d).model(i).coeff(ii).cred_intervals_individual = squeeze(nan([2 ddim(1:end-1)' size(M.model(i).samples(1).individual_betas, 1)]));
                D(d).model(i).coeff(ii).cred_intervals_individual(:, sinz, :) = permute(cred_intervals_individual, [1, 3, 2]);
            end




            
            % LOO and standard error
            D(d).model(i).loo = nan(ddim(1:end-1)');
            D(d).model(i).loo(sinz) = vertcat([M.model(i).samples(:).loo]');
            
            D(d).model(i).loo_se = nan(ddim(1:end-1)');
            D(d).model(i).loo_se(sinz) = vertcat([M.model(i).samples(:).loo_se]');
            
            % Normality test results (standardized)
            D(d).model(i).ktest = nan(ddim(1:end-1)');
            D(d).model(i).ktest(sinz) = vertcat([M.model(i).samples(:).ktest_normality]');
        end
        
        % LME models only
        if strcmp(S.encode.method, 'LMM')
            D(d).model(i).pred = M.model(i).samples(1).pred;
            
            D(d).model(i).r2_adj = nan(ddim(1:end-1)');
            D(d).model(i).r2_adj(sinz) = vertcat([M.model(i).samples(:).r2_adj]');
            
            D(d).model(i).randomdesign = M.model(i).samples(1).randomdesign;
            D(d).model(i).randomnames = M.model(i).samples(1).RandomNames;
            
            if isfield(M.model(i).samples(1), 'failed_s'); D(d).model(i).failed_s = M.model(i).samples(1).failed_s; end
            
            for ii = 1:length(M.model(i).samples(1).Term)
                D(d).model(i).con(ii).term = M.model(i).samples(1).Term{ii};
                
                DF = arrayfun(@(X) X.DF(ii,:), M.model(i).samples, 'UniformOutput', 0);
                D(d).model(i).con(ii).DF = nan([2 ddim(1:end-1)']);
                D(d).model(i).con(ii).DF(:, sinz) = vertcat(DF{:})';
                
                D(d).model(i).con(ii).F = nan(ddim(1:end-1)');
                D(d).model(i).con(ii).F(sinz) = arrayfun(@(X) X.F(ii), M.model(i).samples);
                
                D(d).model(i).con(ii).p = nan(ddim(1:end-1)');
                D(d).model(i).con(ii).p(sinz) = arrayfun(@(X) X.p(ii), M.model(i).samples);
            end
            
            % random effects
            try
                D(d).model(i).random = nan([prod(ddim(1:end-1)') height(D(d).model(i).randomnames)]);
                D(d).model(i).random(sinz, :) = vertcat([M.model(i).samples(:).r]');
                D(d).model(i).random = reshape(D(d).model(i).random, [ddim(1:end-1)' height(D(d).model(i).randomnames)]);
            catch % in case dimension are inconsistent, save as cell
                D(d).model(i).random = {M.model(i).samples(:).r};
            end
            
            % covariance
            D(d).model(i).coeffcov = cell(ddim(1:end-1)');
            D(d).model(i).coeffcov(sinz) = arrayfun(@(X) X.coeffcov, M.model(i).samples, 'UniformOutput', 0);
            for ii = 1:length(M.model(i).samples(1).psi)
                D(d).model(i).psi(ii).cov = cell(ddim(1:end-1)');
                D(d).model(i).psi(ii).cov(sinz) = arrayfun(@(X) X.psi{ii}, M.model(i).samples, 'UniformOutput', 0);
            end

            D(d).model(i).ktest = nan(ddim(1:end-1)');
            D(d).model(i).swtest = nan(ddim(1:end-1)');
            D(d).model(i).skew = nan(ddim(1:end-1)');
            D(d).model(i).kurt = nan(ddim(1:end-1)');
            try
                D(d).model(i).ktest(sinz) = vertcat([M.model(i).samples(:).hnorm]');
                D(d).model(i).swtest(sinz) = vertcat([M.model(i).samples(:).swtest]');
                D(d).model(i).skew(sinz) = vertcat([M.model(i).samples(:).skew]');
                D(d).model(i).kurt(sinz) = vertcat([M.model(i).samples(:).kurt]');
            end
        end
        
        % BRR models only
        if strcmp(S.encode.method, 'BRR')
            D(d).model(i).waic = nan(ddim(1:end-1)');
            D(d).model(i).waic(sinz) = vertcat([M.model(i).samples(:).waic]');
            D(d).model(i).ktest = nan(ddim(1:end-1)');
            D(d).model(i).ktest(sinz) = vertcat([M.model(i).samples(:).ktest_normality]');
            
            for ii = 1:length(M.model(i).samples(1).CoefficientNames)
                D(d).model(i).coeff(ii).t = nan(ddim(1:end-1)');
                D(d).model(i).coeff(ii).t(sinz) = arrayfun(@(X) X.t(ii), M.model(i).samples);
                
                D(d).model(i).coeff(ii).p = nan(ddim(1:end-1)');
                D(d).model(i).coeff(ii).p(sinz) = arrayfun(@(X) X.p(ii), M.model(i).samples);
                
                D(d).model(i).coeff(ii).df = nan(ddim(1:end-1)');
                D(d).model(i).coeff(ii).df(sinz) = arrayfun(@(X) X.df, M.model(i).samples);
            end
        end
        
        % filenames for resid, fitted, input
        D(d).model(i).resid_file = fullfile(pthData, [save_pref num2str(d) '_resid_mi' num2str(i) '_' S.encode.sname '.mat']);
        D(d).model(i).fitted_file = fullfile(pthData, [save_pref num2str(d) '_fitted_mi' num2str(i) '_' S.encode.sname '.mat']);
        D(d).model(i).input_file = fullfile(pthData, [save_pref num2str(d) '_input_mi' num2str(i) '_' S.encode.sname '.mat']);
        
        if strcmp(S.encode.memory_option, 'memory')
            if S.encode.save_residuals
                cd(pthData)
                
                resid = nan([prod(ddim(1:end-1)') ddim(end)]);
                resid(sinz, :) = vertcat([M.model(i).samples(:).resid]');
                resid = reshape(resid, [ddim(1:end-1)' ddim(end)]);
                
                disp([save_pref ' saving residuals from model ' num2str(i) '/' num2str(LM)]);
                save(D(d).model(i).resid_file, 'resid', '-v7.3');
                cd(pth)
            end

            if S.encode.save_input
                cd(pthData)
                
                input = nan([prod(ddim(1:end-1)') ddim(end)]);
                input(sinz, :) = vertcat([M.model(i).samples(:).input]');
                input = reshape(input, [ddim(1:end-1)' ddim(end)]);
                
                disp([save_pref ' saving input from model ' num2str(i) '/' num2str(LM)]);
                save(D(d).model(i).input_file, 'input', '-v7.3');
                cd(pth)
            end

            if S.encode.save_fitted
                cd(pthData)
                
                fitted = nan([prod(ddim(1:end-1)') ddim(end)]);
                fitted(sinz, :) = vertcat([M.model(i).samples(:).fitted]');
                fitted = reshape(fitted, [ddim(1:end-1)' ddim(end)]);
                
                disp([save_pref ' saving fitted from model ' num2str(i) '/' num2str(LM)]);
                save(D(d).model(i).fitted_file, 'fitted', '-v7.3');
                cd(pth)
            end
        elseif iscondor
            if S.encode.save_residuals
                resid = nan(ddim(end), length(M.model(i).samples));
                for s = 1:length(M.model(i).samples)
                    temp = M.model(i).samples(s).resid;

                    % z-score
                    resid_mean = nanmean(temp);
                    resid_std = nanstd(temp);
                    resid(:, s) = (temp - resid_mean) / resid_std;

                    % clear memory            
                    M.model(i).samples(s).resid = []; 
                end
                disp([save_pref ' saving residuals from model ' num2str(i) '/' num2str(LM)]);
                cd(pthData)
                save(D(d).model(i).resid_file, 'resid', '-v7.3');
                %zip('resid.zip', D(d).model(i).resid_file); 
                cd(pth)
            end
            if S.encode.save_input
                input = nan(ddim(end), length(M.model(i).samples));
                for s = 1:length(M.model(i).samples)
                    input(:, s) = M.model(i).samples(s).input;
                    % clear memory            
                    M.model(i).samples(s).input = []; 
                end
                disp([save_pref ' saving inputs from model ' num2str(i) '/' num2str(LM)]);
                cd(pthData)
                save(D(d).model(i).input_file, 'input', '-v7.3');
                cd(pth)
            end
            if S.encode.save_fitted
                fitted = nan(ddim(end), length(M.model(i).samples));
                for s = 1:length(M.model(i).samples)
                    fitted(:, s) = M.model(i).samples(s).fitted;
                    % clear memory            
                    M.model(i).samples(s).fitted = []; 
                end
                disp([save_pref ' saving fitted from model ' num2str(i) '/' num2str(LM)]);
                cd(pthData)
                save(D(d).model(i).fitted_file, 'fitted', '-v7.3');
                cd(pth)
            end
        end
    end

    if isfield(M, 'model_comp')
        MC = length(M.model_comp);
        for i = 1:MC
            D(d).model_comp(i).pair = M.model_comp(i).samples(1).pair;
            
            if strcmp(S.encode.method, 'LMM')
                D(d).model_comp(i).pval = nan(ddim(1:end-1)');
                D(d).model_comp(i).pval(sinz) = vertcat([M.model_comp(i).samples(:).pval]');
                
                D(d).model_comp(i).LR = nan(ddim(1:end-1)');
                D(d).model_comp(i).LR(sinz) = vertcat([M.model_comp(i).samples(:).LRStat]');
            elseif strcmp(S.encode.method, 'HBM')
                D(d).model_comp(i).elpd_diff = nan(ddim(1:end-1)');
                D(d).model_comp(i).elpd_diff(sinz) = vertcat([M.model_comp(i).samples(:).elpd_diff]');
                
                D(d).model_comp(i).elpd_diff_se = nan(ddim(1:end-1)');
                D(d).model_comp(i).elpd_diff_se(sinz) = vertcat([M.model_comp(i).samples(:).elpd_diff_se]');
            end
        end
    end
end

%% save results
disp('saving stats')
mkdir(fullfile(S.encode.path.outputs))
save(fullfile(S.encode.path.outputs, [S.encode.sname '.mat']), 'D', 'S', '-v7.3');
if iscondor
    disp('zipping outputs and inputs')
    zip('outputs.zip', 'output*.mat');
    zip('inputs.zip', 'input*.mat');
    disp('deleting individual outputs and inputs')
    delete('output*.mat');
    delete('input*.mat');

    quit
end


function [S,C] = load_condor_data
pth=pwd;
pthData=fullfile(pwd,'Data');

% add toolbox paths
addpath(genpath(fullfile(pth, 'dependencies_supp')))
addpath(genpath(pthData))

% monitor for Condor job submission to complete
eegstats_condor_monitor_outputs;

% get file info
load('output0.mat','S','chunk_info')
S.temp.pth=pth;
S.temp.pthData=pthData;

if 0%S.encode.delete_input
    try
        copyfile input0.mat example_input.mat
        delete('input*.mat');
    end
end

if 0%S.encode.delete_logs
    try
        delete('*.Z');
        delete('*.log');
        delete('*.err');
        delete('*.out');
    end
end
% try
% %     delete(fullfile(pthData,'resid*.zip'));
%     delete('output.zip');
%     delete('input.zip');
% end

nD=chunk_info.nD;
n_chunks_d=chunk_info.n_chunks_d;
index_length = nD*n_chunks_d;
C=struct;
for d = 1:nD % subject
    for nc = 1:n_chunks_d % chunks per subject
        % index
        c = (d-1)*n_chunks_d +nc;
        
        outname = ['output' num2str(c-1) '.mat'];
        if exist(outname,'file')
            disp(['gathering output from file ' num2str(c) '/' num2str(index_length)])
            fileout = load(outname);
            C(fileout.chunk_info.c).M = fileout.M;
            C(fileout.chunk_info.c).chunk_info = fileout.chunk_info;

%             if S.encode.save_input
%                 filein = load(['input' num2str(c-1) '.mat']);
%                 C(filein.chunk_info.c).Y = filein.Y;
%             end
        else
            disp(['output file ' num2str(c) ' does not exist'])
        end
    end
end

function [in,Z] = eegstats_condor_monitor_outputs


nowtime = now;
in = dir('input*.mat');

fin=0;
while fin==0
    Z = dir('output*.mat');
    complete = [Z(:).bytes]>0;
    if length(Z)>=length(in) && all(complete)
        fin=1;
    end
    disp([num2str(sum(complete)) ' complete out of ' num2str(length(complete)) ' generated from ' num2str(length(in)) ' inputs' ])
    pause(300)
end