% SCRIPT for eegstats function: test case on CORE data with LMM encoding
% and BIEM decoding.

% Primary aim of test-case: graded feature selection with BIEM and comparison across single-subject and group-defined clusters. 
% Currently BIEM does not load multi-level LMM data (requires selection of contrasts, etc.), t-tests (although does allow Bayesreg categorical inputs), 
% nor allow feature selection based on group-level encoding stats (although could easily provide T/FDR/TFCE masks from single subjects). 
% Hierarchical LMM (or Bayes) would be ideal so that subject-level coefficients are regularised by the group, leading to greater differentiation 
% of groups after HGF modelling. Currently can only obtain decoded ODD/EPSI from single subjects or group-level fixed effects (both unreliable).

%% Preliminaries
close all
clear all
dbstop if error % optional instruction to stop at a breakpoint if there is an error - useful for debugging
restoredefaultpath
S.parallel.option = 'local'; % options: local (if multiple cores), condor, none 
S.parallel.chunksize = 2500; % for condor only: max number of datapoints (channels*time) to process on each run. Auto re-adjusts into equal portions.

%% Data preparation: eegstats_dataprep
% EEG data inputs and outputs
S.prep.path.main = 'Q:\Projects\CORE\eeg\ana';
S.prep.path.inputs.subjects = 'Q:\Projects\CORE\Participants\Participant_data_age_balanced.xlsx'; % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any group-level covariates of interest
S.prep.path.inputs.eeg = 'C:\Data\CORE\eeg\ana\prep\cleaned\part2';
S.prep.path.inputs.chanlocs = 'Q:\Projects\CORE\eeg\ana\prep\chanlocs.mat';
S.prep.path.inputs.GSNlocs = 'Q:\Projects\CORE\eeg\GSN-HydroCel-128-Flipmap.mat';
S.prep.path.outputs.eeg = 'C:\Data\CORE\eeg\ana\prep\cleaned\part2\eegstats_dataprep'; 
S.prep.path.outputs.stats = 'Q:\Projects\CORE\eeg\ana\stats'; 
S.prep.fname.eeg.parts = {'subject','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'study','subject','session','block','cond'};
S.prep.fname.eeg.ext = {'set'}; 
S.prep.select.subjects = {}; % either a single subject, or leave blank to process all subjects in folder
S.prep.select.sessions = {};
S.prep.select.blocks = {}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.prep.select.conds = {}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.prep.select.suffixes = {'2_merged_cleaned'}; 
S.prep.original_samples = -200:899; % actual samples in loaded EEG data
S.prep.select.samples = -200:899; % select subset of timepoints in each trial, or leave blank for all
S.prep.select.markers = {[1:24]};
% add toolbox paths
if strcmp(S.parallel.option,'condor')
    addpath(genpath(fullfile(S.prep.path.main, 'dependencies')))
    addpath(genpath(fullfile(S.prep.path.main, 'Data')))
else
    run('Q:\MATLAB\projects\CORE\CORE_addpaths')
end
% EEG data operations / transformations / calculations
S.prep.calc.eeg.smooth_samples = 10; % value is the "span" of a moving average smoother
S.prep.calc.eeg.dsample = 20; % best smooth first
S.prep.calc.eeg.flipchan = [3 4]; % rows of S.cond_idx containing trial types to flip channels right to left 
S.prep.calc.eeg.transform = 'arcsinh'; % EEG transform: arcsinh or notrans
S.prep.calc.eeg.zscore = 1; % transform to z-score over trials (unit variance) so that statistical coefficients are normalised
% S.prep.calc.eeg.ndec=8; % trim data to a number of decimal places
% Predictors: factors from EEG markers - ASSUMES THEY ARE NUMERIC IN EEG DATA FILE but this can be updated if needed
S.prep.pred.factor.markers = { % markers in EEG, organised into factors
    {[1 2 9 10 17 18 5 6 13 14 21 22],[3 4 11 12 19 20 7 8 15 16 23 24]} % fac1: odd,stand
    {[1:2:23],[2:2:24]} % fac2: CD
    {[1:8],[9:16],[17:24]} % fac3: CP
    {[1 2 9 10 17 18 3 4 11 12 19 20],[5 6 13 14 21 22 7 8 15 16 23 24]} % fac4: hand: left, right
    };
S.prep.pred.factor.label = {'ODD','CD','CP','HAND'};
S.prep.pred.factor.ordinal = [3]; % name ordinal conditions
% Covariates: data sources
S.prep.path.covariate(1).type = 'HGF'; % type of data input file
S.prep.path.covariate(1).input = 'Q:\Projects\CORE\behaviour\hgf\fitted\D_fit_r1_it9_n30r.mat'; % path and file
S.prep.path.covariate(1).level = 'trial'; % subject, trial
S.prep.path.covariate(1).supp = 'Q:\Projects\CORE\behaviour\hgf\fitted\CORE_fittedparameters_percmodel2_respmodel4_fractrain0_20190211T074650.mat'; % supplementary file if needed
S.prep.path.covariate(1).def = { % definition of covariates within the file
    {'traj','PL','epsi',[0],[]}; % precision-weighted prediction errors
%     {'RT'}
%     {'choice'}
    };
% S.prep.path.covariate(2).type = 'datfile'; % type of data input file
% S.prep.path.covariate(2).input = 'Q:\Projects\CORE\Participants\Participant_data_age_balanced.xlsx'; % path and file
% S.prep.path.covariate(2).level = 'subject'; % subject, trial
% S.prep.path.covariate(2).def = { % definition of covariates within the file
% %     {'age'}; 
%     };

% continuous covariate predictor transforms
S.prep.calc.pred.transform = 'notrans'; % options: arcsinh, rank or notrans
S.prep.calc.pred.zscore = 1;
% PCA on conds/covariates that are multicollinear. PCA results in common variance component and max unique variance components for each
S.prep.calc.pred.PCA_cov = {'ODD','epsi'}; % will include all covariates that include this text for PCA
S.prep.calc.pred.PCA_minvar = 0.9; % minimum fraction of variance explained by sparse PCA relative to standard PCA. Standard PCA values returned if not meeting this criterion.
S.prep.calc.pred.PCA_models = [1:4]; % specifies how many models to run (columns) and how many PCs to estimate for each
% separate EEG/predictor data into training and testing fractions, only for testing decoders
S.prep.traintest.type = 'random'; % method to split data into training and testing. options: 'random' selection of trials (S.trainfrac must be >0); 'cond': select by condition
S.prep.traintest.frac = 1; % fraction of trials to use as training data (taken at random). Remainder used for testing. E.g. training for encoding vs. testing for decoding. Or, MPVA training vs. MVPA testing.
S.prep.traintest.balance_conds =1; % random training trials balanced across conditions
S.prep.traintest.num_runs = 1; % run analysis this many times, each time with a different random choice of training/test data
% output format
S.prep.output.format = 'group'; % options: 'subjects' separately or 'group' data combined
S.prep.output.save = 0; % save data to disk or not
% RUN FUNCTION
[S,D,stats] = eegstats_dataprep(S);

%% Encoding
S.encode.path.inputs = '';
S.encode.type='LMM'; % options: LMM, BRR, MR, etc. - see function 
S.encode.save_residuals=1;
S.encode.save_fitted=1;
S.encode.save_input_images=1;
S.encode.test_collinearity = 0.8; % Pearson's r value threshold to test for collinearity, or 0 to turn off
% LME settings: Specifiying fixed and random factors: https://www.theanalysisfactor.com/specifying-fixed-and-random-factors-in-mixed-models/
S.encode.lme.type = 'hierarchical'; % 'hierarchical' or 'subject'
S.encode.lme.model = {
    ['data ~ 1 + group + group:M1PC1 +',...
    '(1+M1PC1|ID)'];
    ['data ~ 1 + group + group:M2PC1 + group:M2PC2 +',...
    '(1+M2PC1+M2PC2|ID)'];
    ['data ~ 1 + group + group:M3PC1 + group:M3PC2 + group:M3PC3 +',...
    '(1+M3PC1+M3PC2+M3PC3|ID)'];
    ['data ~ 1 + group + group:M4PC1 + group:M4PC2 + group:M4PC3 + group:M4PC4 +',...
    '(1+M4PC1+M4PC2+M4PC3+M4PC4|ID)'];
};
S.encode.lme.model_compare = {
    };
S.encode.lme.fitmethod = 'ML'; % It is recommended that you fit using the maximum likelihood (ML) method prior to model comparison. If you use the restricted maximum likelihood (REML) method, then both models must have the same fixed-effects design matrix.
% RUN FUNCTION
S.encode.sname=datestr(now,30); % savename
[S,D] = eegstats_encode(S);

%% Decoding
S.decode.path.inputs = '';
S.decode.path.inputs.cfglayout = 'Q:\MATLAB\projects\CORE\cosmo\modified functions\GSN92.sfp'; % layouts for cosmo
S.decode.path.inputs.cfgoutput = 'Q:\MATLAB\projects\CORE\cosmo\modified functions\GSN92.lay'; % layouts for cosmo
S.decode.type='biem';
S.decode.prior = 'subject_training'; % options: 'group_training', 'subject_training', 'uniform'
S.decode.groupweights = '';%'stats_grp_MR_all_chan_cond_arcsinh_20180804T144006.mat'; % repeat of MR mismatch, with decoding';%'stats_grp_RR_all_chan_RT_arcsinh_20180802T130052.mat'; % file containing beta and sigma values from a previous encoding run, group averaged.
S.decode.pred = 1; % specify which predictor variable(s) to correlate with input
% RUN FUNCTION
S.decode.sname=datestr(now,30); % savename
[S,D,stats] = eegstats_decode(S);