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
% add toolbox paths
% condor
%     addpath(genpath(fullfile(S.prep.path.main, 'dependencies')))
%     addpath(genpath(fullfile(S.prep.path.main, 'Data')))
% local
    run('C:\Users\cab79\Google Drive\Docs temp\Matlab\CORE\CORE_addpaths')


%% Data preparation: eegstats_dataprep
% EEG data inputs and outputs
S.prep.path.main = 'C:\Data\CORE\eeg\ana';
S.prep.path.inputs.subjects = 'C:\Users\cab79\Google Drive\Docs temp\Matlab\Participant_data_age_balanced.xlsx'; % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any group-level covariates of interest
S.prep.path.inputs.eeg = 'C:\Data\CORE\eeg\ana\prep\cleaned\part2';
S.prep.path.inputs.chanlocs = 'C:\Users\cab79\Google Drive\Docs temp\Matlab\chanlocs.mat';
S.prep.path.inputs.GSNlocs = 'C:\Users\cab79\Google Drive\Docs temp\Matlab\GSN-HydroCel-128-Flipmap.mat';
S.prep.path.outputs = 'C:\Data\CORE\eeg\ana\prep\cleaned\part2\eegstats_dataprep'; 
S.prep.fname.parts = {'subject','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'study','subject','session','block','cond'};
S.prep.fname.ext = {'set'}; 
S.prep.select.subjects = {}; % either a single subject, or leave blank to process all subjects in folder
S.prep.select.sessions = {};
S.prep.select.blocks = {}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.prep.select.conds = {}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.prep.select.suffixes = {'2_merged_cleaned'}; 
S.prep.original_samples = -200:899; % actual samples in loaded EEG data
S.prep.select.samples = -200:899; % select subset of timepoints in each trial, or leave blank for all
S.prep.select.markers = {[1:24]};
S.prep.parallel.option = 'none'; % options: local (if multiple cores), condor, none 
% EEG data operations / transformations / calculations
S.prep.calc.eeg.smooth_samples = 10; % value is the "span" of a moving average smoother
S.prep.calc.eeg.dsample = 20; % best smooth first
S.prep.calc.eeg.flipchan = [5 6 13 14 21 22 7 8 15 16 23 24]; % markers in EEG
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
S.prep.pred.covariate(1).type = 'HGF'; % type of data input file
S.prep.pred.covariate(1).input = 'C:\Users\cab79\Google Drive\Docs temp\Matlab\D_fit_r1_it9_n30r.mat'; % path and file
S.prep.pred.covariate(1).level = 'trial'; % subject, trial
S.prep.pred.covariate(1).supp = 'C:\Users\cab79\Google Drive\Docs temp\Matlab\CORE_fittedparameters_percmodel2_respmodel4_fractrain0_20190211T074650.mat'; % supplementary file if needed
S.prep.pred.covariate(1).def = { % definition of covariates within the file
    {{'traj'},{'PL'},{'epsi'},{[0]},{[]}}; % precision-weighted prediction errors
%     {'RT'}
%     {'choice'}
    };
% S.prep.pred.covariate(2).type = 'datfile'; % type of data input file
% S.prep.pred.covariate(2).input = 'Q:\Projects\CORE\Participants\Participant_data_age_balanced.xlsx'; % path and file
% S.prep.pred.covariate(2).level = 'subject'; % subject, trial
% S.prep.pred.covariate(2).def = { % definition of covariates within the file
% %     {'age'}; 
%     };

% continuous covariate predictor transforms
S.prep.calc.pred.transform = 'notrans'; % options: arcsinh, rank or notrans
S.prep.calc.pred.zscore = 1;
% PCA on conds/covariates that are multicollinear. PCA results in common variance component and max unique variance components for each
S.prep.calc.pred.PCA_cov = {};%{'ODD','epsi'}; % will include all covariates that include this text for PCA
S.prep.calc.pred.PCA_minvar = 0.9; % minimum fraction of variance explained by sparse PCA relative to standard PCA. Standard PCA values returned if not meeting this criterion.
S.prep.calc.pred.PCA_models = [1:4]; % specifies how many models to run (columns) and how many PCs to estimate for each
S.prep.calc.pred.test_collinearity = 0;%0.8; % best to turn on if not doing PCA or prior to deciding whether to do PCA. value is e.g. 0.8.
% separate EEG/predictor data into training and testing fractions, only for testing decoders
S.prep.traintest.type = 'random'; % method to split data into training and testing. options: 'random' selection of trials (S.trainfrac must be >0); 'cond': select by condition
S.prep.traintest.frac = 1; % fraction of trials to use as training data (taken at random). Remainder used for testing. E.g. training for encoding vs. testing for decoding. Or, MPVA training vs. MVPA testing.
S.prep.traintest.balance_conds =1; % random training trials balanced across conditions
S.prep.traintest.num_runs = 1; % run analysis this many times, each time with a different random choice of training/test data
% output format
S.prep.output.format = 'subject'; % options: 'subject' separately or 'group' data combined
S.prep.output.save = 1; % save data to disk or not
% RUN FUNCTION
S.prep.sname=datestr(now,30); % savename
D = eegstats_dataprep(S);

%% Encoding
S=struct;
S.encode.path.inputs = 'C:\Data\CORE\eeg\ana\prep\cleaned\part2\eegstats_dataprep\eegstats_dtab20191230T153645.mat'; % full path and filename of prepped data to load - or feed in D if in workspace
S.encode.path.outputs = 'C:\Data\CORE\eeg\ana\stats';
S.encode.method='BRR'; % options: LMM, BRR - see functions
% LME settings: Specifiying fixed and random factors: https://www.theanalysisfactor.com/specifying-fixed-and-random-factors-in-mixed-models/
S.encode.model = {
%     ['data ~ 1 + M1PC1 + group + group:M1PC1 +',...
%     '(1+M1PC1|ID)'];
%     ['data ~ 1 + M2PC1 + M2PC2 + group + group:M2PC1 + group:M2PC2 +',...
%     '(1+M2PC1+M2PC2|ID)'];
%     ['data ~ 1 + M3PC1 + M3PC2 + M3PC3 + group + group:M3PC1 + group:M3PC2 + group:M3PC3 +',...
%     '(1+M3PC1+M3PC2+M3PC3|ID)'];
%     ['data ~ 1 + M4PC1 + M4PC2 + M4PC3 + M4PC4 + group + group:M4PC1 + group:M4PC2 + group:M4PC3 + group:M4PC4 +',...
%     '(1+M4PC1+M4PC2+M4PC3+M4PC4|ID)'];
%     ['data ~ 1 + group*ODD*CP +',...
%     '(1+ODD*CP|ID)'];
       ['data ~ 1 + ODD*CP'];
};
S.encode.model_compare = {
    };
S.encode.lmm.fitmethod = 'ML'; % It is recommended that you fit using the maximum likelihood (ML) method prior to model comparison. If you use the restricted maximum likelihood (REML) method, then both models must have the same fixed-effects design matrix.
S.encode.lmm.coding = 'effects'; % 'effects', 'reference', 'full'
S.encode.lmm.contrasts.Term = 'anova'; % specify 'anova' for all main effects/interactions
% S.encode.lmm.contrasts(1).Term = 'ODD';
% S.encode.lmm.contrasts(1).H = [];
S.encode.brr.folds = 0;            % number of folds in the traindata. Set to 0 to not conduct predictive cross-validation.
S.encode.brr.model = 't';   % error distribution - string, one of {'gaussian','laplace','t','binomial'}
S.encode.brr.prior = 'HS+';        %- string, one of {'g','ridge','lasso','horseshoe','horseshoe+'}
S.encode.brr.nsamples = 100;   %- number of posterior MCMC samples (Default: 1000)  
S.encode.brr.burnin = 100;     %- number of burnin MCMC samples (Default: 1000)
S.encode.brr.thin = 5;       %- level of thinning (Default: 5)
S.encode.parallel.option = 'local'; % options: local (if multiple cores), condor, none 
S.encode.parallel.chunksize = inf;%100; % max number of datapoints (channels*time) to process on each run. Larger chunks speed up processing at the cost of memory.
% RUN FUNCTIONS
S.encode.sname=datestr(now,30); % savename
S.encode.save_residuals=1; % REQUIRES LARGE AMOUNT OF MEMORY - use "disk" memory option and adjust chunk size
S.encode.save_fitted=1; % REQUIRES LARGE AMOUNT OF MEMORY - use "disk" memory option and adjust chunk size
S.encode.save_input=1; % REQUIRES LARGE AMOUNT OF MEMORY - use "disk" memory option and adjust chunk size
S.encode.memory_option='disk';
eegstats_encode(S);

%% Create topographic images from statistics outputs
S=struct;
S.img.sname = datestr(now,30);
S.img.path.inputs = 'C:\Data\CORE\eeg\ana\stats';
S.img.file.inputs = 'stats_20191231T164634.mat';
S.img.path.outputs = strrep(fullfile(S.img.path.inputs,S.img.file.inputs),'.mat',['_' S.img.sname]); % a folder name
S.img.path.SPMdata = 'C:\Data\CORE\eeg\ana\spm\SPMdata';
S.img.file.SPMdata = 'spm12_CORE008_2_merged_cleaned.mat'; 
S.img.file.SPMimg = 'C:\Data\CORE\eeg\ana\stats\images\stats_LME_all_chan_condHGF_arcsinh_20191206T115520_20191206T192157\model_1_coeff_1_b.nii'; 
% define parameters of analysis: which models? which F contrasts?
S.img.model.index = [1]; % which models? leave blank for all
S.img.model.contrast = {}; % for each model, which contrasts? leave blank for all
S.img.model_comp.index = [0]; % any model comparisons? leave blank for all
S.img.mask_img = ''; % path/filename to a image to use as a mask
S.img.imgsize = 32;
% RUN
D=eegstats_createimages(S);

%% Multiple comparisons correction of p value images
S.MCC.path.inputs = ''; % leave blank to use paths already specified within D
S.MCC.path.outputs = ''; % leave blank to use paths already specified within D
S.MCC.path.SPMdata = 'Q:\Projects\CORE\eeg\ana\spm\SPMdata';
S.MCC.file.SPMdata = 'mspm12_flip_CORE002_2_merged_cleaned_stats_BRR_all_chan_condHGF_notrans_20190221T154622_pred1.mat'; 
S.MCC.file.SPMimg = 'Q:\Projects\CORE\eeg\ana\spm\SPMdata\sensorimages\t-200_899_b-200_0_mspm12_flip_Sideavg_CORE032_2_merged_cleaned\scondition_1.nii'; 
S.MCC.model.index = [1]; % which models? leave blank for all
S.MCC.model.contrast = {}; % for each model, which contrasts? leave blank for all
S.MCC.model_comp.index = [0]; % any model comparisons? leave blank for all
S.MCC.method = 'pTFCE';%'pTFCE';%'FDR'; % multiple comparisons correction method: FDR, pTFCE
S.MCC.thresh = 0.05;
S.MCC.smooth_est_all = 1; % estimate smoothness for every model and contrast? 0 = uses smoothness estimate from 1st model/contrast only (faster).
S.MCC.DF_type = 'LME_output'; % 'residual', 'LME_output' 
S.MCC.overwrite = 0; % overwrite previous pTFCE outputs if in same folder
S.MCC.mask_img_post = ''; % path/filename to a image to use as a mask
S.MCC.mask_img = ''; % path/filename to a image to use as a mask
D=eegstats_MCC(S,D);

%% Cluster extraction
% IDENTIFY CLUSTERS FROM MCC CORRECTED T/F IMAGES using connected
% components analysis and summarising using median or mean
S.clus.path.inputs = 'C:\Data\CORE\eeg\ana\stats\stats_20191231T164634_20200102T140010'; % leave blank to use paths already specified within D
S.clus.path.outputs = S.clus.path.inputs; % leave blank to use paths already specified within D
S.clus.model.index = [1]; % which models? leave blank for all
S.clus.model.contrast = {}; % for each model, which contrasts? leave blank for all
S.clus.model_comp.index = [0]; % any model comparisons? leave blank for all
S.clus.image = 'F';
S.clus.connected_clusters.connectivity=6;
S.clus.connected_clusters.maximise = 'prodsize'; % find threshold that maximises cluster number and/or product of cluster sizes: 'num', 'prodsize', 'num_prodsize'
S.clus.connected_clusters.ClusExtentThresh = 10; % min number of voxels within a cluster
S.clus.connected_clusters.ClusMaxNum = 10; % maximum number of clusters (only the largest clusters)
S.clus.summary_types= {'mean','median'}; % no longer using MAX: residuals tend to be biased away from zero line
S.clus.save_vols=1;
D=eegstats_extractclusters(S);

%% Create thresholded images across the range of input values
% THRESHOLD T/F IMAGES, BETA IMAGES, ETC.
S.thresh.path.inputs = 'C:\Data\CORE\eeg\ana\stats\stats_20191231T164634_20200102T140010'; % leave blank to use paths already specified within D
S.thresh.path.outputs = S.thresh.path.inputs; % leave blank to use paths already specified within D
S.thresh.model.index = [1]; % which models? leave blank for all
S.thresh.model.coeff = {}; % for each model, which contrasts? leave blank for all
S.thresh.model.contrast = {0}; % for each model, which contrasts? leave blank for all
S.thresh.model_comp.index = [0]; % any model comparisons? leave blank for all
S.thresh.image = 'b'; % Options: 'F','T','p' images from con; 'b' images from coeff; 'logl'/'waic'/'r2' from model, 'LR' from model_comp.
S.thresh.Thresh_method='range'; % options: 'range' or 'values'
S.thresh.Thresh=5; % either a single number (e.g. 5 = 5 thresholds across the full range) or a scalar/vector of actual values.
D=eegstats_threshimages(S);

%% Encoding model diagnostics from clusters
S.diag.path.inputs = ''; % leave blank to use paths already specified within D
S.diag.path.outputs = ''; % leave blank to use paths already specified within D
S.diag.model.index = [1]; % which models? leave blank for all
S.diag.model.contrast = {}; % for each model, which contrasts? leave blank for all
S.diag.model_comp.index = [0]; % any model comparisons? leave blank for all
S.diag.pred = {'ID','group','ODD','CP','CD'}; % which categorical predictors to use to assess residual variance?
S.diag.summary_types={'median'}; % {'mean','median'}; no longer using MAX: residuals tend to be biased away from zero line
D=eegstats_diagnostics(S,D);

%% Decoding with BIEM
S=struct;
S.biem.path.inputs.prep = 'C:\Data\CORE\eeg\ana\prep\cleaned\part2\eegstats_dataprep\eegstats_dtab20191230T153645.mat'; % full path and filename of prepped data, including predictors
S.biem.path.inputs.encode = 'C:\Data\CORE\eeg\ana\stats\stats_20191231T164634_20200102T140010'; % path of D.mat containing encoding model coefficients
S.biem.path.outputs = S.biem.path.inputs.encode; % leave blank to use paths already specified within D
S.biem.model.index = []; % which models? leave blank for all
S.biem.prior = 'subject_pred'; % options: 'group_pred' (use variance in whole group's predictors), 'subject_pred', 'uniform'
S.biem.randomeffects = {};%{'ID'}; % include random effects in betas?
S.biem.maskimg = '';%'C:\Data\CORE\eeg\ana\stats\stats_20191228T232523_20191229T110421\group_1_model_1_con_2_F_thresh1.nii';
S.biem.path.SPMdata = 'C:\Data\CORE\eeg\ana\spm\SPMdata'; % SPM data example file: only needed if using a masking image
S.biem.file.SPMdata = 'spm12_CORE008_2_merged_cleaned.mat'; 
S.biem.file.SPMimg = 'C:\Data\CORE\eeg\ana\stats\images\stats_LME_all_chan_condHGF_arcsinh_20191206T115520_20191206T192157\model_1_coeff_1_b.nii'; 
S.biem.plot_correlations=1;
% RUN FUNCTION
S.biem.sname=datestr(now,30); % savename
D = eegstats_decode_BIEM(S);

%% Decoding with MVPA
S=struct;
S.mvpa.path.inputs.prep = 'C:\Data\CORE\eeg\ana\prep\cleaned\part2\eegstats_dataprep\eegstats_dtab20191228T225925.mat'; % full path and filename of prepped data, including predictors
S.mvpa.path.outputs = 'C:\Data\CORE\eeg\ana\stats\'; % leave blank to use paths already specified within D
S.mvpa.path.inputs.cfglayout = 'Q:\MATLAB\projects\CORE\cosmo\modified functions\GSN92.sfp'; % layouts for cosmo
S.mvpa.path.inputs.cfgoutput = 'Q:\MATLAB\projects\CORE\cosmo\modified functions\GSN92.lay'; % layouts for cosmo
S.mvpa.maskimg = '';%'C:\Data\CORE\eeg\ana\stats\stats_20191228T232523_20191229T110421\group_1_model_1_con_2_F_thresh1.nii';
S.mvpa.path.SPMdata = 'C:\Data\CORE\eeg\ana\spm\SPMdata'; % SPM data example file: only needed if using a masking image
S.mvpa.file.SPMdata = 'spm12_CORE008_2_merged_cleaned.mat'; 
S.mvpa.file.SPMimg = 'C:\Data\CORE\eeg\ana\stats\images\stats_LME_all_chan_condHGF_arcsinh_20191206T115520_20191206T192157\model_1_coeff_1_b.nii'; 
S.mvpa.type='class'; %'class' or 'regress'
S.mvpa.pred='ODD'; % predictor variable - name of header in dtab
S.mvpa.SL_type = 'time';
S.mvpa.select_timewin = []; % empty: all time samples
S.mvpa.search_radius = 25; % data points, not ms.  Set to Inf to analyse all data points at once.
S.mvpa.use_measure = 'crossvalidation';
S.mvpa.balance_dataset_and_partitions =1; % turn on for classification, off for regression
S.mvpa.parti='take-one-out'; % 'take-one-out' (for regression), 'splithalf', 'oddeven', 'nchunks' (for classification)
S.mvpa.use_classifier = 'GP'; 
S.mvpa.use_chunks = 'none'; % 'balance_targets' (over chunks) or 'none'. Not needed for 'take-one-out' option.
S.mvpa.nchunks=0; % for classification, 10 chunks; for regression, 0.
S.mvpa.average_train_count = 1;
S.mvpa.average_train_resamplings = 1;
S.mvpa.plot_correlations=1;
% RUN FUNCTION
S.mvpa.sname=datestr(now,30); % savename
D = eegstats_decode_MVPA(S);