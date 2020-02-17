% function eegstats_encode_compile_samples % condor version
function eegstats_encode_compile_samples(varargin) % local version
% Compiles samples from chunks and reshapes back to EEG data dimensions

% Format: eegstats_encode_compile_samples(S,C)
% C: structure of size c, containing Y (size s), chunk_info and samples (size s). s = number of samples per chunk.

%% find outputted data
if ~exist('varargin','var') %|| isempty(varargin)
    % assume using condor
   iscondor=1
   [S,C] = load_condor_data;
    pth=S.temp.pth;
    pthData=S.temp.pthData;
    
else
    iscondor=0
    S = varargin{1};
    C = varargin{2};
    pth=S.encode.path.outputs;
    pthData=pth;
end

% initiate main output, D
D=struct;

% data info
nD=C(1).chunk_info.nD;
ddim = C(1).chunk_info.dim;% EEG data dimensions for reshaping
if nD>1
    save_pref = 'subject_';
else
    save_pref = 'group_';
end

% loop through chunks of samples
for d = 1:nD % subject
    M=struct; % temporary structure
    %% compile samples from chunks
    % assume same number of chunks per subject
    n_chunks_d=C(1).chunk_info.n_chunks_d;
    % chunks per subject
    for nc = 1:n_chunks_d 
        % chunk index, c
        c = (d-1)*n_chunks_d +nc;
    
        if ~isfield(M,'model')
            M = C(c).M;
        else
            for i=1:length(M.model)
                M.model(i).samples = [M.model(i).samples,C(c).M.model(i).samples];
            end
            if isfield(M,'model_comp')
                for i=1:length(M.model_comp)
                    M.model_comp(i).samples = [M.model_comp(i).samples,C(c).M.model_comp(i).samples];
                end
            end
        end
        C(c).M=[]; % save memory
    end

    %% reshape to EEG data dimensions
    % loop through chunks of samples
    LM=length(M.model);
    % for each model
    for i = 1:LM
        % items per model
        D(d).model(i).def = M.model(i).samples(1).def;
        D(d).model(i).fixeddesign = M.model(i).samples(1).fixeddesign;
        D(d).model(i).s = reshape(vertcat([M.model(i).samples(:).mse]'),ddim(1),ddim(2));
        D(d).model(i).logl = reshape(vertcat([M.model(i).samples(:).logl]'),ddim(1),ddim(2));
        D(d).model(i).r2_ord = reshape(vertcat([M.model(i).samples(:).r2_ord]'),ddim(1),ddim(2));
        % items per model/coefficient
        for ii = 1:length(M.model(i).samples(1).CoefficientNames)
            D(d).model(i).coeff(ii).name = M.model(i).samples(1).CoefficientNames{ii};
            D(d).model(i).coeff(ii).b = reshape(arrayfun(@(X) X.b(ii), M.model(i).samples),ddim(1),ddim(2));
        end
        % LME models only
        if strcmp(S.encode.method,'LMM')
            D(d).model(i).pred = M.model(i).samples(1).pred;
            D(d).model(i).r2_adj = reshape(vertcat([M.model(i).samples(:).r2_adj]'),ddim(1),ddim(2));
            D(d).model(i).randomdesign = M.model(i).samples(1).randomdesign;
            D(d).model(i).randomnames = M.model(i).samples(1).RandomNames;
    %         D(d).model(i).skew=reshape(vertcat([M.model(i).samples(:).skew]'),ddim(1),ddim(2));
    %         D(d).model(i).kurt=reshape(vertcat([M.model(i).samples(:).kurt]'),ddim(1),ddim(2));
    %         D(d).model(i).hnorm=reshape(vertcat([M.model(i).samples(:).hnorm]'),ddim(1),ddim(2));
            % items per model/effect
            for ii = 1:length(M.model(i).samples(1).Term)
                D(d).model(i).con(ii).term = M.model(i).samples(1).Term{ii};
                DF=arrayfun(@(X) X.DF(ii,:), M.model(i).samples,'UniformOutput',0);
                D(d).model(i).con(ii).DF = reshape(vertcat(DF{:})',2,ddim(1),ddim(2));
                D(d).model(i).con(ii).F = reshape(arrayfun(@(X) X.F(ii), M.model(i).samples),ddim(1),ddim(2));
                D(d).model(i).con(ii).p = reshape(arrayfun(@(X) X.p(ii), M.model(i).samples),ddim(1),ddim(2));
            end
            % random effects
            D(d).model(i).random = reshape(vertcat([M.model(i).samples(:).r]'),ddim(1),ddim(2),[]);
            % covariance
            D(d).model(i).coeffcov = reshape(arrayfun(@(X) X.coeffcov, M.model(i).samples,'UniformOutput',0),ddim(1),ddim(2),[]);
            for ii = 1:length(M.model(i).samples(1).psi)
                D(d).model(i).psi(ii).cov = reshape(arrayfun(@(X) X.psi{ii}, M.model(i).samples,'UniformOutput',0),ddim(1),ddim(2),[]);
            end
        end
        
        % BRR models only
        if strcmp(S.encode.method,'BRR')
            D(d).model(i).waic = reshape(vertcat([M.model(i).samples(:).waic]'),ddim(1),ddim(2));
        end
        
        % filenames for resid, fitted, input
        D(d).model(i).resid_file = fullfile(pthData,[save_pref num2str(d) '_resid_mi' num2str(i) '_' S.encode.sname '.mat']);
        D(d).model(i).fitted_file = fullfile(pthData,[save_pref num2str(d) '_fitted_mi' num2str(i) '_' S.encode.sname '.mat']);
        D(d).model(i).input_file = fullfile(pthData,[save_pref num2str(d) '_input_mi' num2str(i) '_' S.encode.sname '.mat']);
        % save
        if strcmp(S.encode.memory_option,'memory')
            if S.encode.save_residuals
                cd(pthData)
                resid = reshape(vertcat([M.model(i).samples(:).resid]'),ddim(1),ddim(2),[]);
                disp([save_pref ' saving residuals from model ' num2str(i) '/' num2str(LM)]);
                save(D(d).model(i).resid_file,'resid','-v7.3');
                cd(pth)
            end

            if S.encode.save_input
                cd(pthData)
                input = reshape(vertcat([M.model(i).samples(:).input]'),ddim(1),ddim(2),[]);
                disp([save_pref ' saving input from model ' num2str(i) '/' num2str(LM)]);
                save(D(d).model(i).input_file,'input','-v7.3');
                cd(pth)
            end

            if S.encode.save_fitted
                cd(pthData)
                fitted = reshape(vertcat([M.model(i).samples(:).fitted]'),ddim(1),ddim(2),[]);
                disp([save_pref ' saving fitted from model ' num2str(i) '/' num2str(LM)]);
                save(D(d).model(i).fitted_file,'fitted','-v7.3');
                cd(pth)
            end
        end
    end

    if isfield(M,'model_comp')
        MC=length(M.model_comp);
        for i = 1:MC
            D(d).model_comp(i).pair = M.model_comp(i).samples(1).pair;
            D(d).model_comp(i).pval = reshape(vertcat([M.model_comp(i).samples(:).pval]'),ddim(1),ddim(2));
            D(d).model_comp(i).LR = reshape(vertcat([M.model_comp(i).samples(:).LRStat]'),ddim(1),ddim(2));
        end
    end
end

%% save results
% varargout = {D};
try
    disp('saving stats')
    save(fullfile(S.encode.path.outputs,['stats_' S.encode.sname '.mat']),'D','S');
catch
    error('cannot save results')
end
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

if ~S.encode.save_input
    try
        copyfile input0.mat example_input.mat
        delete('input*.mat');
    end
end

try
    delete('*.Z');
    delete('*.log');
    delete('*.err');
    delete('*.out');
end
try
%     delete(fullfile(pthData,'resid*.zip'));
    delete('output.zip');
    delete('input.zip');
end

nD=chunk_info.nD;
n_chunks_d=chunk_info.n_chunks_d;
index_length = nD*n_chunks_d;
C=struct;
for d = 1:nD % subject
    for nc = 1:n_chunks_d % chunks per subject
        % index
        c = (d-1)*n_chunks_d +nc;

        disp(['compiling output from file ' num2str(c) '/' num2str(index_length)])
        fileout = load(['output' num2str(c-1) '.mat']);
        C(fileout.chunk_info.c).M = fileout.M;
        C(fileout.chunk_info.c).chunk_info = fileout.chunk_info;

        if S.encode.save_input
            filein = load(['input' num2str(c-1) '.mat']);
            C(filein.chunk_info.c).Y = filein.Y;
        end
    end
end

function [in,Z] = eegstats_condor_monitor_outputs

nowtime = now;
in = dir('input*.mat');

fin=0;
while fin==0
    pause(10)
    Z = dir('output*.mat');
    if length(Z)==length(in) || length(Z)>length(in)
        complete = [Z(:).bytes]>0; %& [Z(:).datenum]>nowtime;
        if all(complete)
            fin=1;
        end
        disp(['number of outputs complete: ' num2str(sum(complete)) '/' num2str(length(complete))])
    end
end