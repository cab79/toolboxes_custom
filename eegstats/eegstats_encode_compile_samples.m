function eegstats_encode_compile_samples(varargin) % local version
% Compiles samples from chunks and reshapes back to EEG data dimensions

% Format: eegstats_encode_compile_samples(S,C)
% C: structure of size c, containing Y (size s), chunk_info and samples (size s). s = number of samples per chunk.

%% find outputted data
if ~isempty(varargin)
    iscondor=0
    S = varargin{1};
    C = varargin{2};
    pth=S.encode.path.outputs;
    pthData=pth;
    
else
    % assume using condor
   iscondor=1
   [S,C] = load_condor_data;
    pth=S.temp.pth;
    pthData=S.temp.pthData;
end

% initiate main output, D
D=struct;

% data info
cidx = find(~cellfun(@isempty,{C(:).chunk_info}));
nD=C(cidx(1)).chunk_info.nD;
ddim = C(cidx(1)).chunk_info.dim;% EEG data dimensions for reshaping
if nD>1
    save_pref = 'subject_';
else
    save_pref = 'group_';
end

% loop through chunks of samples

for d = 1:nD % subject

    D(d).ID = S.ID{d};

    M=struct; % temporary structure
    si = [];sinz = [];
    
    %% compile samples from chunks
    
    % assume same number of chunks per subject
    n_chunks_d=C(cidx(1)).chunk_info.n_chunks_d;
    
    % chunks for all subjects - THIS VERSION WORKS FOR LMMs, BUT NOT FOR
    % BRRS
    %for c = cidx
        
    % chunks per subject - THIS VERSION WORKS FOR BRRs, UNSURE IF WORKS FOR
    % LMMs
    for nc = 1:n_chunks_d 
        
        % chunk index, c
        c = (d-1)*n_chunks_d +nc;
        
        if c>length(C); continue; end
        
        disp(['compiling output from file ' num2str(c)])
        
        si = [si C(c).chunk_info.sample_index]; % samples
        if isfield(C(c).chunk_info,'sample_index_nonzero')
            sinz = [sinz C(c).chunk_info.sample_index_nonzero]; % samples
        else
            sinz = si; % samples
        end
    
        if ~isfield(M,'model')
            M = C(c).M;
        else
            if isfield(C(c).M,'model')
                for i=1:length(M.model)
                    M.model(i).samples = [M.model(i).samples,C(c).M.model(i).samples];
                end
            end
            if isfield(C(c).M,'model_comp')
                for i=1:length(M.model_comp)
                    M.model_comp(i).samples = [M.model_comp(i).samples,C(c).M.model_comp(i).samples];
                end
            end
        end
        if d==nD; C(c).M=[]; end % save memory
    end
    

    %% reshape to EEG data dimensions
    % loop through chunks of samples
    LM=length(M.model);
    % for each model
    for i = 1:LM
        % items per model
        D(d).model(i).def = M.model(i).samples(1).def;
        D(d).model(i).fixeddesign = M.model(i).samples(1).fixeddesign;
        
        %OLD: D(d).model(i).s = reshape(vertcat([M.model(i).samples(:).mse]'),ddim(1),ddim(2));
        D(d).model(i).s = nan(ddim(1:end-1)');
        D(d).model(i).s(sinz) = vertcat([M.model(i).samples(:).mse]');
        
        % D(d).model(i).logl = reshape(vertcat([M.model(i).samples(:).logl]'),ddim(1),ddim(2));
        D(d).model(i).logl = nan(ddim(1:end-1)');
        D(d).model(i).logl(sinz) = vertcat([M.model(i).samples(:).logl]');
        
        % D(d).model(i).r2_ord = reshape(vertcat([M.model(i).samples(:).r2_ord]'),ddim(1),ddim(2));
        D(d).model(i).r2_ord = nan(ddim(1:end-1)');
        D(d).model(i).r2_ord(sinz) = vertcat([M.model(i).samples(:).r2_ord]');
        
        % items per model/coefficient
        for ii = 1:length(M.model(i).samples(1).CoefficientNames)
            D(d).model(i).coeff(ii).name = M.model(i).samples(1).CoefficientNames{ii};
            
            D(d).model(i).coeff(ii).b = nan(ddim(1:end-1)');
            D(d).model(i).coeff(ii).b(sinz) = arrayfun(@(X) X.b(ii), M.model(i).samples);
        end
        % LME models only
        if strcmp(S.encode.method,'LMM')
            D(d).model(i).pred = M.model(i).samples(1).pred;
            
            D(d).model(i).r2_adj = nan(ddim(1:end-1)');
            D(d).model(i).r2_adj(sinz) = vertcat([M.model(i).samples(:).r2_adj]');
            
            D(d).model(i).randomdesign = M.model(i).samples(1).randomdesign;
            D(d).model(i).randomnames = M.model(i).samples(1).RandomNames;
            
            if isfield(M.model(i).samples(1),'failed_s'); D(d).model(i).failed_s = M.model(i).samples(1).failed_s; end
    %         D(d).model(i).skew=reshape(vertcat([M.model(i).samples(:).skew]'),ddim(1),ddim(2));
    %         D(d).model(i).kurt=reshape(vertcat([M.model(i).samples(:).kurt]'),ddim(1),ddim(2));
    %         D(d).model(i).hnorm=reshape(vertcat([M.model(i).samples(:).hnorm]'),ddim(1),ddim(2));
            % items per model/effect
            for ii = 1:length(M.model(i).samples(1).Term)
                D(d).model(i).con(ii).term = M.model(i).samples(1).Term{ii};
                
                %D(d).model(i).con(ii).DF = reshape(vertcat(DF{:})',2,ddim(1),ddim(2));
                DF=arrayfun(@(X) X.DF(ii,:), M.model(i).samples,'UniformOutput',0);
                D(d).model(i).con(ii).DF = nan([2 ddim(1:end-1)']);
                D(d).model(i).con(ii).DF(:,sinz) = vertcat(DF{:})';
                
                D(d).model(i).con(ii).F  = nan(ddim(1:end-1)');
                D(d).model(i).con(ii).F(sinz) = arrayfun(@(X) X.F(ii), M.model(i).samples);
                
                D(d).model(i).con(ii).p = nan(ddim(1:end-1)');
                D(d).model(i).con(ii).p(sinz) = arrayfun(@(X) X.p(ii), M.model(i).samples);
            end
            % random effects
            try
                D(d).model(i).random = nan([prod(ddim(1:end-1)') height(D(d).model(i).randomnames)]);
                D(d).model(i).random(sinz,:) = vertcat([M.model(i).samples(:).r]');
                D(d).model(i).random = reshape(D(d).model(i).random,[ddim(1:end-1)' height(D(d).model(i).randomnames)]);
            catch % in case dimension are inconsistent, save as cell
                D(d).model(i).random = {M.model(i).samples(:).r};
            end
            % covariance
            D(d).model(i).coeffcov = cell(ddim(1:end-1)');
            D(d).model(i).coeffcov(sinz) = arrayfun(@(X) X.coeffcov, M.model(i).samples,'UniformOutput',0);
            for ii = 1:length(M.model(i).samples(1).psi)
                D(d).model(i).psi(ii).cov = cell(ddim(1:end-1)');
                D(d).model(i).psi(ii).cov(sinz) = arrayfun(@(X) X.psi{ii}, M.model(i).samples,'UniformOutput',0);
            end
        end
        
        % BRR models only
        if strcmp(S.encode.method,'BRR')
            D(d).model(i).waic = nan(ddim(1:end-1)');
            D(d).model(i).waic(sinz) = vertcat([M.model(i).samples(:).waic]');
            D(d).model(i).ktest = nan(ddim(1:end-1)');
            D(d).model(i).ktest(sinz) = vertcat([M.model(i).samples(:).ktest_normality]');

            % items per model/coefficient
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
        D(d).model(i).resid_file = fullfile(pthData,[save_pref num2str(d) '_resid_mi' num2str(i) '_' S.encode.sname '.mat']);
        D(d).model(i).fitted_file = fullfile(pthData,[save_pref num2str(d) '_fitted_mi' num2str(i) '_' S.encode.sname '.mat']);
        D(d).model(i).input_file = fullfile(pthData,[save_pref num2str(d) '_input_mi' num2str(i) '_' S.encode.sname '.mat']);
        % save
        if strcmp(S.encode.memory_option,'memory')
            if S.encode.save_residuals
                cd(pthData)
                
                resid = nan([prod(ddim(1:end-1)') ddim(end)]);
                resid(sinz,:) = vertcat([M.model(i).samples(:).resid]');
                resid = reshape(resid,[ddim(1:end-1)' ddim(end)]);
                
                disp([save_pref ' saving residuals from model ' num2str(i) '/' num2str(LM)]);
                save(D(d).model(i).resid_file,'resid','-v7.3');
                cd(pth)
            end

            if S.encode.save_input
                cd(pthData)
                
                input = nan([prod(ddim(1:end-1)') ddim(end)]);
                input(sinz,:) = vertcat([M.model(i).samples(:).input]');
                input = reshape(input,[ddim(1:end-1)' ddim(end)]);
                
                disp([save_pref ' saving input from model ' num2str(i) '/' num2str(LM)]);
                save(D(d).model(i).input_file,'input','-v7.3');
                cd(pth)
            end

            if S.encode.save_fitted
                cd(pthData)
                
                fitted = nan([prod(ddim(1:end-1)') ddim(end)]);
                fitted(sinz,:) = vertcat([M.model(i).samples(:).fitted]');
                fitted = reshape(fitted,[ddim(1:end-1)' ddim(end)]);
                
                disp([save_pref ' saving fitted from model ' num2str(i) '/' num2str(LM)]);
                save(D(d).model(i).fitted_file,'fitted','-v7.3');
                cd(pth)
            end
        elseif iscondor
            if S.encode.save_residuals
                resid = nan(ddim(end),length(M.model(i).samples));
                for s = 1:length(M.model(i).samples)
                    temp = M.model(i).samples(s).resid;

                    % z-score
                    resid_mean=nanmean(temp);
                    resid_std=nanstd(temp);
                    resid(:,s)=(temp-resid_mean)/resid_std;

                    % clear memory            
                    M.model(i).samples(s).resid = []; 
                end
                disp([save_pref ' saving residuals from model ' num2str(i) '/' num2str(LM)]);
                cd(pthData)
                save(D(d).model(i).resid_file,'resid','-v7.3');
                %zip('resid.zip', D(d).model(i).resid_file); 
                cd(pth)
            end
            if S.encode.save_input
                input = nan(ddim(end),length(M.model(i).samples));
                for s = 1:length(M.model(i).samples)
                    input(:,s) = M.model(i).samples(s).input;
                    % clear memory            
                    M.model(i).samples(s).input = []; 
                end
                disp([save_pref ' saving inputs from model ' num2str(i) '/' num2str(LM)]);
                cd(pthData)
                save(D(d).model(i).input_file,'input','-v7.3');
                cd(pth)
            end
            if S.encode.save_fitted
                fitted = nan(ddim(end),length(M.model(i).samples));
                for s = 1:length(M.model(i).samples)
                    fitted(:,s) = M.model(i).samples(s).fitted;
                    % clear memory            
                    M.model(i).samples(s).fitted = []; 
                end
                disp([save_pref ' saving fitted from model ' num2str(i) '/' num2str(LM)]);
                cd(pthData)
                save(D(d).model(i).fitted_file,'fitted','-v7.3');
                cd(pth)
            end
            
            
%             if S.encode.save_residuals
%                 disp([save_pref num2str(C(c).chunk_info.d) ' saving residuals from chunk ' num2str(c) ' model ' num2str(i)]);
%                 resid = vertcat([M.model(i).samples(:).resid]');
%                 m = matfile(fullfile(S.encode.path.outputs,[save_pref num2str(C(c).chunk_info.d) '_resid_mi' num2str(i) '_' S.encode.sname '.mat']),'Writable',true);
%                 varIsInMat = @(name) ~isempty(who(m, name));
%                 if varIsInMat('resid')
%                     [x,~] = size(m.resid);
%                     [~,ys] = size(m.resid_samples);
%                 else
%                     x=0;
%                     ys=0;
%                 end
%                 m.resid(x+1:x+size(resid,1),1:size(resid,2))=resid;
%                 m.resid_samples(1,ys+1:ys+size(resid,1)) = C(c).chunk_info.sample_index_nonzero;
%                 [M.model(i).samples(:).resid] = deal([]);
%               
%             end
%             if S.encode.save_fitted
%                 disp([save_pref num2str(C(c).chunk_info.d) ' saving fitted from chunk ' num2str(c) ' model ' num2str(i)]);
%                 fitted = vertcat([M.model(i).samples(:).fitted]');
%                 m = matfile(fullfile(S.encode.path.outputs,[save_pref num2str(C(c).chunk_info.d) '_fitted_mi' num2str(i) '_' S.encode.sname '.mat']),'Writable',true);
%                 varIsInMat = @(name) ~isempty(who(m, name));
%                 if varIsInMat('fitted')
%                     [x,~] = size(m.fitted);
%                     [~,ys] = size(m.fitted_samples);
%                 else
%                     x=0;
%                     ys=0;
%                 end
%                 m.fitted(x+1:x+size(fitted,1),1:size(fitted,2))=fitted;
%                 m.fitted_samples(1,ys+1:ys+size(fitted,1)) = C(c).chunk_info.sample_index_nonzero;
%                 [M.model(i).samples(:).fitted] = deal([]);
%              
%             end
%             if S.encode.save_input
%                 disp([save_pref num2str(C(c).chunk_info.d) ' saving input from chunk ' num2str(c) ' model ' num2str(i)]);
%                 input = vertcat([M.model(i).samples(:).input]');
%                 m = matfile(fullfile(S.encode.path.outputs,[save_pref num2str(C(c).chunk_info.d) '_input_mi' num2str(i) '_' S.encode.sname '.mat']),'Writable',true);
%                 varIsInMat = @(name) ~isempty(who(m, name));
%                 if varIsInMat('input')
%                     [x,~] = size(m.input);
%                     [~,ys] = size(m.input_samples);
%                 else
%                     x=0;
%                     ys=0;
%                 end
%                 m.input(x+1:x+size(input,1),1:size(input,2))=input;
%                 m.input_samples(1,ys+1:ys+size(input,1)) = C(c).chunk_info.sample_index_nonzero;
%                 [M.model(i).samples(:).input] = deal([]);
%                
%             end
        end
            
    end

    if isfield(M,'model_comp')
        MC=length(M.model_comp);
        for i = 1:MC
            D(d).model_comp(i).pair = M.model_comp(i).samples(1).pair;
            
            %D(d).model_comp(i).pval = reshape(vertcat([M.model_comp(i).samples(:).pval]'),ddim(1),ddim(2));
            D(d).model_comp(i).pval = nan(ddim(1:end-1)');
            D(d).model_comp(i).pval(sinz) = vertcat([M.model_comp(i).samples(:).pval]');
            
            %D(d).model_comp(i).LR = reshape(vertcat([M.model_comp(i).samples(:).LRStat]'),ddim(1),ddim(2));
            D(d).model_comp(i).LR = nan(ddim(1:end-1)');
            D(d).model_comp(i).LR(sinz) = vertcat([M.model_comp(i).samples(:).LRStat]');
        end
    end
end

%% save results
% varargout = {D};
%try
    disp('saving stats')
    save(fullfile(S.encode.path.outputs,[S.encode.sname '.mat']),'D','S','-v7.3');
%catch
%    error('cannot save results')
%end
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