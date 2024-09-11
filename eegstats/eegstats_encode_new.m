function eegstats_encode_new(S,varargin)
% Wrapper function; calls subfunctions for different decoding methods, or
% generates input files and job submission file for condor.
% Compiles all samples from the outputs and saves them.

%% find prepared data
if isempty(varargin)
   % try loading the data if not supplied in varargin
   try
       disp('loading data...');
       load(S.encode.path.inputs,'D') 
   catch
       error('Supply a D structure or path/file of the data')
   end
else
    D= varargin{1};
end

C=struct;
if isfield(S.encode,'nprep')
    np = S.encode.nprep;
else
    np=1;
end
cum_nc=0; % cumulative chunk number
cvalid = []; % accumulate c indices that contain data
NumWorkers = [];

% this needs creating first so all IDs exist in every Condor output file
S.ID={};
for d = 1:length(D)
    S.ID{d} = unique(D(d).prep.dtab.ID,'stable');
end

for d = 1:length(D)
    
    % NEW select whether to load the data Y from D or from a separate file
    % containing Y (if Y is too large to keep in memory entirely)
    if ~isempty(S.encode.path.inputsY)
        matY = matfile(S.encode.path.inputsY);
        [~,lenY] = size(matY,'Y');
        dim = [reshape(matY.dim,[],1); height(D(d).prep.dtab)];
    else
        Yin = D(d).prep(np).Y;
        lenY = length(Yin);
        dim = reshape(D(d).prep(np).dim,[],1);
    end
    
    %% Chunking
    % chunk the data - ideal run time is
    % about 10 mins for Condor. runs at about 250 samples per
    % minutes, per processor, so 2500 sample is
    % ideal.
    n_chunks_d = max(1,ceil(lenY/S.encode.parallel.chunksize));
    chunksize = ceil(lenY/n_chunks_d);

    %% Prepare to save residuals and other data
    if length(D)>1
        save_pref = 'subject_';
    else
        save_pref = 'group_';
    end
%     if strcmp(S.encode.memory_option,'disk')
%         LM=length(S.encode.model);
%         for i=1:LM
%             if S.encode.save_residuals
%                 disp(['creating resid file, model ' num2str(i)])
%                 m = matfile(fullfile(S.encode.path.outputs,[save_pref num2str(d) '_resid_mi' num2str(i) '_' S.encode.sname '.mat']),'Writable',true);  
%                 for nc = 1:n_chunks_d
%                     si = chunksize*(nc-1)+1 : min(chunksize*nc,lenY);
%                     if length(dim)==4
%                         m.resid(1:dim(1),1:dim(2),1:dim(3),si)=nan([reshape(dim(1:end-1),1,[]) length(si)]);
%                     elseif length(dim)==4
%                         m.resid(1:dim(1),1:dim(2),si)=nan([reshape(dim(1:end-1),1,[]) length(si)]);
%                     end
%                 end
%                 %save(fullfile(S.encode.path.outputs,[save_pref num2str(d) '_resid_mi' num2str(i) '_' S.encode.sname '.mat']),'resid','-v7.3');
%             end
%             % fitted
%             if S.encode.save_fitted
%                 disp(['creating fitted file, model ' num2str(i)])
%                 m = matfile(fullfile(S.encode.path.outputs,[save_pref num2str(d) '_fitted_mi' num2str(i) '_' S.encode.sname '.mat']),'fitted','-v7.3');
%                 for ti = dim(end):-1:1
%                     if length(dim)==4
%                         m.fitted(1:dim(1),1:dim(2),1:dim(3),ti)=nan(reshape(dim(1:end-1),1,[]));
%                     elseif length(dim)==4
%                         m.fitted(1:dim(1),1:dim(2),ti)=nan(reshape(dim(1:end-1),1,[]));
%                     end
%                 end
%             end
%             % input
%             if S.encode.save_input
%                 disp(['creating input file, model ' num2str(i)])
%                 m = matfile(fullfile(S.encode.path.outputs,[save_pref num2str(d) '_input_mi' num2str(i) '_' S.encode.sname '.mat']),'input','-v7.3');
%                 for ti = dim(end):-1:1
%                     if length(dim)==4
%                         m.input(1:dim(1),1:dim(2),1:dim(3),ti)=nan(reshape(dim(1:end-1),1,[]));
%                     elseif length(dim)==4
%                         m.input(1:dim(1),1:dim(2),ti)=nan(reshape(dim(1:end-1),1,[]));
%                     end
%                 end
%             end
%         end
%     end
    
    for nc = 1:n_chunks_d
        si = chunksize*(nc-1)+1 : min(chunksize*nc,lenY);

        % chunk info
        chunk_info.S=lenY;
        chunk_info.nD=length(D);
        chunk_info.n_chunks_d = n_chunks_d;
        chunk_info.sample_index=si;
        chunk_info.d=d;
        chunk_info.dim=dim;

        % chunk index
        c = (d-1)*n_chunks_d +nc;
        index_length = chunk_info.nD*n_chunks_d;
        chunk_info.c=c;
        
        % setup
        if isempty(NumWorkers)
            NumWorkers = set_workers(S);
        end
        
        Y=struct;
        
        % combine design matrix with data
        if ~isempty(S.encode.path.inputsY)
            disp('loading samples from matfile...');
            tempY = matY.Y(:,si); 
        else
            Yy = Yin(si);
            tempY=zeros(height(Yy(1).dtab),length(si));
            for s = 1:length(si)
                tempY(:,s) = Yy(s).dtab.data;
            end
        end
            
        % find zero values
        iszero = find(all(tempY==0 | isnan(tempY),1));
        nonzero = si;
        nonzero(iszero)=[];
        tempY(:,iszero) = [];
        chunk_info.sample_index_nonzero=nonzero;

        if ~isfield(S.encode,'variables')
            S.encode.variables = D(d).prep(np).dtab.Properties.VariableNames;
        end
        if ~strcmp(S.encode.parallel.option,'condor')
            % otherwise, condor input files should not contain tables,
            % so reduce file sizes
            for s = 1:size(tempY,2)
                temp = array2table(tempY(:,s),'VariableNames',{'data'});
                Y(s).dtab = horzcat(D(d).prep(np).dtab,temp);
                Y(s).dtab = Y(s).dtab(:,ismember(Y(s).dtab.Properties.VariableNames,vertcat(S.encode.variables',{'data'})));
            end
        else
            dtab = D(d).prep(np).dtab(:,ismember(D(d).prep(np).dtab.Properties.VariableNames,S.encode.variables));
        end
        
        if ~isempty(tempY)
            cvalid = [cvalid c];
        else 
            disp(['chunk ' num2str(c) ' contains no nonzero samples - moving on to the next one'])
            continue
        end
        
        % save
        if strcmp(S.encode.parallel.option,'condor')
            disp(['creating input file ' num2str(length(cvalid)) '/' num2str(index_length)])
            save(fullfile(S.encode.path.outputs,['input' num2str(length(cvalid)-1) '.mat']),'tempY','dtab','S','chunk_info');
            continue
        else
            disp(['creating chunk ' num2str(c) '/' num2str(index_length)])
            C(c).Y=Y;
            C(c).M=struct;
            C(c).chunk_info=chunk_info;
        end
        
        model_run = 0;
        %if ~isempty(S.encode.path.inputsY)
            cum_nc=cum_nc+1;
            % run on a certain number of chunks and clear memory
            if cum_nc==NumWorkers
                model_run = 1;
                disp(['running model on ' num2str(NumWorkers) 'workers'])
                Csub(1:NumWorkers) = C(cvalid(end-NumWorkers+1:end));
                C(cvalid(end-NumWorkers+1:end)) = run_model(Csub,S,NumWorkers,save_pref);
                cum_nc=0;
            end
        %end
        
    end
    
end

% run models
if strcmp(S.encode.parallel.option,'condor')
    setup_condor(S,length(cvalid));
    
elseif strcmp(S.encode.parallel.option,'local') || strcmp(S.encode.parallel.option,'none')
    
    % run all samples if not run within previous loop
    if model_run==0
        %C = run_model(C,S,min(length(C),NumWorkers),save_pref);
        indices = find(~cellfun(@isempty, {C.Y}));
        Csub = C(indices);
        C(indices) = run_model(Csub,S,length(indices),save_pref);
    end
    
    save('temp.mat','S','C') % save in case this next function goes wrong!
    eegstats_encode_compile_samples(S,C);
    delete('temp.mat') % tidy up
end

%% Prepare Condor job submission, or setup parallel processing
function setup_condor(S,index_length)

methodfunc=['eegstats_encode_' S.encode.method];

% create_job_submission_file
pth=S.encode.path.outputs;
nf=index_length;
disp('creating job submission file')
A = {
    ['executable=' methodfunc '.exe']
    'indexed_input_files=input.mat'
    'indexed_output_files=output.mat'
    ['indexed_stdout=' methodfunc '.out']
    ['indexed_stderr=' methodfunc '.err']
    ['indexed_log=' methodfunc '.log']
    'runtime=240'
    ['total_jobs=' num2str(nf)]
};

fid = fopen(fullfile(pth, [methodfunc '_run.sub']),'w');
for i = 1:size(A,1)
    fprintf(fid,'%s\n',A{i});
end
fclose(fid);
quit

function NumWorkers = set_workers(S)

if strcmp(S.encode.parallel.option,'condor')
    NumWorkers = 0;
elseif strcmp(S.encode.parallel.option,'local')
    
    % Set the desired number of workers
    if ~isfield(S.encode.parallel,'Nworkers')
        clusterInfo = parcluster;
        maxNumWorkers = clusterInfo.NumWorkers;
        desiredNumWorkers = maxNumWorkers; 
    else
        desiredNumWorkers = S.encode.parallel.Nworkers; 
    end
    
    % Parallel processing
    checkp = gcp('nocreate');  % Check if a parallel pool already exists
    if isempty(checkp)
        myPool = parpool(desiredNumWorkers);  % Create a new pool with the specified number of workers
        NumWorkers = myPool.NumWorkers;
    else
        NumWorkers = checkp.NumWorkers;  % Use the existing pool
        if NumWorkers ~= desiredNumWorkers
            % If the existing pool has a different number of workers, delete it and create a new one
            delete(gcp);
            myPool = parpool(desiredNumWorkers);  % Create a new pool with the specified number of workers
            NumWorkers = myPool.NumWorkers;
        end
    end

else
    % process on single core
    NumWorkers = 1;
end


%% run analysis locally
function C=run_model(C,S,NumWorkers,save_pref)
methodfunc=['eegstats_encode_' S.encode.method];
func_model=str2func(methodfunc);
LC=NumWorkers;

if strcmp(S.encode.parallel.option,'local')
    
    % parallel processing
    nanpad=nan(1,LC);
    nanpad(1:LC)=1:LC;
    new_ind = reshape(nanpad,LC,[])';
    for ni = 1:size(new_ind,1)
        cc=new_ind(ni,:);
        cc=cc(~isnan(cc));
        parfor (c = cc)
            disp(['running model: chunk ' num2str(c) '/' num2str(LC)])
            C(c).M=func_model(S,C(c).Y);
            C(c).Y=[]; % save memory
        end
        if strcmp(S.encode.memory_option,'disk')
            for c=cc
                C=save_matfiles(S,C,c,save_pref);
            end
        end
    end
else
    % process on single core
    for c = 1:LC
        disp(['running model: chunk ' num2str(c) '/' num2str(LC)])
        C(c).M=func_model(S,C(c).Y);
        C(c).Y=[]; % save memory
        if strcmp(S.encode.memory_option,'disk')
            C=save_matfiles(S,C,c,save_pref);
        end

    end
end

function C=save_matfiles(S,C,c,save_pref)
% DON'T REALLY NEED THIS CODE AS RESIDUALS ARE ALREADY IN THE M struct? but
% memory may be too large
if S.encode.save_residuals
    for i = 1:length(C(c).M.model)
        disp([save_pref num2str(C(c).chunk_info.d) ' saving residuals from chunk ' num2str(c) ' model ' num2str(i)]);
        resid = vertcat([C(c).M.model(i).samples(:).resid]');
        m = matfile(fullfile(S.encode.path.outputs,[save_pref num2str(C(c).chunk_info.d) '_resid_mi' num2str(i) '_' S.encode.sname '.mat']),'Writable',true);
        varIsInMat = @(name) ~isempty(who(m, name));
        if varIsInMat('resid')
            [x,~] = size(m.resid);
            [~,ys] = size(m.resid_samples);
        else
            x=0;
            ys=0;
        end
        m.resid(x+1:x+size(resid,1),1:size(resid,2))=resid;
        m.resid_samples(1,ys+1:ys+size(resid,1)) = C(c).chunk_info.sample_index_nonzero;
        [C(c).M.model(i).samples(:).resid] = deal([]);
    end
end
if S.encode.save_fitted
    for i = 1:length(C(c).M.model)
        disp([save_pref num2str(C(c).chunk_info.d) ' saving fitted from chunk ' num2str(c) ' model ' num2str(i)]);
        fitted = vertcat([C(c).M.model(i).samples(:).fitted]');
        m = matfile(fullfile(S.encode.path.outputs,[save_pref num2str(C(c).chunk_info.d) '_fitted_mi' num2str(i) '_' S.encode.sname '.mat']),'Writable',true);
        varIsInMat = @(name) ~isempty(who(m, name));
        if varIsInMat('fitted')
            [x,~] = size(m.fitted);
            [~,ys] = size(m.fitted_samples);
        else
            x=0;
            ys=0;
        end
        m.fitted(x+1:x+size(fitted,1),1:size(fitted,2))=fitted;
        m.fitted_samples(1,ys+1:ys+size(fitted,1)) = C(c).chunk_info.sample_index_nonzero;
        [C(c).M.model(i).samples(:).fitted] = deal([]);
    end
end
if S.encode.save_input
    for i = 1:length(C(c).M.model)
        disp([save_pref num2str(C(c).chunk_info.d) ' saving input from chunk ' num2str(c) ' model ' num2str(i)]);
        input = vertcat([C(c).M.model(i).samples(:).input]');
        m = matfile(fullfile(S.encode.path.outputs,[save_pref num2str(C(c).chunk_info.d) '_input_mi' num2str(i) '_' S.encode.sname '.mat']),'Writable',true);
        varIsInMat = @(name) ~isempty(who(m, name));
        if varIsInMat('input')
            [x,~] = size(m.input);
            [~,ys] = size(m.input_samples);
        else
            x=0;
            ys=0;
        end
        m.input(x+1:x+size(input,1),1:size(input,2))=input;
        m.input_samples(1,ys+1:ys+size(input,1)) = C(c).chunk_info.sample_index_nonzero;
        [C(c).M.model(i).samples(:).input] = deal([]);
    end
end

% function C=save_matfiles(S,C,c,save_pref)
% [x,y,z] = ind2sub(C(c).chunk_info.dim(1:end-1),C(c).chunk_info.sample_index_nonzero(1):C(c).chunk_info.sample_index_nonzero(end));
% if S.encode.save_residuals
%     for i = 1:length(C(c).M.model)
%         disp([save_pref num2str(C(c).chunk_info.d) ' saving residuals from chunk ' num2str(c) ' model ' num2str(i)]);
%         resid_temp = vertcat([C(c).M.model(i).samples(:).resid]');
%         resid = nan(length(C(c).chunk_info.sample_index_nonzero(1):C(c).chunk_info.sample_index_nonzero(end)), size(resid_temp,2));
%         resid(C(c).chunk_info.sample_index_nonzero-C(c).chunk_info.sample_index_nonzero(1)+1,:) = resid_temp;
%         m = matfile(fullfile(S.encode.path.outputs,[save_pref num2str(C(c).chunk_info.d) '_resid_mi' num2str(i) '_' S.encode.sname '.mat']),'Writable',true);
%         %resid=permute(resid,[1 3 2]);
%         for ci = 1:length(resid)
%             try
%                 if length(C(c).chunk_info.dim)==4
%                     m.resid(x(ci),y(ci),z(ci),:)=resid(ci,:);
%                 elseif length(C(c).chunk_info.dim)==3
%                     m.resid(x(ci),y(ci),:)=resid(ci,:);
%                 end
%             catch
%                 resid;
%             end
%         end
%         [C(c).M.model(i).samples(:).resid] = deal([]);
%     end
% end
% if S.encode.save_fitted
%     for i = 1:length(C(c).M.model)
%         disp([save_pref num2str(C(c).chunk_info.d) ' saving fitted from chunk ' num2str(c) ' model ' num2str(i)]);
%         fitted = vertcat([C(c).M.model(i).samples(:).fitted]');
%         m = matfile(fullfile(S.encode.path.outputs,[save_pref num2str(C(c).chunk_info.d) '_fitted_mi' num2str(i) '_' S.encode.sname '.mat']),'Writable',true);
%         fitted=permute(fitted,[1 3 2]);
%         for ci = unique(col)
%             m.fitted(row(col==ci),unique(col(col==ci)),:)=fitted(col==ci,:,:);
%         end
%         [C(c).M.model(i).samples(:).fitted] = deal([]);
%     end
% end
% if S.encode.save_input
%     for i = 1:length(C(c).M.model)
%         disp([save_pref num2str(C(c).chunk_info.d) ' saving input from chunk ' num2str(c) ' model ' num2str(i)]);
%         input = vertcat([C(c).M.model(i).samples(:).input]');
%         m = matfile(fullfile(S.encode.path.outputs,[save_pref num2str(C(c).chunk_info.d) '_input_mi' num2str(i) '_' S.encode.sname '.mat']),'Writable',true);
%         input=permute(input,[1 3 2]);
%         for ci = unique(col)
%             m.input(row(col==ci),unique(col(col==ci)),:)=input(col==ci,:,:);
%         end
%         [C(c).M.model(i).samples(:).input] = deal([]);
%     end
% end
