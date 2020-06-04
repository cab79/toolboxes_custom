function eegstats_encode(S,varargin)
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
%% Chunking
% chunk the data - ideal run time is
% about 10 mins for Condor. runs at about 250 samples per
% minutes, per processor, so 2500 sample is
% ideal.
C=struct;
if isfield(S.encode,'nprep')
    np = S.encode.nprep;
else
    np=1;
end
for d = 1:length(D)
    Yin = D(d).prep(np).Y;
    n_chunks_d = max(1,ceil(length(Yin)/S.encode.parallel.chunksize));
    chunksize = ceil(length(Yin)/n_chunks_d);

    for nc = 1:n_chunks_d
        si = chunksize*(nc-1)+1 : min(chunksize*nc,length(Yin));

        % chunk info
        chunk_info.S=length(Yin);
        chunk_info.nD=length(D);
        chunk_info.n_chunks_d = n_chunks_d;
        chunk_info.sample_index=si;
        chunk_info.d=d;
        chunk_info.dim=reshape(D(d).prep(np).dim,[],1);

        % chunk index
        c = (d-1)*n_chunks_d +nc;
        index_length = chunk_info.nD*n_chunks_d;
        chunk_info.c=c;
        
        % combine design matrix with data
        Yy = Yin(si);
        Y=struct;
        for s = 1:length(Yy)
            Y(s).dtab = horzcat(D(d).prep(np).dtab,Yy(s).dtab);
            Y(s).data_mean = Yy(s).data_mean;
            Y(s).data_std = Yy(s).data_std;
        end

        % save
        if strcmp(S.encode.parallel.option,'condor')
            disp(['creating input file ' num2str(c) '/' num2str(index_length)])
            save(fullfile(S.encode.path.outputs,['input' num2str(c-1) '.mat']),'Y','S','chunk_info');
        else
            disp(['creating chunk ' num2str(c) '/' num2str(index_length)])
            C(c).Y=Y;
            C(c).chunk_info=chunk_info;
        end
    end
    %% Prepare to save residuals and other data
    if strcmp(S.encode.memory_option,'disk')
        if length(D)>1
            save_pref = 'subject_';
        else
            save_pref = 'group_';
        end
        LM=length(S.encode.model);
        for i=1:LM
            if S.encode.save_residuals
                disp(['creating resid file, model ' num2str(i)])
                resid=nan(reshape(D(d).prep(np).dim,1,[]));
                save(fullfile(S.encode.path.outputs,[save_pref num2str(d) '_resid_mi' num2str(i) '_' S.encode.sname '.mat']),'resid','-v7.3');
            end
            % fitted
            if S.encode.save_fitted
                disp(['creating fitted file, model ' num2str(i)])
                fitted=nan(reshape(D(d).prep(np).dim,1,[]));
                save(fullfile(S.encode.path.outputs,[save_pref num2str(d) '_fitted_mi' num2str(i) '_' S.encode.sname '.mat']),'fitted','-v7.3');
            end
            % input
            if S.encode.save_input
                disp(['creating input file, model ' num2str(i)])
                input=nan(reshape(D(d).prep(np).dim,1,[]));
                save(fullfile(S.encode.path.outputs,[save_pref num2str(d) '_input_mi' num2str(i) '_' S.encode.sname '.mat']),'input','-v7.3');
            end
        end
    end
end

%% Prepare Condor job submission, or run analysis locally
methodfunc=['eegstats_encode_' S.encode.method];
LC=length(C);
func_model=str2func(methodfunc);
if strcmp(S.encode.parallel.option,'condor')
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
elseif strcmp(S.encode.parallel.option,'local')
    % parallel processing
    checkp = gcp('nocreate');
    if isempty(checkp)
        myPool = parpool;
        NumWorkers = myPool.NumWorkers;
    else
        NumWorkers = checkp.NumWorkers;
    end

    nanpad=nan(1,ceil(LC/NumWorkers)*NumWorkers);
    nanpad(1:LC)=1:LC;
    new_ind = reshape(nanpad,NumWorkers,[])';
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
if strcmp(S.encode.parallel.option,'local') || strcmp(S.encode.parallel.option,'none')
    save('temp.mat','S','C') % save in case this next function goes wrong!
    eegstats_encode_compile_samples(S,C);
    delete('temp.mat') % tidy up
end

function C=save_matfiles(S,C,c,save_pref)
[row,col] = ind2sub(C(c).chunk_info.dim(1:2),C(c).chunk_info.sample_index);
if S.encode.save_residuals
    for i = 1:length(C(c).M.model)
        disp([save_pref num2str(C(c).chunk_info.d) ' saving residuals from chunk ' num2str(c) ' model ' num2str(i)]);
        resid = vertcat([C(c).M.model(i).samples(:).resid]');
        m = matfile(fullfile(S.encode.path.outputs,[save_pref num2str(C(c).chunk_info.d) '_resid_mi' num2str(i) '_' S.encode.sname '.mat']),'Writable',true);
        resid=permute(resid,[1 3 2]);
        for ci = unique(col)
            m.resid(row(col==ci),unique(col(col==ci)),:)=resid(col==ci,:,:);
        end
        [C(c).M.model(i).samples(:).resid] = deal([]);
    end
end
if S.encode.save_fitted
    for i = 1:length(C(c).M.model)
        disp([save_pref num2str(C(c).chunk_info.d) ' saving fitted from chunk ' num2str(c) ' model ' num2str(i)]);
        fitted = vertcat([C(c).M.model(i).samples(:).fitted]');
        m = matfile(fullfile(S.encode.path.outputs,[save_pref num2str(C(c).chunk_info.d) '_fitted_mi' num2str(i) '_' S.encode.sname '.mat']),'Writable',true);
        fitted=permute(fitted,[1 3 2]);
        for ci = unique(col)
            m.fitted(row(col==ci),unique(col(col==ci)),:)=fitted(col==ci,:,:);
        end
        [C(c).M.model(i).samples(:).fitted] = deal([]);
    end
end
if S.encode.save_input
    for i = 1:length(C(c).M.model)
        disp([save_pref num2str(C(c).chunk_info.d) ' saving input from chunk ' num2str(c) ' model ' num2str(i)]);
        input = vertcat([C(c).M.model(i).samples(:).input]');
        m = matfile(fullfile(S.encode.path.outputs,[save_pref num2str(C(c).chunk_info.d) '_input_mi' num2str(i) '_' S.encode.sname '.mat']),'Writable',true);
        input=permute(input,[1 3 2]);
        for ci = unique(col)
            m.input(row(col==ci),unique(col(col==ci)),:)=input(col==ci,:,:);
        end
        [C(c).M.model(i).samples(:).input] = deal([]);
    end
end
