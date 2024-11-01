function D=eegstats_MCCv3(S,varargin)
dbstop if error

%% find D saved from createimages function
if isempty(varargin)
   % try loading the data if not supplied in varargin
   try
       load(fullfile(S.MCC.path.inputs,'D.mat'),'D') 
   catch
       error('Supply a D structure or path/file of the data')
   end
else
    D= varargin{1};
end

% set paths
S.path.code = {
    1, 'Q:\MATLAB\toolboxes_external\spm12' % SPM
    0, 'Q:\MATLAB\toolboxes_custom\various' % FDR
    0, 'Q:\MATLAB\toolboxes_external\spm12\toolbox\pTFCE' % pTFCE
    };
set_paths(S)

if length(D)>1
    save_pref = 'subject_';
else
    save_pref = 'group_';
end

for d=1:length(D)
    
    % get model and contrast indices
    if isempty(S.MCC.model.index)
        S.MCC.model.index = 1:length(D(d).model);
    end
    if isempty(S.MCC.model.contrast) && isfield(D(d).model(1),'con')
        for i = S.MCC.model.index
            S.MCC.model.contrast{i} = 1:length(D(d).model(i).con);
        end
    end
    if isempty(S.MCC.model_comp.index)
        try
            S.MCC.model_comp.index = 1:length(D(d).model_comp);
        catch
            S.MCC.model_comp.index = 0;
        end
    end
    
    % initialise V for saving images later
    if ~isempty(S.img.file.coord)
        V=makeV;
    else
        V=makeV3; % assume MNI volumn
    end

    % mask
    if ~isempty(S.MCC.mask_img)
        mask_img = double(spm_read_vols(spm_vol(S.MCC.mask_img))>0);
        mask_img(mask_img==0) = nan;
    else
        try
            mii=find(~cellfun(@isempty,{D(d).model(:).mask_img_file}));
            mask_img = spm_read_vols(spm_vol(D(d).model(mii(1)).mask_img_file));
        catch
            mii=find(~cellfun(@isempty,{D(d).model_comp(:).mask_img_file}));
            mask_img = spm_read_vols(spm_vol(D(d).model_comp(mii(1)).mask_img_file));
        end
    end
    V.dim(1:length(size(mask_img))) = size(mask_img);
    
    if S.MCC.model.index
        for i = S.MCC.model.index

            %- Filename of mapped mask image
            if ~isempty(S.MCC.mask_img)
                VM    = S.MCC.mask_img;
            else
                VM    = D(d).model(i).mask_img_file;
            end
                
            if S.MCC.reestimate || (S.MCC.estimate && (S.MCC.use_pTFCE || strcmp(S.MCC.method,'FWE_RF')) && (~isfield(D(d).model(i),'smooth') || isempty(D(d).model(i).smooth)))
                D(d).model(i).resid_vol_file = strrep(D(d).model(i).resid_file,'.mat','_vol.mat');
                D(d).model(i).resid_vol_dir = strrep(D(d).model(i).resid_vol_file,'.mat','_dir');

                if ~exist(D(d).model(i).resid_vol_dir,'dir') || length(dir(fullfile(D(d).model(i).resid_vol_dir,'*.nii'))) ~= length(D.model(i).fixeddesign)
                    % mkdir(D(d).model(i).resid_vol_dir)
                    
                    vr=[];
                    if exist(D(d).model(i).resid_vol_file,'file')
                        disp(['MCC: for model ' num2str(i) ', loading resid vol file'])
                        load(D(d).model(i).resid_vol_file);
                        
                    elseif exist(D(d).model(i).resid_file,'file')
                        disp(['MCC: for model ' num2str(i) ', loading resid file'])
                        
                        
                        if S.MCC.load_resid_as_matfile
                            vrmat = matfile(D(d).model(i).resid_file);
                        else
                            vr = load(D(d).model(i).resid_file);
                            % % TEMPORARY CODE - REMOVE
                            % try
                            %     vr.resid=vr.fitted;
                            %     vr = rmfield(vr,'fitted');
                            % end
                        end
                        
                        % if a 3D resid file of channels x timepoints x samples
                        if ~exist('resid_vol','var') && ~isempty(S.img.file.coord) && ~isempty(vr) && ndims(vr.resid)>2
                            % % chunk it to reduce memory load
                            % if S.MCC.save_vols_chunksize
                            %     chunksize = S.clus.save_vols_chunksize;
                            % else 
                            %     chunksize = size(vr.resid,3);
                            % end
                            % n_chunks = max(1,ceil(size(vr.resid,3)/chunksize));
                            % resid_vol=nan(S.img.imgsize,S.img.imgsize,size(vr.resid,2),size(vr.resid,3));
                            % D(d).model(i).resid_vol_file = strrep(D(d).model(i).resid_file,'.mat','_vol.mat');
                            % save(D(d).model(i).resid_vol_file,'resid_vol','-v7.3');
                            % clear resid_vol;
                            % m = matfile(D(d).model(i).resid_vol_file,'Writable',true);
                            % for nc = 1:n_chunks
                            %     si = chunksize*(nc-1)+1 : min(chunksize*nc,size(vr.resid,3));
                            %     disp(['MCC: for model ' num2str(i) ', creating resid vol, chunk ' (num2str(nc))])
                            %     m.resid_vol(1:S.img.imgsize,1:S.img.imgsize,1:size(vr.resid,2),si)=topotime_3D(vr.resid(:,:,si),S);
                            % end
                            % load(D(d).model(i).resid_vol_file);
                            vr.resid = reshape(vr.resid, [], size(vr.resid,3));
                        elseif ~isempty(vr) && ndims(vr.resid)>2
                            D(d).model(i).resid_vol_file = D(d).model(i).resid_file;
                            if ~exist('resid_vol','var')
                                resid_vol=vr.resid;
                            end
                        else
                            disp('using 2D residuals')
                        end
                        %delete(D(d).model(i).resid_file)
                    elseif isfield(D(d).model(i),'resid') && ~isempty(D(d).model(i).resid)
                        disp('using residuals within data structure')
                    end

                    
                    disp('MCC: standardise residuals') % chunked to reduce memory load
                    if exist('resid_vol','var')
                        ndim_resid = ndims(resid_vol);
                        for nc = 1:size(resid_vol,3)
                            resid_vol_mean(:,:,nc)=mean(resid_vol(:,:,nc,:),ndim_resid);
                            resid_vol_std(:,:,nc)=std(resid_vol(:,:,nc,:),[],ndim_resid);
                            resid_vol(:,:,nc,:)= bsxfun(@rdivide,bsxfun(@minus,resid_vol(:,:,nc,:),resid_vol_mean(:,:,nc)),resid_vol_std(:,:,nc));
                        end
                        szR = size(resid_vol,ndim_resid);

                    elseif (isfield(D(d).model(i),'resid') && ~isempty(D(d).model(i).resid)) || (~isempty(vr) || exist('vrmat','var'))
                        
                        szR = size(D(d).model(i).fixeddesign,1);
                        dim_resid = [size(D(d).model(i).s), szR];
                        ndim_resid = length(dim_resid);
                        
                        if isfield(D(d).model(i),'resid') && ~isempty(D(d).model(i).resid) % if resids are part of D
                            maxs = length(D(d).model(i).resid);
                            lasts=0;
                            for s = 1:maxs

                                % monitor progress
                                prog=floor(s*100/maxs);
                                if prog>lasts
                                    disp(['standardising: ' num2str(prog) '%']);
                                    lasts = prog;
                                end

                                resid_mean=nanmean(D(d).model(i).resid(s).samples);
                                resid_std=nanstd(D(d).model(i).resid(s).samples);
                                D(d).model(i).resid(s).samples=(D(d).model(i).resid(s).samples-resid_mean)/resid_std;
                            end
                        elseif ~isempty(vr) && ndims(vr.resid)==2
                            % separate resid file
                            % check dimensions
                            obs_dim = find(size(vr.resid)==szR);
                            if obs_dim==1
                                maxs = size(vr.resid,2);

                            elseif obs_dim==2
                                maxs = size(vr.resid,1);

                            end
                            

                            % check if z-scoring needed - likely already z-scored on
                            % condor (in eegstats_encode_compile_samples.m)
                            check_mean = mean(vr.resid,obs_dim);
                            check_std = std(vr.resid,[],obs_dim);

                            % Tolerances (you can adjust these based on how close you want to be to zero/one)
                            meanTolerance = 0.1;  % tolerance for how close to zero the means should be
                            stdTolerance = 0.1;   % tolerance for how close to one the stds should be
                            
                            % Check if all means are approximately zero
                            meansCheck = all(abs(check_mean) < meanTolerance);
                            
                            % Check if all standard deviations are approximately one
                            stdsCheck = all(abs(check_std - 1) < stdTolerance);
                            
                            % Output results
                            if meansCheck && stdsCheck
                                disp('The data are approximately standardised.');
                            else
                                disp('Standardising...');
                                vr.resid = bsxfun(@rdivide, bsxfun(@minus, vr.resid, check_mean), check_std);
                            end

                            
                        elseif S.MCC.load_resid_as_matfile
                            % separate resid file - already z-scored on
                            % condor (in eegstats_encode_compile_samples.m)
                            maxs = size(vrmat,'resid');
                            maxs = maxs(2);
                            
                        end
                        
                        % empty volume 
                        resid_empty = nan(dim_resid(1:end-1));
                    end

                    %% Smoothing

                    % Initialize variables and structures
                    pcc = 0;  % Progress counter
                    headsize = 200; % mm
                    timeres = 1;
                    desired_fwhm = S.MCC.desired_fwhm;  % Desired smoothing in mm, mm, ms
                    use_resid_vol = exist('resid_vol', 'var');  % Flag for using resid_vol
                    
                    % Determine voxel size for smoothing
                    if use_resid_vol
                        szhead = size(mask_img, 1:3);  % Assume resid_vol has consistent dimensions
                    else
                        szhead = size(mask_img);
                    end
                    voxel_size = [headsize / szhead(1), headsize / szhead(2), timeres]; % in mm, mm, ms 
                    fwhm_voxels = desired_fwhm ./ voxel_size;
                    
                    % Initialize 4D array for storing residuals if using method 2
                    if ~use_resid_vol
                        resid_vol = zeros([szhead, szR]);
                    end
                    
                    % Precompute sample indices once for all methods
                    if ~isempty(vr) || S.MCC.load_resid_as_matfile
                        sidx = find(~isnan(D(d).model(i).s(:)))'; % sample indices within image
                        if length(sidx) ~= maxs % check
                            error('wrong number of samples in resid data');
                        end
                    end
                    
                    
                    
                    % Loop through each sample (no need for parallel processing if only storing in 4D array)
                    
                    for tr = 1:szR
                        % Monitor progress
                        pc = floor(tr * 100 / szR);
                        if pcc < pc
                            pcc = pc;
                            disp(['processing resid volumes: ' num2str(pc) '%, time: ' num2str(toc) 's']);
                            tic
                        end
                    
                        % Initialize residual image for current iteration
                        resid = resid_empty;
                    
                        % Method 1: Use precomputed residual volumes (resid_vol)
                        if use_resid_vol
                            if ndim_resid == 4
                                resid = resid_vol(:, :, :, tr);  % Use data directly from resid_vol
                            elseif ndim_resid == 3
                                resid = resid_vol(:, :, tr);  % Use data directly from resid_vol
                            end
                    
                            % Apply smoothing
                            if any(fwhm_voxels)
                                resid_smooth = zeros(size(resid));
                                if isfield(S.MCC,'fwhm_freq_pnts') && ~isempty(S.MCC.fwhm_freq_pnts)
                                    % smooth within frequency bands separately
                                    cumPnts = cumsum(S.MCC.fwhm_freq_pnts);
                                    for f = 1:length(S.MCC.fwhm_freq_pnts)
                                        if f == 1
                                            startIndex = 1;
                                        else
                                            startIndex = cumPnts(f-1) + 1;
                                        end
                                        endIndex = cumPnts(f);
                                        indices = startIndex:endIndex;
                                        resid_temp = resid(:,:,indices,:);
                                        resid_smooth_temp = zeros(size(resid_temp));
                                        spm_smooth(resid_temp, resid_smooth_temp, fwhm_voxels);
                                        resid_smooth(:,:,indices,:) = resid_smooth_temp;
                                    end
                                else
                                    spm_smooth(resid, resid_smooth, fwhm_voxels);
                                end
                                resid = resid_smooth .* double(mask_img); % Mask it
                                resid(isnan(resid)) = 0; % No NaN for smoothness function
                            end
                    
                            % Update resid_vol with smoothed data
                            resid_vol(:, :, :, tr) = resid;  % Store the smoothed volume
                    
                        % Method 2: Directly use residuals and store in 4D array
                        elseif ~isempty(vr) && ndims(vr.resid) == 2
                            if obs_dim==1
                                resid(sidx) = vr.resid(tr, :); % Load residuals into image
                            elseif obs_dim==2
                                resid(sidx) = vr.resid(:, tr)'; % Load residuals into image
                            end

                            % resid_vol = reshape(vr.resid',[size(resid) size(vr.resid,1)]);
                            % tic
                            % resid_vol = topotime_3D(resid_vol, S);
                            % toc
                    
                            % If topotime data
                            if ndims(resid) == 2
                                resid = topotime_3D(resid, S);
                            end
                    
                            % Apply smoothing
                            if any(fwhm_voxels)
                                resid_smooth = zeros(size(resid));
                                if isfield(S.MCC,'fwhm_freq_pnts') && ~isempty(S.MCC.fwhm_freq_pnts)
                                    % smooth within frequency bands separately
                                    cumPnts = cumsum(S.MCC.fwhm_freq_pnts);
                                    for f = 1:length(S.MCC.fwhm_freq_pnts)
                                        if f == 1
                                            startIndex = 1;
                                        else
                                            startIndex = cumPnts(f-1) + 1;
                                        end
                                        endIndex = cumPnts(f);
                                        indices = startIndex:endIndex;
                                        resid_temp = resid(:,:,indices,:);
                                        resid_smooth_temp = zeros(size(resid_temp));
                                        spm_smooth(resid_temp, resid_smooth_temp, fwhm_voxels);
                                        resid_smooth(:,:,indices,:) = resid_smooth_temp;
                                    end
                                else
                                    spm_smooth(resid, resid_smooth, fwhm_voxels);
                                end
                                resid = resid_smooth .* double(mask_img); % Mask it
                                resid(isnan(resid)) = 0; % No NaN for smoothness function
                            end
                    
                            % Store in 4D array
                            resid_vol(:, :, :, tr) = resid;  % Store the processed volume
                    
                        % Method 3: Load residuals from MATFILE (large datasets)
                        elseif S.MCC.load_resid_as_matfile

                            if tr==1
                                % Define directory for saving residual files for Method 3
                                resid_vol_dir = fullfile(D(d).model(i).resid_vol_dir);
                                if ~exist(resid_vol_dir, 'dir')
                                    mkdir(resid_vol_dir);
                                end
                            end
                            % Determine chunk to load based on current tr
                            chunk_checkpoints = [1, ceil(szR / 2)];
                            if ismember(tr, chunk_checkpoints)
                                if tr == 1
                                    tempmat = vrmat.resid(1:ceil(szR / 2) - 1, :);
                                    tr_sub = 0;
                                elseif tr == ceil(szR / 2)
                                    clear tempmat
                                    tempmat = vrmat.resid(ceil(szR / 2):szR, :);
                                    tr_sub = ceil(szR / 2) - 1;
                                end
                            end
                    
                            resid(sidx) = tempmat(tr - tr_sub, :); % Load residuals into image
                    
                            % If topotime data
                            if ndims(resid) == 2
                                resid = topotime_3D(resid, S);
                            end
                    
                            % Apply smoothing
                            resid_smooth = zeros(size(resid));
                            if isfield(S.MCC,'fwhm_freq_pnts') && ~isempty(S.MCC.fwhm_freq_pnts)
                                % smooth within frequency bands separately
                                cumPnts = cumsum(S.MCC.fwhm_freq_pnts);
                                for f = 1:length(S.MCC.fwhm_freq_pnts)
                                    if f == 1
                                        startIndex = 1;
                                    else
                                        startIndex = cumPnts(f-1) + 1;
                                    end
                                    endIndex = cumPnts(f);
                                    indices = startIndex:endIndex;
                                    resid_temp = resid(:,:,indices,:);
                                    resid_smooth_temp = zeros(size(resid_temp));
                                    spm_smooth(resid_temp, resid_smooth_temp, fwhm_voxels);
                                    resid_smooth(:,:,indices,:) = resid_smooth_temp;
                                end
                            else
                                spm_smooth(resid, resid_smooth, fwhm_voxels);
                            end
                            resid = resid_smooth .* double(mask_img); % Mask it
                            resid(isnan(resid)) = 0; % No NaN for smoothness function
                    
                            % Write image to disk
                            V_local = V;  % Create a local copy of V structure
                            V_local.fname = fullfile(resid_vol_dir, sprintf('resid%06d.nii', tr));  % Set filename
                            spm_write_vol(V_local, resid);  % Write smoothed data to file
                        end
                    end


                    
                    % % save as multiple volumes
                    % pcc=0;
                    % for tr = 1:szR
                    % 
                    %     % monitor progress
                    %     pc = floor(tr*100/szR);
                    %     if pcc<pc
                    %         pcc=pc;
                    %         disp(['writing resid volumes: ' num2str(pc) '%'])
                    %     end
                    % 
                    %     % filename
                    %     trn = sprintf( '%06d', tr );
                    %     V.fname = fullfile(D(d).model(i).resid_vol_dir,['resid' trn '.nii']);
                    % 
                    %     % method
                    %     if exist('resid_vol','var')
                    %         if ndim_resid ==4
                    %             spm_write_vol(V,resid_vol(:,:,:,tr));
                    %         elseif ndim_resid ==3
                    %             spm_write_vol(V,resid_vol(:,:,tr));
                    %         end
                    %     elseif ~isempty(vr) && ndims(vr.resid)==2
                    %         % SOURCE DATA from condor: create empty NaN image template and it's indices; copy it; add indices of all samples within
                    %         % image N as lookup matrix; find indices of rows of
                    %         % resid (from resid_samples) corresponding to image N
                    %         % indices (intersect indices?); load each in turn and index into the empty image template. 
                    %         resid = resid_empty;
                    %         if tr==1 
                    %             sidx = find(~isnan(D(d).model(i).s(:)))'; % sample indices within image
                    %             if length(sidx)~= maxs % check
                    %                 error('wrong number of samples in resid data');
                    %             end
                    %         end
                    %         resid(sidx) = vr.resid(tr,:); % 
                    % 
                    %         % if topotime data
                    %         if ndims(resid) == 2
                    %             resid = topotime_3D(resid,S);
                    %         end
                    % 
                    %         % smooth
                    %         headsize=200; % mm
                    %         timeres=1;
                    %         szhead = size(resid);
                    %         voxel_size = [headsize/szhead(1) headsize/szhead(2) timeres]; %in mm, mm, ms 
                    %         desired_fwhm = [8 8 0];  % Desired smoothing in mm, mm, ms
                    %         fwhm_voxels = desired_fwhm ./ voxel_size;
                    %         resid_smooth=zeros(size(resid));
                    %         spm_smooth(resid,resid_smooth,fwhm_voxels);
                    % 
                    %         resid = resid.*double(mask_img); % mask it
                    %         resid(isnan(resid))= 0; % must be no NaN for smoothness function
                    % 
                    %         % write image
                    %         spm_write_vol(V,resid.*double(mask_img));
                    %     elseif S.MCC.load_resid_as_matfile % MATFILE VERSION
                    %         % SOURCE DATA from condor
                    %         resid = resid_empty;
                    %         if tr==1 
                    %             sidx = find(~isnan(D(d).model(i).s(:)))'; % sample indices within image
                    %             if length(sidx)~= maxs % check
                    %                 error('wrong number of samples in resid data');
                    %             end
                    %         end
                    % 
                    %         % for now, chunk into two halves. Later,
                    %         % chunking size can be a setting.
                    %         chunk_checkpoints = [1,ceil(szR/2)];
                    %         if ismember(tr, chunk_checkpoints)
                    %             if tr==1
                    %                 tempmat = vrmat.resid(1:ceil(szR/2)-1,:);
                    %                 tr_sub = 0;
                    %             elseif tr==ceil(szR/2)
                    %                 clear tempmat
                    %                 tempmat = vrmat.resid(ceil(szR/2):szR,:);
                    %                 tr_sub = ceil(szR/2)-1;
                    %             end
                    %         end
                    % 
                    %         resid(sidx) = tempmat(tr-tr_sub,:); % 
                    % 
                    %         % if topotime data
                    %         if ndims(resid) == 2
                    %             resid = topotime_3D(resid,S);
                    %         end
                    % 
                    %         resid = resid.*double(mask_img); % mask it
                    %         resid(isnan(resid))= 0; % must be no NaN for smoothness function
                    % 
                    %         % write image
                    %         spm_write_vol(V,resid.*double(mask_img));
                    %     end
                    % 
                    % 
                    % end
                end
                clear vr vrmat tempmat
                if ~exist('resid_vol', 'var')
                    fnames = dir(fullfile(D(d).model(i).resid_vol_dir,'*.nii'));
                    resid_vol = fullfile(D(d).model(i).resid_vol_dir,{fnames.name});
                end
                
                % smoothness estimation
                disp('MCC: smoothness estimation')
                
%                 % run one or multiple estimations?
%                 mask_dims = 1:ndims(mask_img);
%                 rep_dim = mask_dims(~ismember(mask_dims,S.MCC.dim));
%                 if isempty(rep_dim)
                    [L,Ve,ssq,Dx,Ix,Iy,Iz] = resid_smoothness(resid_vol,VM);
                    D(d).model(i).smooth.L = L;
                    D(d).model(i).smooth.Ve = Ve;
                    D(d).model(i).smooth.ssq = ssq;
                    D(d).model(i).smooth.Dx = Dx;
                    D(d).model(i).smooth.Ix = Ix;
                    D(d).model(i).smooth.Iy = Iy;
                    D(d).model(i).smooth.Iz = Iz;
%                 else
%                     nrep = size(mask_img,rep_dim);
%                     for nr = 1:nrep
%                         [L,Ve,ssq,Dx,Ix,Iy,Iz] = resid_smoothness(resid_vol,VM);
%                     end
%                 end
            end

            corrZ = nan(1,length(S.MCC.model.contrast{i}));
            for c = S.MCC.model.contrast{i}
                disp(['MCC: model ' num2str(i) 'contrast ' num2str(c)])

                % image names
                D(d).model(i).con(c).MCC.mask_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_MCC_mask.nii']);
                D(d).model(i).con(c).MCC.p_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_MCC_p.nii']);
                D(d).model(i).con(c).MCC.Z_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_MCC_Z.nii']);
                D(d).model(i).con(c).MCC.mask_p_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_MCC_mask_p.nii']);
                D(d).model(i).con(c).MCC.mask_Z_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_MCC_mask_Z.nii']);

                % if already created, load images from disk
                if exist(D(d).model(i).con(c).MCC.mask_Z_img_file,'file') && ~S.MCC.estimate
                    mask_img=spm_read_vols(spm_vol(D(d).model(i).con(c).MCC.mask_img_file));
                    Z_img=spm_read_vols(spm_vol(D(d).model(i).con(c).MCC.Z_img_file));
                    p_img=spm_read_vols(spm_vol(D(d).model(i).con(c).MCC.p_img_file));
                    mask_Z_img=spm_read_vols(spm_vol(D(d).model(i).con(c).MCC.mask_Z_img_file));
                    mask_p_img=spm_read_vols(spm_vol(D(d).model(i).con(c).MCC.mask_p_img_file));
                elseif ~S.MCC.estimate
                    error('images do not exist - estimation required')
                else
                    DF = D(d).model(i).con(c).DF;
                    if strcmp (S.MCC.DF_type,'residualP')
                        DF(2) = size(D(d).model(i).fixeddesign,1)-size(D(d).model(i).fixeddesign,2);
                    elseif strcmp (S.MCC.DF_type,'residualPR')
                        DF(2) = size(D(d).model(i).fixeddesign,1)-size(D(d).model(i).fixeddesign,2)-size(D(d).model(i).randomdesign,2);
                    elseif strcmp (S.MCC.DF_type,'halfway')
                        DF(2) = (size(D(d).model(i).fixeddesign,1)-size(D(d).model(i).fixeddesign,2) + DF(2)) / 2;
                    end
                    p_img = spm_read_vols(spm_vol(D(d).model(i).con(c).p_img_file));
                    szR = length(D(d).model(i).fixeddesign);
                    if S.MCC.use_pTFCE || strcmp(S.MCC.method,'FWE_RF')
                        [FWHM,VRpv,R] = smoothness_df_correct(...
                            D(d).model(i).smooth.L,...
                            D(d).model(i).smooth.Ve,...
                            D(d).model(i).smooth.ssq,...
                            D(d).model(i).smooth.Dx,...
                            D(d).model(i).smooth.Ix,...
                            D(d).model(i).smooth.Iy,...
                            D(d).model(i).smooth.Iz,...
                            VM,[szR DF(2)]);
                        D(d).model(i).con(c).MCC.FWHM = FWHM;
                        D(d).model(i).con(c).MCC.VRpv = VRpv;
                        D(d).model(i).con(c).MCC.R = R;

                        % calculate smoothness without any adjustment for
                        % DF, for reference
                        [FWHM_noDF,VRpv_noDF,R_noDF] = smoothness_df_correct(...
                            D(d).model(i).smooth.L,...
                            D(d).model(i).smooth.Ve,...
                            D(d).model(i).smooth.ssq,...
                            D(d).model(i).smooth.Dx,...
                            D(d).model(i).smooth.Ix,...
                            D(d).model(i).smooth.Iy,...
                            D(d).model(i).smooth.Iz,...
                            VM,[szR szR]);
                        D(d).model(i).con(c).MCC.FWHM_noDF = FWHM_noDF;
                        D(d).model(i).con(c).MCC.VRpv_noDF = VRpv_noDF;
                        D(d).model(i).con(c).MCC.R_noDF = R_noDF;

                        % ensure resel count is between 1 and num voxels
                        nv = nansum(mask_img(:)); % Number of voxels in mask (e.g. SPM.xVol.S)
                        % R:  Resel count, e.g. SPM.xVol.R(4).
                        minR = min(max(1,R(4)),nv);
                        D(d).model(i).con(c).MCC.minR = minR;
                    end
                    
                    
                    % convert it to Z-score (from T or F)
                    % Z_img: Z-score image to enhance
                    if isfield(D(d).model(i).con(c),'T')
                        T_img = spm_read_vols(spm_vol(D(d).model(i).con(c).T_img_file));
                        % this might be wrong - see: https://github.com/spisakt/pTFCE/issues/9
                        Z_img = img_t2z(T_img, DF(2), S.MCC.tails);
                    elseif isfield(D(d).model(i).con(c),'F')
                        F_img = spm_read_vols(spm_vol(D(d).model(i).con(c).F_img_file));
                        if any(F_img(:)<0)
                            F_img(F_img(:)<0)=0;
                        end
                        Z_img = img_f2z(F_img, DF(1), DF(2), 1);
                    else
                        error('Statistical parameter type not supported');
                    end
                    
                    if S.MCC.use_pTFCE
                        [pTFCE_Z, pTFCE_p] = pTFCE(Z_img,mask_img,minR,nv);
                        pTFCE_Z(mask_img==0 | isnan(mask_img))=NaN;
                        pTFCE_p(mask_img==0 | isnan(mask_img))=NaN;
%                         figure; scatter(1:numel(pTFCE_Z),sort(pTFCE_Z(:)),'r'); hold on; scatter(1:numel(Z_img),sort(Z_img(:)),'b'); ylabel('Z scores');
%                         figure; scatter(Z_img(:),pTFCE_Z(:)); xlabel('original Z');ylabel('pTFCE Z');
                        Z_img=pTFCE_Z;
                        p_img=pTFCE_p;
                    end
                    zsize=size(Z_img);

                    switch S.MCC.method
                        case 'FWE_RF'
                            % FWE corrected critical height threshold at specified significance level
                            corrZ(c) = spm_uc_RF(S.MCC.thresh,DF,'Z',minR,1);
                            D(d).model(i).con(c).MCC.zthresh = corrZ(c);
                            switch length(zsize)
                                case 2
                                    MCCmask_img = reshape(Z_img>=corrZ(c),zsize(1),zsize(2));
                                case 3
                                    MCCmask_img = reshape(Z_img>=corrZ(c),zsize(1),zsize(2),zsize(3));
                            end

                        case 'FDR'
                            if strcmp(S.MCC.FDR_type,'para')
                                [fdr_p] = FDR(p_img(:),S.MCC.thresh);
                            elseif strcmp(S.MCC.FDR_type,'nonpara')
                                [~,fdr_p] = FDR(p_img(:),S.MCC.thresh);
                            elseif strcmp(S.MCC.FDR_type,'twostep')
                                % currently only implements two-sided test 
                                mask_img0=mask_img(:);
                                mask_img0(isnan(mask_img0))=0;
                                [~,~,fdr_p,~] = FDR2(p_img(:),ones(size(p_img(:))),find(mask_img0),S.MCC.thresh,0); % pval,sgn,maskvtx,rate,tail
                            end
                            corrZ(c) = norminv(1 - fdr_p); % corrected threshold
                            psize=size(p_img);
                            switch length(psize)
                                case 1
                                    MCCmask_img = p_img<fdr_p;
                                case 2
                                    MCCmask_img = reshape(p_img<fdr_p,psize(1),psize(2));
                                case 3
                                    MCCmask_img = reshape(p_img<fdr_p,psize(1),psize(2),psize(3));
                            end  
                        case 'none'
                            psize=size(p_img);
                            switch length(psize)
                                case 1
                                    MCCmask_img = p_img<S.MCC.thresh;
                                case 2
                                    MCCmask_img = reshape(p_img<S.MCC.thresh,psize(1),psize(2));
                                case 3
                                    MCCmask_img = reshape(p_img<S.MCC.thresh,psize(1),psize(2),psize(3));
                            end
                    end

                    % cluster extent masking
                    % pTFCE code: do connected component analysis for the whole set of threesholds
                    cc = arrayfun(@(x) bwconncomp(bsxfun(@ge,Z_img,x),6), corrZ(c));
                    CLUST=zeros(size(Z_img));
                    ccc = cc(1);
                    voxpercc = cellfun(@numel,ccc.PixelIdxList);
                    for ci = 1:ccc.NumObjects
                        CLUST(ind2sub(size(Z_img), ccc.PixelIdxList{ci})) = voxpercc(ci);
                    end
                    MCCmask_img(CLUST<S.MCC.extent_thresh) = 0;

                    if ~isempty(S.MCC.mask_img_post)
                        temp = spm_read_vols(spm_vol(S.MCC.mask_img_post))>0;
                        MCCmask_img=MCCmask_img.*temp;
                    end
                    mask_Z_img = Z_img.*MCCmask_img;
                    pMCCmask_img=double(MCCmask_img);
                    pMCCmask_img(MCCmask_img==0)=NaN;
                    mask_p_img = p_img.*pMCCmask_img;

                    % save images
                    V.fname = D(d).model(i).con(c).MCC.mask_img_file;
                    spm_write_vol(V,MCCmask_img);
                    V.fname = D(d).model(i).con(c).MCC.p_img_file;
                    spm_write_vol(V,p_img);
                    V.fname = D(d).model(i).con(c).MCC.Z_img_file;
                    spm_write_vol(V,Z_img);
                    V.fname = D(d).model(i).con(c).MCC.mask_p_img_file;
                    spm_write_vol(V,mask_p_img);
                    V.fname = D(d).model(i).con(c).MCC.mask_Z_img_file;
                    spm_write_vol(V,mask_Z_img);

                    % convert back to original statistics (T or F)
                    if isfield(D(d).model(i).con(c),'T')
                        T_img = tinv(1-p_img, DF);
                        D(d).model(i).con(c).MCC.T_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_MCC_T.nii']);
                        V.fname = D(d).model(i).con(c).MCC.T_img_file;
                        spm_write_vol(V,T_img);
                        D(d).model(i).con(c).MCC.mask_T_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_MCC_mask_T.nii']);
                        V.fname = D(d).model(i).con(c).MCC.mask_T_img_file;
                        spm_write_vol(V,T_img.*MCCmask_img);
                    elseif isfield(D(d).model(i).con(c),'F')
                        F_img = finv(1-p_img, DF(1), DF(2));
                        D(d).model(i).con(c).MCC.F_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_MCC_F.nii']);
                        V.fname = D(d).model(i).con(c).MCC.F_img_file;
                        spm_write_vol(V,F_img);
                        D(d).model(i).con(c).MCC.mask_F_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_MCC_mask_F.nii']);
                        V.fname = D(d).model(i).con(c).MCC.mask_F_img_file;
                        spm_write_vol(V,F_img.*MCCmask_img);
                    end
                end

                % plot data
                plotdatY(:,c) = Z_img(:);
                plotdatX(:,c) = c*double(CLUST(:)>0);
                ClustX(:,c) = CLUST(:);
                xtl{c} = D(d).model(i).con(c).term;

                % save clus info
                D(d).model(i).con(c).MCCclus_nvox = sum(unique(CLUST(CLUST>=S.MCC.extent_thresh)));
            end

            % plot
            figname = ['model ' num2str(i) ' post-threshold Z values'];
            xlab = 'contrast';
            ylab = 'Z values';
            uncorrz = norminv(1 - S.MCC.thresh/S.MCC.tails);
            try
                % f=jitterplot(figname,plotdatX,plotdatY,xlab,ylab,xtl,[uncorrz, corrZ]);
                f=jitterplotclust(figname,plotdatX,plotdatY,xlab,ylab,xtl,[uncorrz, corrZ],ClustX,S.MCC.extent_thresh);
                saveas(f,fullfile(S.MCC.path.inputs,['Z_values_model' num2str(i) '.png']))
            catch
                disp(['no significant effects for model ' num2str(i)])
            end
            clear plotdatX plotdatY xtl resid_vol ClustX
        end
    end

    if S.MCC.model_comp.index
        for i = S.MCC.model_comp.index

            % image names
            D(d).model_comp(i).MCC.mask_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_modelcomp_' num2str(i) '_MCC_mask.nii']);
            D(d).model_comp(i).MCC.p_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_modelcomp_' num2str(i) '_MCC_p.nii']);
            D(d).model_comp(i).MCC.Z_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_modelcomp_' num2str(i) '_MCC_Z.nii']);
            D(d).model_comp(i).MCC.mask_p_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_modelcomp_' num2str(i) '_MCC_mask_p.nii']);
            D(d).model_comp(i).MCC.mask_Z_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_modelcomp_' num2str(i) '_MCC_mask_Z.nii']);
            D(d).model_comp(i).MCC.mask_LR_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_modelcomp_' num2str(i) '_MCC_mask_LR.nii']);

            % if already created, load images from disk
            if exist(D(d).model_comp(i).MCC.mask_LR_img_file,'file') && ~S.MCC.estimate
                mask_img=spm_read_vols(spm_vol(D(d).model_comp(i).MCC.mask_img_file));
                p_img=spm_read_vols(spm_vol(D(d).model_comp(i).MCC.p_img_file));
                mask_p_img=spm_read_vols(spm_vol(D(d).model_comp(i).MCC.mask_p_img_file));
                mask_LR_img=spm_read_vols(spm_vol(D(d).model_comp(i).MCC.mask_LR_img_file));
            else

                % mask
                contrast_Z_img=[];
                for p = 1:2
                    for c = 1:length(D(d).model(D(d).model_comp(i).pair(p)).con)
                        contrast_Z_img(:,:,:,p,c)=spm_read_vols(spm_vol(D(d).model(D(d).model_comp(i).pair(p)).con(c).MCC.mask_Z_img_file));
                    end
                end
                contrast_mask = any(contrast_Z_img,[4 5]);
                MCCmask_img = mask_img.*contrast_mask;

                % CONVERT TO Z-SCORE
                p_img = spm_read_vols(spm_vol(D(d).model_comp(i).pval_img_file));
                LR_img = spm_read_vols(spm_vol(D(d).model_comp(i).LR_img_file));
                
                % use LR threshold
                LR_img_thresh = LR_img>S.MCC.LR_thresh | LR_img<(1/S.MCC.LR_thresh);
                psize=size(p_img);
                switch length(psize)
                    case 2
                        MCCmask_img = reshape(LR_img_thresh,psize(1),psize(2));
                    case 3
                        MCCmask_img = reshape(LR_img_thresh,psize(1),psize(2),psize(3));
                end
                    
                if ~isempty(S.MCC.mask_img_post)
                    temp = spm_read_vols(spm_vol(S.MCC.mask_img_post))>0;
                    MCCmask_img=MCCmask_img.*temp;
                end
                mask_LR_img = LR_img.*MCCmask_img;
                %mask_Z_img = Z_img.*MCCmask_img;
                
                pMCCmask_img=double(MCCmask_img);
                pMCCmask_img(MCCmask_img==0)=NaN;
                mask_p_img = p_img.*pMCCmask_img;

                % save images
                V.fname = D(d).model_comp(i).MCC.mask_img_file;
                spm_write_vol(V,MCCmask_img);
                V.fname = D(d).model_comp(i).MCC.p_img_file;
                spm_write_vol(V,p_img);
                V.fname = D(d).model_comp(i).MCC.mask_p_img_file;
                spm_write_vol(V,mask_p_img);
                V.fname = D(d).model_comp(i).MCC.mask_LR_img_file;
                spm_write_vol(V,mask_LR_img);
            end
            % plot data
            plotdatY(:,i) = mask_LR_img(:);
            plotdatX(:,i) = i*double(plotdatY(:,i)>0);
            xtl{i} = num2str(D(d).model_comp(i).pair);
        end
        figname = 'model comparison post-threshold LR values';
        xlab = 'models compared';
        ylab = 'likelihood ratios';
        try
            jitterplot(figname,plotdatX,plotdatY,xlab,ylab,xtl,[1/20 20])
        catch
            disp(['no significant effects for model comp ' num2str(i)])
        end
        ylab = 'log likelihood ratios';
        try
            jitterplot(figname,plotdatX,log(plotdatY),xlab,ylab,xtl,[1/20 20])
        catch
            disp(['no significant effects for model comp ' num2str(i)])
        end
        clear plotdatX plotdatY xtl
    end

end
save(fullfile(S.MCC.path.outputs, 'D.mat'),'D','-v7.3');

function set_paths(S)
for p = 1:size(S.path.code,1)
    if S.path.code{p,1}
        addpath(genpath(S.path.code{p,2}));
    else
        addpath(S.path.code{p,2});
    end
end

function f = jitterplotclust(figname, plotdatX, plotdatY, xlab, ylab, xtl, thresh, CLUST, cthresh)
    % Ensure plotdatX, plotdatY, and CLUST are matrices of the same size
    [m, n] = size(plotdatY);
    assert(all(size(plotdatX) == [m, n]), 'plotdatX must be the same size as plotdatY');
    assert(all(size(CLUST) == [m, n]), 'CLUST must be the same size as plotdatY');
    
    meanX = [];
    meanY = [];
    clusterSizes = [];
    clusterIDs = []; % Store cluster IDs
    xGroups_used = [];
    
    % For each column (group on the x-axis)
    for xi = 1:n
        % Extract xGroup value from plotdatX (assumes same value in column for valid data)
        xGroup_col = plotdatX(:, xi);
        % Valid indices where plotdatX and CLUST are non-zero and not NaN
        valid_indices = (xGroup_col ~= 0) & ~isnan(xGroup_col) & (CLUST(:, xi) ~= 0) & ~isnan(CLUST(:, xi));
        if ~any(valid_indices)
            continue; % No valid data in this column, skip to next column
        end
        
        xGroup = xGroup_col(find(valid_indices, 1)); % Get xGroup value (should be same for all valid data)
        xGroups_used = [xGroups_used; xGroup]; % Store used xGroups
        
        % Get valid data for this column
        plotdatY_col = plotdatY(valid_indices, xi);
        CLUST_col = CLUST(valid_indices, xi);
        
        % Find unique clusters within this column
        clustersInXGroup = unique(CLUST_col);
        
        for ci = 1:length(clustersInXGroup)
            clusterID = clustersInXGroup(ci);
            cluster_indices = (CLUST_col == clusterID);
            
            % Calculate mean of plotdatY for this cluster
            meanYVal = mean(plotdatY_col(cluster_indices));
            clusterSize = sum(cluster_indices);
            
            % Store the results
            meanX = [meanX; -xGroup]; % Negate xGroup to match original plotting style
            meanY = [meanY; meanYVal];
            clusterSizes = [clusterSizes; clusterSize];
            clusterIDs = [clusterIDs; clusterID]; % Store clusterID
        end
    end
    
    % Check if there are any data points to plot
    if isempty(meanX)
        warning('No data points to plot after processing.');
        f = figure('name', figname); % Return an empty figure
        return;
    end
    
    % Normalize clusterSizes to determine marker sizes
    maxMarkerSize = 1000; % Adjust this value as needed
    markerSizes = (clusterSizes / max(clusterSizes)) * maxMarkerSize;
    
    % Determine colors based on clusterSizes and cthresh
    numPoints = length(meanX);
    colors = repmat([0 0 0], numPoints, 1); % Default color is black
    red_indices = clusterSizes < cthresh;
    colors(red_indices, :) = repmat([1 0 0], sum(red_indices), 1); % Set red color where clusterSize < cthresh
    
    % Plotting
    f = figure('name', figname);
    hold on
    s = scatter(meanX, meanY, markerSizes, colors, 'filled', 'jitter', 'on', 'jitterAmount', 0.3);
    s.MarkerFaceAlpha = 0.5;
    
    % Add threshold lines
    lineColors = {'red', 'green'};
    xGroups_used = unique(xGroups_used);
    xLimits = [-(max(xGroups_used) + 0.5), -0.5];
    line(xLimits, [thresh(1) thresh(1)], 'Color', lineColors{1}, 'LineStyle', '--');
    for t = 2:length(thresh)
        if ~isnan(thresh(t)) && ~isinf(thresh(t))
            line(xLimits, [thresh(t) thresh(t)], 'Color', lineColors{2}, 'LineStyle', '--');
        end
    end
    hold off
    
    % Axis settings
    xlim(xLimits)
    xticks(sort(-xGroups_used))
    % Map xGroups_used to indices for xtl
    xtl_indices = xGroups_used;
    % Ensure xtl_indices are valid indices for xtl
    xtl_indices(xtl_indices > length(xtl)) = length(xtl);
    xticklabels(fliplr(xtl(xtl_indices)))
    xlabel(xlab)
    ylabel(ylab)
    view([90 -90])



function f=jitterplot(figname,plotdatX,plotdatY,xlab,ylab,xtl,thresh)
plotdatX=-plotdatX;
plotdatX(plotdatX==0)=NaN;
f=figure('name',figname);
minmax = -[size(plotdatX,2)+0.5,0.5];
hold on
s=scatter(plotdatX(:), plotdatY(:),10,'k','filled', 'jitter','on', 'jitterAmount',0.3);
% s=scatter(plotdatX(:), plotdatY(:),10,'k','filled');
s.MarkerFaceAlpha=0.5;
%if strcmp(ylab,'Z values')
    colors = {'red','green'};
    % uncorrected
    line(minmax,[thresh(1) thresh(1)],'Color',colors{1},'LineStyle','--') 
    for t = 2:length(thresh)
        if ~isnan(thresh(t)) && ~isinf(thresh(t))
            minmaxt = -[(t-1)+0.5,0.5];
            line(minmaxt,[thresh(t) thresh(t)],'Color',colors{2},'LineStyle','--') 
        end
    end
%end
hold off
xlim(minmax)
xticks(-size(plotdatX,2):-1)
xticklabels(fliplr(xtl))
xlabel(xlab)
ylabel(ylab)
view([90 -90])

function V = makeV
V.fname = '';
V.dim = [1 1 1];
V.dt = [16,0];
V.pinfo = [1;0;352];
V.mat = [-4.25000000000000,0,0,68;0,5.37500000000000,0,-100;0,0,-1,-201;0,0,0,1];
V.n = [1,1];
V.descrip = '';
% V.private

function V = makeV3
% makes 3D volume in MNI space
V.fname = '';
V.dim = [1 1 1];
V.dt = [16,0];
V.pinfo = [1;0;32];
V.mat = [-2,0,0,92;0,2,0,-128;0,0,2,-74;0,0,0,1];
V.n = [1,1];
V.descrip = '';
% V.private
