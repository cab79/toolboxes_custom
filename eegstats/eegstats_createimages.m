function D=eegstats_createimages(S,varargin)
dbstop if error

%% find encoding models
if isempty(varargin)
   % try loading the data if not supplied in varargin
   try
       load(fullfile(S.img.path.inputs,S.img.file.inputs),'D') 
   catch
       error('Supply a D structure or path/file of the data')
   end
else
    D= varargin{1};
end

% set paths
S.path.code = {
    1, 'E:\Matlab_backup_moved_to_Q\spm12' % SPM
    };
set_paths(S)
if ~exist(S.img.path.outputs,'dir')
    mkdir(S.img.path.outputs)
end
cd(S.img.path.outputs)

if length(D)>1
    save_pref = 'subject_';
else
    save_pref = 'group_';
end

if ~isfield(S.img,'desired_fwhm')
    S.img.desired_fwhm = [0 0 0];  % Desired smoothing in mm, mm, ms
end
if ~isfield(S.img,'headsize')
    S.img.headsize = 200;
end

for d=1:length(D)
    % get model and contrast indices
    if isempty(S.img.model.index)
        S.img.model.index = 1:length(D(d).model);
    end
    if isempty(S.img.model.contrast) && isfield(D(d).model(1),'con')
        for i = S.img.model.index
            S.img.model.contrast{i} = 1:length(D(d).model(i).con);
        end
    end
    if isempty(S.img.model_comp.index)
        try
            S.img.model_comp.index = 1:length(D(d).model_comp);
        catch
            S.img.model_comp.index = 0;
        end
    end

%     % design
%     if isfield(in.stats.LME,'design')
%         D(d).design = in.stats.LME.design;
%     else
%         temp=load(fullfile(S.img.path.stats_load,'example_input.mat'));
%         D(d).design=temp.Y.dtab;
%     end

    % model stats
    if isfield(D(d).model(1),'con')
        if S.img.model.index
            for i = S.img.model.index
    %             D(d).model(i).def = in.S.lme_model(i); 
                disp(['model: ' D(d).model(i).def])

                % take max of DF
                for c = S.img.model.contrast{i}
                    DF(1) = max(max(squeeze(D(d).model(i).con(c).DF(1,:,:))));
                    DF(2) = max(max(squeeze(D(d).model(i).con(c).DF(2,:,:))));
                    D(d).model(i).con(c).DF = DF;
                end
            end
        end
    end

    % model comparison stats
    if S.img.model_comp.index
        for i = S.img.model_comp.index
            D(d).model_comp(i).LR_p = D(d).model_comp(i).pval<=S.img.pthresh .* D(d).model_comp(i).LR;
            %disp(['model comparison: ' num2str(D(d).model_comp(i).pair{:})])
        end
    end

    if ~isempty(S.img.mask_img)
        mask_img = double(spm_read_vols(spm_vol(S.img.mask_img))>0);
        mask_img(mask_img==0) = nan;
    end
    
    % change file locations for resid, fitted, input
    
    if S.img.model.index
        for i = S.img.model.index
            [~,residfile,ext]=fileparts(D(d).model(i).resid_file);
            D(d).model(i).resid_file = fullfile(S.img.path.inputs,[residfile ext]);
            [~,fittedfile,~]=fileparts(D(d).model(i).fitted_file);
            D(d).model(i).fitted_file = fullfile(S.img.path.inputs,[fittedfile ext]);
            [~,inputfile,~]=fileparts(D(d).model(i).input_file);
            D(d).model(i).input_file = fullfile(S.img.path.inputs,[inputfile ext]);
        end
    end

    %% convert to images
    if S.img.model.index
        for i = S.img.model.index
            disp(['converting to images for ' save_pref num2str(d) ', model ' num2str(i)])
            
            
            
            %% betas
            for b = 1:length(D(d).model(i).coeff)
                % create
                if ~isempty(S.img.file.coord)
                    if ~isempty(S.img.mask_img)
                        [b_img] = topotime_3D(D(d).model(i).coeff(b).b,S);
                    else
                        [b_img, mask_img] = topotime_3D(D(d).model(i).coeff(b).b,S);
                    end
                    b_img = smooth_images(b_img,mask_img,S.img.desired_fwhm,S.img.headsize,1);
                    %b_img = topotime_3D(D(d).model(i).coeff(b).b,S);
                    V = makeV;
                else
                    b_img=D(d).model(i).coeff(b).b;
                    if isempty(S.img.mask_img)
                        mask_img=ones(size(b_img));
                    end
                    V = makeV3; % assume MNI volumne
                end
                % save images
                D(d).model(i).coeff(b).b_img_file = fullfile(S.img.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_coeff_' num2str(b) '_b.nii']);
                V.fname = D(d).model(i).coeff(b).b_img_file;
                V.dim(1:length(size(mask_img))) = size(mask_img);
                spm_write_vol(V,b_img.*mask_img);
            end
            
            %% random effects - assumes full random slopes model
            if isfield(D(d).model(i),'random')
                randmat = D(d).model(i).random;
                rm_dim = ndims(randmat);
                D(d).model(i).random = struct;
                rU = unique(D(d).model(i).randomnames.Name);
                for b = 1:length(rU)
                    [U,~,Ui] = unique(D(d).model(i).randomnames.Level,'stable');
                    for u = 1:length(U)
                        idx = find(Ui==u & strcmp(D(d).model(i).randomnames.Name,rU{b}));
                        if length(idx)>1; error('wrong random count'); end
                        % create
                        if ~isempty(S.img.file.coord)
                            if ~isempty(S.img.mask_img)
                                [r_img] = topotime_3D(randmat(:,:,idx),S);
                            else
                                [r_img, mask_img] = topotime_3D(randmat(:,:,idx),S);
                            end
                            r_img = smooth_images(r_img,mask_img,S.img.desired_fwhm,S.img.headsize,1);
                            V = makeV;
                        else
                            if rm_dim==3 % sensor
                                r_img=randmat(:,:,idx);
                                V = makeV;
                            elseif rm_dim==4 % source
                                r_img=randmat(:,:,:,idx);
                                V = makeV3; % assume MNI volumne
                            end
                        end
                        % add to file
                        D(d).model(i).random(idx).Level = D(d).model(i).randomnames.Level{idx};
                        D(d).model(i).random(idx).Name = D(d).model(i).randomnames.Name{idx};
                        D(d).model(i).random(idx).b = r_img;
                        % save images
                        D(d).model(i).random(idx).b_img_file = fullfile(S.img.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_randcoeff_' num2str(b) '_' num2str(u) '.nii']);
                        V.fname = D(d).model(i).random(idx).b_img_file;
                        V.dim(1:length(size(mask_img))) = size(mask_img);
                        spm_write_vol(V,r_img.*mask_img);
                    end
                end
            end
            
            %% model quantities
            Q = {'logl','r2_ord','waic'};
            for q = 1:length(Q)
                if isfield(D(d).model(i),Q{q})
                    % create
                    if ~isempty(S.img.file.coord)
                        img = topotime_3D(D(d).model(i).(Q{q}),S);
                    else
                        img=D(d).model(i).(Q{q});
                    end
                    
                    % save
                    V.dim(1:length(size(mask_img))) = size(mask_img);
                    D(d).model(i).([Q{q} '_img_file']) = fullfile(S.img.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_' Q{q} '.nii']);
                    V.fname = D(d).model(i).([Q{q} '_img_file']);
                    spm_write_vol(V,img);
                end
            end
            
            %% contrasts
            if isfield(D(d).model(i),'con')
                for c = S.img.model.contrast{i}
                    % create images
                    if ~isempty(S.img.file.coord)
                        [F_img] = topotime_3D(D(d).model(i).con(c).F,S);
                        [p_img] = topotime_3D(D(d).model(i).con(c).p,S);
                        
                        % correct interpolation problem: creates negative p
                        % values
                        p_img(p_img<0)=0;
                    else
                        F_img=D(d).model(i).con(c).F;
                        p_img=D(d).model(i).con(c).p;
                    end
                    F_img = smooth_images(F_img,mask_img,S.img.desired_fwhm,S.img.headsize,1);
                    p_img = smooth_images(p_img,mask_img,S.img.desired_fwhm,S.img.headsize,1);
                    V.dim(1:length(size(mask_img))) = size(mask_img);
                    % save images
                    D(d).model(i).con(c).F_img_file = fullfile(S.img.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_F.nii']);
                    V.fname = D(d).model(i).con(c).F_img_file;
                    spm_write_vol(V,F_img.*mask_img);
                    D(d).model(i).con(c).p_img_file = fullfile(S.img.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_p.nii']);
                    V.fname = D(d).model(i).con(c).p_img_file;
                    spm_write_vol(V,p_img.*mask_img);
                end
            end

            % save image mask
            D(d).model(i).mask_img_file = fullfile(S.img.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_mask.nii']);
            V.fname = D(d).model(i).mask_img_file;
%             spm_write_vol(V,double(~isnan(D(d).model(i).s)) .* mask_img);
            spm_write_vol(V,double(mask_img));
        end
    end

    %% model comparisons
    if S.img.model_comp.index
        for i = S.img.model_comp.index
            % create
            if ~isempty(S.img.file.coord)
                [p_img] = topotime_3D(D(d).model_comp(i).pval,S);
                
                % correct interpolation problem: creates negative p
                % values
                p_img(p_img<0)=0;
                
                if ~isempty(S.img.mask_img)
                    [LR_img] = topotime_3D(D(d).model_comp(i).LR,S);
                else
                    [LR_img,mask_img] = topotime_3D(D(d).model_comp(i).LR,S);
                end
                [LR_p_img] = topotime_3D(D(d).model_comp(i).LR_p,S);
            else
                LR_img=D(d).model_comp(i).LR;
                p_img=D(d).model_comp(i).pval;
                LR_p_img=D(d).model_comp(i).LR_p;
            end
            p_img = smooth_images(p_img,mask_img,S.img.desired_fwhm,S.img.headsize,1);
            LR_img = smooth_images(LR_img,mask_img,S.img.desired_fwhm,S.img.headsize,1);
            
            % save
            V.dim(1:length(size(mask_img))) = size(mask_img);
            D(d).model_comp(i).pval_img_file = fullfile(S.img.path.outputs, [save_pref num2str(d) '_modelcomp_' num2str(i) '_p.nii']);
            V.fname = D(d).model_comp(i).pval_img_file;
            spm_write_vol(V,p_img.*mask_img);
            D(d).model_comp(i).LR_img_file = fullfile(S.img.path.outputs, [save_pref num2str(d) '_modelcomp_' num2str(i) '_LR.nii']);
            V.fname = D(d).model_comp(i).LR_img_file;
            spm_write_vol(V,LR_img.*mask_img);
            D(d).model_comp(i).LR_p_img_file = fullfile(S.img.path.outputs, [save_pref num2str(d) '_modelcomp_' num2str(i) '_LR_p.nii']);
            V.fname = D(d).model_comp(i).LR_p_img_file;
            spm_write_vol(V,LR_p_img.*mask_img);

            % save image mask
            D(d).model_comp(i).mask_img_file = fullfile(S.img.path.outputs, [save_pref num2str(d) '_modelcomp_' num2str(i) '_mask.nii']);
            V.fname = D(d).model_comp(i).mask_img_file;
            spm_write_vol(V,mask_img);
        end

    end
    
end
save(fullfile(S.img.path.outputs, 'D.mat'),'D','-v7.3');

function set_paths(S)
for p = 1:size(S.path.code,1)
    if S.path.code{p,1}
        addpath(genpath(S.path.code{p,2}));
    else
        addpath(S.path.code{p,2});
    end
end

function V = makeV
V.fname = '';
V.dim = [1 1 1];
V.dt = [16,0];
%V.pinfo = [0;0;0];
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

function image_smooth = smooth_images(image,mask_img,desired_fwhm,headsize,timeres)

if ~any(desired_fwhm>0)
    image_smooth = image;
end

% Determine voxel size for smoothing
szhead = size(mask_img);
voxel_size = [headsize / szhead(1), headsize / szhead(2), timeres]; % in mm, mm, ms 
fwhm_voxels = desired_fwhm ./ voxel_size;

% Apply smoothing
image_smooth = zeros(size(image));
spm_smooth(image, image_smooth, fwhm_voxels);
image_smooth = image_smooth .* double(mask_img); % Mask it
image_smooth(isnan(image_smooth)) = 0; % No NaN for smoothness function