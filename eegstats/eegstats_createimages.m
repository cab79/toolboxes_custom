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
        mask_img = spm_read_vols(spm_vol(S.img.mask_img))>0;
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
    V=spm_vol(S.img.file.SPMimg);
    if S.img.model.index
        for i = S.img.model.index
            disp(['converting to images for ' save_pref num2str(d) ', model ' num2str(i)])
            
            
            
            %% betas
            for b = 1:length(D(d).model(i).coeff)
                % create
                if ~isempty(S.img.mask_img)
                    [b_img] = topotime_3D(D(d).model(i).coeff(b).b,S);
                else
                    [b_img, mask_img] = topotime_3D(D(d).model(i).coeff(b).b,S);
                end
                b_img = topotime_3D(D(d).model(i).coeff(b).b,S);
                % save images
                D(d).model(i).coeff(b).b_img_file = fullfile(S.img.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_coeff_' num2str(b) '_b.nii']);
                V.fname = D(d).model(i).coeff(b).b_img_file;
                V.dim = size(mask_img);
                spm_write_vol(V,b_img.*mask_img);
            end
            
            %% model quantities
            Q = {'logl','r2_ord','waic'};
            for q = 1:length(Q)
                if isfield(D(d).model(i),Q{q})
                    % create
                    img = topotime_3D(D(d).model(i).(Q{q}),S);
                    % save
                    V.dim = size(mask_img);
                    D(d).model(i).([Q{q} '_img_file']) = fullfile(S.img.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_' Q{q} '.nii']);
                    V.fname = D(d).model(i).([Q{q} '_img_file']);
                    spm_write_vol(V,img);
                end
            end
            
            %% contrasts
            if isfield(D(d).model(i),'con')
                for c = S.img.model.contrast{i}
                    % create images
                    [F_img] = topotime_3D(D(d).model(i).con(c).F,S);
                    [p_img] = topotime_3D(D(d).model(i).con(c).p,S);
                    V.dim = size(mask_img);
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
            spm_write_vol(V,mask_img);
        end
    end

    %% model comparisons
    if S.img.model_comp.index
        for i = S.img.model_comp.index
            % create
            [p_img] = topotime_3D(D(d).model_comp(i).pval,S);
            if ~isempty(S.img.mask_img)
                [LR_img] = topotime_3D(D(d).model_comp(i).LR,S);
            else
                [LR_img,mask_img] = topotime_3D(D(d).model_comp(i).LR,S);
            end
            [LR_p_img] = topotime_3D(D(d).model_comp(i).LR_p,S);
            % save
            V.dim = size(mask_img);
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
save(fullfile(S.img.path.outputs, 'D.mat'),'D');

function set_paths(S)
for p = 1:size(S.path.code,1)
    if S.path.code{p,1}
        addpath(genpath(S.path.code{p,2}));
    else
        addpath(S.path.code{p,2});
    end
end


