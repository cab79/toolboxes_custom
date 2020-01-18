function D=eegstats_MCC(S,varargin)
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
    V=spm_vol(S.MCC.file.SPMimg);

    % mask
    if ~isempty(S.MCC.mask_img)
        mask_img = spm_read_vols(spm_vol(S.MCC.mask_img))>0;
    else
        try
            mii=find(~cellfun(@isempty,{D(d).model(:).mask_img_file}));
            mask_img = spm_read_vols(spm_vol(D(d).model(mii(1)).mask_img_file));
        catch
            mii=find(~cellfun(@isempty,{D(d).model_comp(:).mask_img_file}));
            mask_img = spm_read_vols(spm_vol(D(d).model_comp(mii(1)).mask_img_file));
        end
    end
    V.dim = size(mask_img);

    % method
    switch S.MCC.method
        case 'FDR' % can apply to 2D chan*time arrays, or 3D topo-time images.
                % rows = data to correct over; cols = different corrections
            disp('FDR...')
            % loop
            if S.MCC.model.index
                for i = S.MCC.model.index
                    for c = S.MCC.model.contrast{i}
                        p=D(d).model(i).con(c).p .* mask_img;
                        F=D(d).model(i).con(c).F .* mask_img;
                        [~,fdr_p] = FDR(p(:),S.MCC.thresh);
                        psize=size(p);
                        switch length(psize)
                            case 2
                                FDRmask = p<=fdr_p;
                            case 3
                                FDRmask = reshape(p<=fdr_p,psize(1),psize(2));
                            case 4
                                FDRmask = reshape(p<=fdr_p,psize(1),psize(2),psize(3));
                        end
                        % save images
                        [FDRmask_img] = topotime_3D(FDRmask,S);
                        [FDRmask_p_img] = topotime_3D(FDRmask.*p,S);
                        [FDRmask_F_img] = topotime_3D(FDRmask.*F,S);

                        D(d).model(i).con(c).MCC.FDRmask_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_FDRmask.nii']);
                        V.fname = D(d).model(i).con(c).MCC.FDRmask_img_file;
                        spm_write_vol(V,FDRmask_img);
                        D(d).model(i).con(c).MCC.FDRmask_p_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_FDRmask_p.nii']);
                        V.fname = D(d).model(i).con(c).MCC.FDRmask_p_img_file;
                        spm_write_vol(V,FDRmask_p_img);
                        D(d).model(i).con(c).MCC.FDRmask_F_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_FDRmask_F.nii']);
                        V.fname = D(d).model(i).con(c).MCC.FDRmask_F_img_file;
                        spm_write_vol(V,FDRmask_F_img);
                    end
                end
            end
            for i = S.MCC.model_comp.index
                p=D(d).model_comp(i).p .* mask_img;
                LR=D(d).model_comp(i).LR .* mask_img;
                [~,fdr_p] = FDR(p(:),S.MCC.thresh);
                psize=size(p);
                switch length(psize)
                    case 2
                        FDRmask = p<=fdr_p;
                    case 3
                        FDRmask = reshape(p<=fdr_p,psize(1),psize(2));
                    case 4
                        FDRmask = reshape(p<=fdr_p,psize(1),psize(2),psize(3));
                end
                % save images
                [FDRmask_img] = topotime_3D(FDRmask,S);
                [FDRmask_p_img] = topotime_3D(FDRmask.*p,S);
                [FDRmask_LR_img] = topotime_3D(FDRmask.*LR,S);

                D(d).model_comp(i).MCC.FDRmask_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_modelcomp_' num2str(i) '_FDRmask.nii']);
                V.fname = D(d).model_comp(i).MCC.FDRmask_img_file;
                spm_write_vol(V,FDRmask_img);
                D(d).model_comp(i).MCC.FDRmask_p_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_modelcomp_' num2str(i) '_FDRmask_p.nii']);
                V.fname = D(d).model_comp(i).MCC.FDRmask_p_img_file;
                spm_write_vol(V,FDRmask_p_img);
                D(d).model_comp(i).MCC.FDRmask_LR_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_modelcomp_' num2str(i) '_FDRmask_LR.nii']);
                V.fname = D(d).model_comp(i).MCC.FDRmask_LR_img_file;
                spm_write_vol(V,FDRmask_LR_img);
            end

        case 'pTFCE' % below is based on SPM and pTCFE code
            % First, estimate of smoothness based on [residual] images
            % Inputs:
            % loop
            if S.MCC.model.index
                for i = S.MCC.model.index
                    for c = S.MCC.model.contrast{i}

                        % image names
                        D(d).model(i).con(c).MCC.pTFCEmask_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_pTFCE_mask.nii']);
                        D(d).model(i).con(c).MCC.pTFCE_p_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_pTFCE_p.nii']);
                        D(d).model(i).con(c).MCC.pTFCE_Z_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_pTFCE_Z.nii']);
                        D(d).model(i).con(c).MCC.pTFCEmask_p_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_pTFCEmask_p.nii']);
                        D(d).model(i).con(c).MCC.pTFCEmask_Z_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_pTFCEmask_Z.nii']);

                        % if already created, load images from disk
                        if exist(D(d).model(i).con(c).MCC.pTFCEmask_Z_img_file,'file') && ~S.MCC.overwrite
                            pTFCEmask_img=spm_read_vols(spm_vol(D(d).model(i).con(c).MCC.pTFCEmask_img_file));
                            pTFCE_Z_img=spm_read_vols(spm_vol(D(d).model(i).con(c).MCC.pTFCE_Z_img_file));
                            pTFCE_p_img=spm_read_vols(spm_vol(D(d).model(i).con(c).MCC.pTFCE_p_img_file));
                            pTFCEmask_Z_img=spm_read_vols(spm_vol(D(d).model(i).con(c).MCC.pTFCEmask_Z_img_file));
                            pTFCEmask_p_img=spm_read_vols(spm_vol(D(d).model(i).con(c).MCC.pTFCEmask_p_img_file));
                        else

                            if ~exist('resid_vol','var')
                                %load residuals)
                                D(d).model(i).resid_vol_file = strrep(D(d).model(i).resid_file,'.mat','_vol.mat');
                                if exist(D(d).model(i).resid_vol_file,'file')
                                    disp(['pTFCE: for model ' num2str(i) ', loading resid vol file'])
                                    load(D(d).model(i).resid_vol_file);
                                else
                                    disp(['pTFCE: for model ' num2str(i) ', loading resid file'])
                                    load(D(d).model(i).resid_file);
                                    if ~exist('resid_vol','var')
                                        disp(['pTFCE: for model ' num2str(i) ', creating resid vol'])
                                        resid_vol = topotime_3D(resid,S);
                                        disp(['pTFCE: for model ' num2str(i) ', saving resid vol'])
                                        D(d).model(i).resid_vol_file = strrep(D(d).model(i).resid_file,'.mat','_vol.mat');
                                        save(D(d).model(i).resid_vol_file,'resid_vol','-v7.3');
                                    else
                                        D(d).model(i).resid_vol_file = D(d).model(i).resid_file;
                                    end
                                    %delete(D(d).model(i).resid_file)
                                end
                            end

                            DF = D(d).model(i).con(c).DF;%- A 2-vector, [n df], the original n & dof of the linear model     
                            if S.MCC.smooth_est_all || (i==1 && c==1)
                                %Vol   = D(d).model(i).resid_vol; % spm_vol(char(D(d).model(i).resid_img_file(:).name));%- Filenames of mapped standardized residual images
                                if ~isempty(S.MCC.mask_img)
                                    VM    = S.MCC.mask_img;
                                else
                                    VM    = D(d).model(i).mask_img_file;%- Filename of mapped mask image
                                end
                                % Outputs:
                                %FWHM - estimated FWHM in all image directions
                                %VRpv - handle of Resels per Voxel image
                                %R    - vector of resel counts
                                disp('pTFCE: smoothness estimation')
                                [FWHM,VRpv,R] = spm_est_smoothness_cab(resid_vol,VM,[size(resid_vol,4) DF(2)]);
                                D(d).model(i).con(c).MCC.FWHM = FWHM;
                                D(d).model(i).con(c).MCC.VRpv = VRpv;
                                D(d).model(i).con(c).MCC.R = R;
                            else
                                D(d).model(i).con(c).MCC.FWHM = D(d).model(1).con(1).MCC.FWHM;
                                D(d).model(i).con(c).MCC.VRpv = D(d).model(1).con(1).MCC.VRpv;
                                D(d).model(i).con(c).MCC.R = D(d).model(1).con(1).MCC.R;
                            end

                            % convert it to Z-score (from T or F)
                            % imgZ: Z-score image to enhance
                            if isfield(D(d).model(i).con(c),'T')
                                T_img = spm_read_vols(spm_vol(D(d).model(i).con(c).T_img_file));
                                imgZ = img_t2z(T_img, DF, 1);
                            elseif isfield(D(d).model(i).con(c),'F')
                                F_img = spm_read_vols(spm_vol(D(d).model(i).con(c).F_img_file));
                                imgZ = img_f2z(F_img, DF(1), DF(2), 1);
                            else
                                error('Statistical parameter type not supported');
                            end

                            disp(['pTFCE: model ' num2str(i) 'contrast ' num2str(c)])

                            nv = sum(mask_img(:)); % Number of voxels in mask (e.g. SPM.xVol.S)
                            % R:  Resel count, e.g. SPM.xVol.R(4).
                            minR = max(1,D(d).model(i).con(c).MCC.R(4));
                            [pTFCE_Z, pTFCE_p] = pTFCE(imgZ,mask_img,minR,nv);
                            pTFCE_Z(mask_img==0)=NaN;
                            pTFCE_p(mask_img==0)=NaN;
                            pTFCE_Z_img=pTFCE_Z;
                            pTFCE_p_img=pTFCE_p;
                            psize=size(pTFCE_p);
                            switch length(psize)
            %                     case 2
            %                         pTFCEmask_img = pTFCE_p<=0.05;
                                case 2
                                    pTFCEmask_img = reshape(pTFCE_p<=0.05,psize(1),psize(2));
                                case 3
                                    pTFCEmask_img = reshape(pTFCE_p<=0.05,psize(1),psize(2),psize(3));
                            end
                            if ~isempty(S.MCC.mask_img_post)
                                temp = spm_read_vols(spm_vol(S.MCC.mask_img_post))>0;
                                pTFCEmask_img=pTFCEmask_img.*temp;
                            end

                            % save images
                            V.fname = D(d).model(i).con(c).MCC.pTFCEmask_img_file;
                            spm_write_vol(V,pTFCEmask_img);
                            V.fname = D(d).model(i).con(c).MCC.pTFCE_p_img_file;
                            spm_write_vol(V,pTFCE_p_img);
                            V.fname = D(d).model(i).con(c).MCC.pTFCE_Z_img_file;
                            spm_write_vol(V,pTFCE_Z_img);
                            V.fname = D(d).model(i).con(c).MCC.pTFCEmask_p_img_file;
                            spm_write_vol(V,pTFCE_p_img.*pTFCEmask_img);
                            V.fname = D(d).model(i).con(c).MCC.pTFCEmask_Z_img_file;
                            spm_write_vol(V,pTFCE_Z_img.*pTFCEmask_img);

                            % convert back to original statistics (T or F)
                            if isfield(D(d).model(i).con(c),'T')
                                pTFCE_T_img = tinv(1-pTFCE_p, DF);
                                D(d).model(i).con(c).MCC.pTFCE_T_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_pTFCE_T.nii']);
                                V.fname = D(d).model(i).con(c).MCC.pTFCE_T_img_file;
                                spm_write_vol(V,pTFCE_T_img);
                                D(d).model(i).con(c).MCC.pTFCEmask_T_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_pTFCEmask_T.nii']);
                                V.fname = D(d).model(i).con(c).MCC.pTFCEmask_T_img_file;
                                spm_write_vol(V,pTFCE_T_img.*pTFCEmask_img);
                            elseif isfield(D(d).model(i).con(c),'F')
                                pTFCE_F_img = finv(1-pTFCE_p, DF(1), DF(2));
                                D(d).model(i).con(c).MCC.pTFCE_F_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_pTFCE_F.nii']);
                                V.fname = D(d).model(i).con(c).MCC.pTFCE_F_img_file;
                                spm_write_vol(V,pTFCE_F_img);
                                D(d).model(i).con(c).MCC.pTFCEmask_F_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_pTFCEmask_F.nii']);
                                V.fname = D(d).model(i).con(c).MCC.pTFCEmask_F_img_file;
                                spm_write_vol(V,pTFCE_F_img.*pTFCEmask_img);
                            end
                        end
                    end
                    clear resid_vol
                end
            end

            if S.MCC.model_comp.index
                for i = S.MCC.model_comp.index

                    % image names
                    D(d).model_comp(i).MCC.pTFCEmask_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_modelcomp_' num2str(i) '_pTFCE_mask.nii']);
                    D(d).model_comp(i).MCC.pTFCE_p_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_modelcomp_' num2str(i) '_pTFCE_p.nii']);
                    D(d).model_comp(i).MCC.pTFCE_Z_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_modelcomp_' num2str(i) '_pTFCE_Z.nii']);
                    D(d).model_comp(i).MCC.pTFCEmask_p_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_modelcomp_' num2str(i) '_pTFCEmask_p.nii']);
                    D(d).model_comp(i).MCC.pTFCEmask_Z_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_modelcomp_' num2str(i) '_pTFCEmask_Z.nii']);
                    D(d).model_comp(i).MCC.pTFCEmask_LR_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_modelcomp_' num2str(i) '_pTFCEmask_LR.nii']);

                    % if already created, load images from disk
                    if exist(D(d).model_comp(i).MCC.pTFCEmask_LR_img_file,'file') && ~S.MCC.overwrite
                        pTFCEmask_img=spm_read_vols(spm_vol(D(d).model_comp(i).MCC.pTFCEmask_img_file));
                        pTFCE_Z_img=spm_read_vols(spm_vol(D(d).model_comp(i).MCC.pTFCE_Z_img_file));
                        pTFCE_p_img=spm_read_vols(spm_vol(D(d).model_comp(i).MCC.pTFCE_p_img_file));
                        pTFCEmask_Z_img=spm_read_vols(spm_vol(D(d).model_comp(i).MCC.pTFCEmask_Z_img_file));
                        pTFCEmask_p_img=spm_read_vols(spm_vol(D(d).model_comp(i).MCC.pTFCEmask_p_img_file));
                        pTFCEmask_LR_img=spm_read_vols(spm_vol(D(d).model_comp(i).MCC.pTFCEmask_LR_img_file));
                    else

                        % use smoothness estimates from the original models
                        R4=[];
                        for p = 1:2
                            for c = 1:length(D(d).model(D(d).model_comp(i).pair(p)).con)
                                try
                                  R4=[R4 D(d).model(D(d).model_comp(i).pair(p)).con(c).MCC.R(4)];
                                end
                            end
                        end
                        R4=mean(R4);

                        % convert it to Z-score
                        % X2 is a folded-over and stretched-out normal: https://www.statpower.net/Content/310/Lecture%20Notes/ChiSquareAndF.pdf
            %             imgX = D(d).model_comp(i).LR_img;
            %             imgZ = sqrt(abs(imgX));

                        % CONVERT TO Z-SCORE
                        p_img = spm_read_vols(spm_vol(D(d).model_comp(i).pval_img_file));
                        LR_img = spm_read_vols(spm_vol(D(d).model_comp(i).LR_img_file));
                        pimg=p_img;
                        pimg(pimg==0) = realmin;
                        pimg(pimg==1) = 1-eps(1);
                        imgZ = -norminv(pimg); 
                        % remove negative values
                        imgZ(imgZ<0)=0;

                        disp(['pTFCE: model comp ' num2str(i)])
                        nv = sum(mask_img(:)); % Number of voxels in mask (e.g. SPM.xVol.S)
                        % R4:  Resel count, e.g. SPM.xVol.R(4).
                        minR = max(1,R4);
                        [pTFCE_Z, pTFCE_p] = pTFCE(imgZ,mask_img,R4,nv);
                        pTFCE_Z(mask_img==0)=NaN;
                        pTFCE_p(mask_img==0)=NaN;
                        pTFCE_Z_img=pTFCE_Z;
                        pTFCE_p_img=pTFCE_p;
                        psize=size(pTFCE_p);
                        switch length(psize)
                            case 2
                                pTFCEmask_img = reshape(pTFCE_p<=0.05,psize(1),psize(2));
                            case 3
                                pTFCEmask_img = reshape(pTFCE_p<=0.05,psize(1),psize(2),psize(3));
                        end
                        if ~isempty(S.MCC.mask_img_post)
                            temp = spm_read_vols(spm_vol(S.MCC.mask_img_post))>0;
                            pTFCEmask_img=pTFCEmask_img.*temp;
                        end

                        % save images
                        V.fname = D(d).model_comp(i).MCC.pTFCEmask_img_file;
                        spm_write_vol(V,pTFCEmask_img);
                        V.fname = D(d).model_comp(i).MCC.pTFCE_p_img_file;
                        spm_write_vol(V,pTFCE_p_img);
                        V.fname = D(d).model_comp(i).MCC.pTFCE_Z_img_file;
                        spm_write_vol(V,pTFCE_Z_img);
                        V.fname = D(d).model_comp(i).MCC.pTFCEmask_p_img_file;
                        spm_write_vol(V,pTFCE_p_img.*pTFCEmask_img);
                        V.fname = D(d).model_comp(i).MCC.pTFCEmask_Z_img_file;
                        spm_write_vol(V,pTFCE_Z_img.*pTFCEmask_img);
                        V.fname = D(d).model_comp(i).MCC.pTFCEmask_LR_img_file;
                        spm_write_vol(V,LR_img.*pTFCEmask_img);
                    end
                end
            end

    end
end
save(fullfile(S.MCC.path.outputs, 'D.mat'),'D');

function set_paths(S)
for p = 1:size(S.path.code,1)
    if S.path.code{p,1}
        addpath(genpath(S.path.code{p,2}));
    else
        addpath(S.path.code{p,2});
    end
end

