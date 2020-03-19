function D=eegstats_MCC_old(S,varargin)
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
    switch S.MCC.use_pTFCE
        case 0 % can apply to 2D chan*time arrays, or 3D topo-time images.
                % rows = data to correct over; cols = different corrections
            disp('FDR...')
            % loop
            if S.MCC.model.index
                for i = S.MCC.model.index
                    for c = S.MCC.model.contrast{i}
                        p_img = spm_read_vols(spm_vol(D(d).model(i).con(c).p_img_file)) .* mask_img;
                        F_img=spm_read_vols(spm_vol(D(d).model(i).con(c).F_img_file)) .* mask_img;
                        [~,fdr_p] = FDR(p_img(:),S.MCC.thresh);
                        psize=size(p_img);
                        switch length(psize)
                            case 1
                                FDRmask_img = p_img<=fdr_p;
                            case 2
                                FDRmask_img = reshape(p_img<=fdr_p,psize(1),psize(2));
                            case 3
                                FDRmask_img = reshape(p_img<=fdr_p,psize(1),psize(2),psize(3));
                        end
                        % save images
%                         [FDRmask_img] = topotime_3D(FDRmask,S);
                        [FDRmask_p_img] = FDRmask_img.*p_img;
                        [FDRmask_F_img] = FDRmask_img.*F_img;
                        DF = D(d).model(i).con(c).DF;%- A 2-vector, [n df], the original n & dof of the linear model     
                        FDRmask_Z_img = img_f2z(FDRmask_F_img, DF(1), DF(2), 1);

                        D(d).model(i).con(c).MCC.FDRmask_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_FDRmask.nii']);
                        V.fname = D(d).model(i).con(c).MCC.FDRmask_img_file;
                        spm_write_vol(V,FDRmask_img);
                        D(d).model(i).con(c).MCC.FDRmask_p_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_FDRmask_p.nii']);
                        V.fname = D(d).model(i).con(c).MCC.FDRmask_p_img_file;
                        spm_write_vol(V,FDRmask_p_img);
                        D(d).model(i).con(c).MCC.FDRmask_F_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_FDRmask_F.nii']);
                        V.fname = D(d).model(i).con(c).MCC.FDRmask_F_img_file;
                        spm_write_vol(V,FDRmask_F_img);
                        D(d).model(i).con(c).MCC.FDRmask_Z_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_FDRmask_Z.nii']);
                        V.fname = D(d).model(i).con(c).MCC.FDRmask_Z_img_file;
                        spm_write_vol(V,FDRmask_Z_img)
                    end
                end
            end
            if S.MCC.model_comp.index
                for i = S.MCC.model_comp.index
                    p_img=spm_read_vols(spm_vol(D(d).model_comp(i).pval_img_file)) .* mask_img;
                    LR_img=spm_read_vols(spm_vol(D(d).model_comp(i).LR_img_file)) .* mask_img;
                    [~,fdr_p] = FDR(p_img(:),S.MCC.thresh);
                    psize=size(p_img);
                    switch length(psize)
                        case 1
                            FDRmask_img = p_img<=fdr_p;
                        case 2
                            FDRmask_img = reshape(p_img<=fdr_p,psize(1),psize(2));
                        case 3
                            FDRmask_img = reshape(p_img<=fdr_p,psize(1),psize(2),psize(3));
                    end
                    % save images
    %                 [FDRmask_img] = topotime_3D(FDRmask,S);
                    [FDRmask_p_img] = FDRmask_img.*p_img;
                    [FDRmask_LR_img] = FDRmask_img.*LR_img;

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
            end

        case 1 % below is based on SPM and pTCFE code
            if S.MCC.model.index
                for i = S.MCC.model.index
                    % First, estimate of smoothness based on [residual] images
                    if S.MCC.estimate
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
                        disp('MCC: smoothness estimation')
                        % first standardise the residuals
                        resid_vol_mean=nanmean(resid_vol,4);
                        resid_vol_std=nanstd(resid_vol,[],4);
                        resid_vol_z = bsxfun(@rdivide,bsxfun(@minus,resid_vol,resid_vol_mean),resid_vol_std);
                        clear resid_vol
                        [L,Ve,ssq,Dx,Ix,Iy,Iz] = resid_smoothness(resid_vol_z,VM);
                        szR = size(resid_vol_z,4);
                        clear resid_vol_z
                    end

                    for c = S.MCC.model.contrast{i}
                        disp(['MCC: model ' num2str(i) 'contrast ' num2str(c)])

                        % image names
                        D(d).model(i).con(c).MCC.pTFCEmask_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_pTFCE_mask.nii']);
                        D(d).model(i).con(c).MCC.pTFCE_p_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_pTFCE_p.nii']);
                        D(d).model(i).con(c).MCC.pTFCE_Z_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_pTFCE_Z.nii']);
                        D(d).model(i).con(c).MCC.pTFCEmask_p_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_pTFCEmask_p.nii']);
                        D(d).model(i).con(c).MCC.pTFCEmask_Z_img_file = fullfile(S.MCC.path.outputs, [save_pref num2str(d) '_model_' num2str(i) '_con_' num2str(c) '_pTFCEmask_Z.nii']);

                        % if already created, load images from disk
                        if exist(D(d).model(i).con(c).MCC.pTFCEmask_Z_img_file,'file') && ~S.MCC.estimate
                            pTFCEmask_img=spm_read_vols(spm_vol(D(d).model(i).con(c).MCC.pTFCEmask_img_file));
                            pTFCE_Z_img=spm_read_vols(spm_vol(D(d).model(i).con(c).MCC.pTFCE_Z_img_file));
                            pTFCE_p_img=spm_read_vols(spm_vol(D(d).model(i).con(c).MCC.pTFCE_p_img_file));
                            pTFCEmask_Z_img=spm_read_vols(spm_vol(D(d).model(i).con(c).MCC.pTFCEmask_Z_img_file));
                            pTFCEmask_p_img=spm_read_vols(spm_vol(D(d).model(i).con(c).MCC.pTFCEmask_p_img_file));
                        elseif ~S.MCC.estimate
                            error('images do not exist - estimation required')
                        else
                            
                            DF = D(d).model(i).con(c).DF;   
%                             DF(2)=szR;
%                             if S.MCC.smooth_est_all || (i==1 && c==1)
                                
                                [FWHM,VRpv,R] = smoothness_df_correct(L,Ve,ssq,Dx,Ix,Iy,Iz,VM,[szR DF(2)]);
                                D(d).model(i).con(c).MCC.FWHM = FWHM;
                                D(d).model(i).con(c).MCC.VRpv = VRpv;
                                D(d).model(i).con(c).MCC.R = R;
%                             else
%                                 D(d).model(i).con(c).MCC.FWHM = D(d).model(1).con(1).MCC.FWHM;
%                                 D(d).model(i).con(c).MCC.VRpv = D(d).model(1).con(1).MCC.VRpv;
%                                 D(d).model(i).con(c).MCC.R = D(d).model(1).con(1).MCC.R;
%                             end

                            % convert it to Z-score (from T or F)
                            % imgZ: Z-score image to enhance
                            if isfield(D(d).model(i).con(c),'T')
                                T_img = spm_read_vols(spm_vol(D(d).model(i).con(c).T_img_file));
                                imgZ = img_t2z(T_img, DF(2), 1);
                            elseif isfield(D(d).model(i).con(c),'F')
                                F_img = spm_read_vols(spm_vol(D(d).model(i).con(c).F_img_file));
                                if any(F_img(:)<0)
                                    F_img(F_img(:)<0)=0;
                                end
                                imgZ = img_f2z(F_img, DF(1), DF(2), 1);
                            else
                                error('Statistical parameter type not supported');
                            end


                            nv = sum(mask_img(:)); % Number of voxels in mask (e.g. SPM.xVol.S)
                            % R:  Resel count, e.g. SPM.xVol.R(4).
                            minR = max(1,D(d).model(i).con(c).MCC.R(4));
                            [pTFCE_Z, pTFCE_p] = pTFCE(imgZ,mask_img,minR,nv);
                            pTFCE_Z(mask_img==0)=NaN;
                            pTFCE_p(mask_img==0)=NaN;
                            pTFCE_Z_img=pTFCE_Z;
                            pTFCE_p_img=pTFCE_p;
                            zsize=size(pTFCE_Z);
                            
                            switch S.MCC.method
                                case 'FWE'
                                    % FWE corrected critical height threshold at specified significance level
                                    u = spm_uc_RF(S.MCC.thresh,DF,'Z',R(4),1);
                                    switch length(zsize)
                    %                     case 2
                    %                         pTFCEmask_img = pTFCE_p<=0.05;
                                        case 2
                                            pTFCEmask_img = reshape(pTFCE_Z_img>=u,zsize(1),zsize(2));
                                        case 3
                                            pTFCEmask_img = reshape(pTFCE_Z_img>=u,zsize(1),zsize(2),zsize(3));
                                    end
                                
                                case 'FDR'
                                    [fdr_p] = FDR(pTFCE_p(:),S.MCC.thresh);
                                    psize=size(pTFCE_p);
                                    switch length(psize)
                                        case 1
                                            pTFCEmask_img = pTFCE_p<=fdr_p;
                                        case 2
                                            pTFCEmask_img = reshape(pTFCE_p<=fdr_p,psize(1),psize(2));
                                        case 3
                                            pTFCEmask_img = reshape(pTFCE_p<=fdr_p,psize(1),psize(2),psize(3));
                                    end  
                            end
                            
                            
                            if ~isempty(S.MCC.mask_img_post)
                                temp = spm_read_vols(spm_vol(S.MCC.mask_img_post))>0;
                                pTFCEmask_img=pTFCEmask_img.*temp;
                            end
                            pTFCEmask_p_img = pTFCE_p_img.*pTFCEmask_img;
                            pTFCEmask_Z_img = pTFCE_Z_img.*pTFCEmask_img;

                            % save images
                            V.fname = D(d).model(i).con(c).MCC.pTFCEmask_img_file;
                            spm_write_vol(V,pTFCEmask_img);
                            V.fname = D(d).model(i).con(c).MCC.pTFCE_p_img_file;
                            spm_write_vol(V,pTFCE_p_img);
                            V.fname = D(d).model(i).con(c).MCC.pTFCE_Z_img_file;
                            spm_write_vol(V,pTFCE_Z_img);
                            V.fname = D(d).model(i).con(c).MCC.pTFCEmask_p_img_file;
                            spm_write_vol(V,pTFCEmask_p_img);
                            V.fname = D(d).model(i).con(c).MCC.pTFCEmask_Z_img_file;
                            spm_write_vol(V,pTFCEmask_Z_img);

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
                    
                        % plot data
                        plotdatY(:,c) = pTFCEmask_Z_img(:);
                        plotdatX(:,c) = c*double(plotdatY(:,c)>0);
                        xtl{c} = D(d).model(i).con(c).term;
                    end
                    
                    % plot
                    figname = ['model ' num2str(i) ' post-threshold Z values'];
                    xlab = 'contrast';
                    ylab = 'Z values';
                    try
                        jitterplot(figname,plotdatX,plotdatY,xlab,ylab,xtl)
                    catch
                        disp(['no significant effects for model ' num2str(i)])
                    end
                    clear plotdatX plotdatY xtl
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
                    if exist(D(d).model_comp(i).MCC.pTFCEmask_LR_img_file,'file') && ~S.MCC.estimate
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
                            LR_DF(p)=length(D.model(D(d).model_comp(i).pair(p)).coeff);
                        end
                        R4=mean(R4);
                        LR_DF = LR_DF(2)-LR_DF(1);
                        
                        % mask
                        contrast_Z_img=[];
                        for p = 1:2
                            for c = 1:length(D(d).model(D(d).model_comp(i).pair(p)).con)
                                contrast_Z_img(:,:,:,p,c)=spm_read_vols(spm_vol(D(d).model(D(d).model_comp(i).pair(p)).con(c).MCC.pTFCEmask_Z_img_file));
                            end
                        end
                        contrast_mask = any(contrast_Z_img,[4 5]);
                        mask_img = mask_img.*contrast_mask;

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

                        nv = sum(mask_img(:)); % Number of voxels in mask (e.g. SPM.xVol.S)
                        % R4:  Resel count, e.g. SPM.xVol.R(4).
                        minR = max(1,R4);
                        
                        disp(['pTFCE: model comp ' num2str(i)])
                        [pTFCE_Z, pTFCE_p] = pTFCE(imgZ,mask_img,R4,nv);
                        pTFCE_Z(mask_img==0)=NaN;
                        pTFCE_p(mask_img==0)=NaN;
                        pTFCE_Z_img=pTFCE_Z;
                        pTFCE_p_img=pTFCE_p;
                        zsize=size(pTFCE_Z);
%                         figure;plot(pTFCE_Z(:))
                        
%                         if 1 % use LR threshold
                            LR_img_thresh = LR_img>=S.MCC.LR_thresh | LR_img<=(1/S.MCC.LR_thresh);
                            switch length(zsize)
                                case 2
                                    pTFCEmask_img = reshape(LR_img_thresh,zsize(1),zsize(2));
                                case 3
                                    pTFCEmask_img = reshape(LR_img_thresh,zsize(1),zsize(2),zsize(3));
                            end
                          %% THESE METHODS WILL NOT IDENTIFY EVIDENCE IN FAVOUR OF FIRST MODEL (LR<1)
%                         elseif 0 % correction using LR_stat and X2 distribution (DOESN'T USE TFCE)
%                             % FWE corrected critical height threshold at specified significance level
%                             u = spm_uc_RF(S.MCC.thresh,[1 LR_DF],'X',R4,1);
% 
%                             switch length(zsize)
%                                 case 2
%                                     pTFCEmask_img = reshape(LR_img>=u,zsize(1),zsize(2));
%                                 case 3
%                                     pTFCEmask_img = reshape(LR_img>=u,zsize(1),zsize(2),zsize(3));
%                             end
%                         elseif 0 % correction using Z image
%                             % FWE corrected critical height threshold at specified significance level
%                             u = spm_uc_RF(S.MCC.thresh,[NaN NaN],'Z',R4,1);
% 
%                             switch length(zsize)
%                                 case 2
%                                     pTFCEmask_img = reshape(pTFCE_Z_img>=u,zsize(1),zsize(2));
%                                 case 3
%                                     pTFCEmask_img = reshape(pTFCE_Z_img>=u,zsize(1),zsize(2),zsize(3));
%                             end
%                         elseif 0 % FDR
%                             [fdr_p] = FDR(pTFCE_p(:),S.MCC.thresh);
%                             psize=size(pTFCE_p);
%                             switch length(psize)
%                                 case 1
%                                     pTFCEmask_img = pTFCE_p<=fdr_p;
%                                 case 2
%                                     pTFCEmask_img = reshape(pTFCE_p<=fdr_p,psize(1),psize(2));
%                                 case 3
%                                     pTFCEmask_img = reshape(pTFCE_p<=fdr_p,psize(1),psize(2),psize(3));
%                             end  
%                         end
                        if ~isempty(S.MCC.mask_img_post)
                            temp = spm_read_vols(spm_vol(S.MCC.mask_img_post))>0;
                            pTFCEmask_img=pTFCEmask_img.*temp;
                        end
                        pTFCEmask_Z_img = pTFCE_Z_img.*pTFCEmask_img;
                        pTFCEmask_p_img = pTFCE_p_img.*pTFCEmask_img;
                        pTFCEmask_LR_img = LR_img.*pTFCEmask_img;

                        % save images
                        V.fname = D(d).model_comp(i).MCC.pTFCEmask_img_file;
                        spm_write_vol(V,pTFCEmask_img);
                        V.fname = D(d).model_comp(i).MCC.pTFCE_p_img_file;
                        spm_write_vol(V,pTFCE_p_img);
                        V.fname = D(d).model_comp(i).MCC.pTFCE_Z_img_file;
                        spm_write_vol(V,pTFCE_Z_img);
                        V.fname = D(d).model_comp(i).MCC.pTFCEmask_p_img_file;
                        spm_write_vol(V,pTFCEmask_p_img);
                        V.fname = D(d).model_comp(i).MCC.pTFCEmask_Z_img_file;
                        spm_write_vol(V,pTFCEmask_Z_img);
                        V.fname = D(d).model_comp(i).MCC.pTFCEmask_LR_img_file;
                        spm_write_vol(V,pTFCEmask_LR_img);
                    end
                    % plot data
                    plotdatY(:,i) = pTFCEmask_LR_img(:);
                    plotdatX(:,i) = i*double(plotdatY(:,i)>0);
                    xtl{i} = num2str(D(d).model_comp(i).pair);
                end
                figname = 'model comparison post-threshold LR values';
                xlab = 'models compared';
                ylab = 'likelihood ratios';
                try
                    jitterplot(figname,plotdatX,plotdatY,xlab,ylab,xtl)
                catch
                    disp(['no significant effects for model comp ' num2str(i)])
                end
                ylab = 'log likelihood ratios';
                try
                    jitterplot(figname,plotdatX,log(plotdatY),xlab,ylab,xtl)
                catch
                    disp(['no significant effects for model comp ' num2str(i)])
                end
                clear plotdatX plotdatY xtl
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

function jitterplot(figname,plotdatX,plotdatY,xlab,ylab,xtl)
plotdatX=-plotdatX;
plotdatX(plotdatX==0)=NaN;
figure('name',figname)
minmax = -[size(plotdatX,2)+0.5,0.5];
hold on
s=scatter(plotdatX(:), plotdatY(:),10,'k','filled', 'jitter','on', 'jitterAmount',0.3);
s.MarkerFaceAlpha=0.5;
if strcmp(ylab,'Z values')
    line(minmax,[1.645 1.645],'Color','red','LineStyle','--') % threshold for one-tailed z score, uncorrected
end
hold off
xlim(minmax)
xticks(-size(plotdatX,2):-1)
xticklabels(fliplr(xtl))
xlabel(xlab)
ylabel(ylab)
view([90 -90])
