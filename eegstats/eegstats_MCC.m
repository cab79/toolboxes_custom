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
    V=makeV;

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
    V.dim(1:length(size(mask_img))) = size(mask_img);
    
    if S.MCC.model.index
        for i = S.MCC.model.index
            if S.MCC.estimate && (S.MCC.use_pTFCE || strcmp(S.MCC.method,'FWE_RF'))
            % First, estimate of smoothness based on [residual] images
                D(d).model(i).resid_vol_file = strrep(D(d).model(i).resid_file,'.mat','_vol.mat');
                D(d).model(i).resid_vol_dir = strrep(D(d).model(i).resid_vol_file,'.mat','_dir');
                if ~exist(D(d).model(i).resid_vol_dir,'dir') || length(dir(fullfile(D(d).model(i).resid_vol_dir,'*.nii'))) ~= length(D.model.fixeddesign)
                    mkdir(D(d).model(i).resid_vol_dir)
                    if exist(D(d).model(i).resid_vol_file,'file')
                        disp(['MCC: for model ' num2str(i) ', loading resid vol file'])
                        load(D(d).model(i).resid_vol_file);
                    else
                        disp(['MCC: for model ' num2str(i) ', loading resid file'])
                        load(D(d).model(i).resid_file);
                        if ~exist('resid_vol','var') && ~isempty(S.img.file.coord)
                            disp(['MCC: for model ' num2str(i) ', creating resid vol'])
                            resid_vol = topotime_3D(resid,S);
                            disp(['MCC: for model ' num2str(i) ', saving resid vol'])
                            D(d).model(i).resid_vol_file = strrep(D(d).model(i).resid_file,'.mat','_vol.mat');
                            save(D(d).model(i).resid_vol_file,'resid_vol','-v7.3');
                        else
                            D(d).model(i).resid_vol_file = D(d).model(i).resid_file;
                            if ~exist('resid_vol','var')
                                resid_vol=resid;
                            end
                        end
                        %delete(D(d).model(i).resid_file)
                    end
                    
                    disp('MCC: standardise residuals')
                    ndim_resid = ndims(resid_vol);
                    resid_vol_mean=mean(resid_vol,ndim_resid);
                    resid_vol_std=std(resid_vol,[],ndim_resid);
                    resid_vol_z = bsxfun(@rdivide,bsxfun(@minus,resid_vol,resid_vol_mean),resid_vol_std);
                    szR = size(resid_vol_z,ndim_resid);
                    clear resid_vol
                
                    % save as multiple volumes
                    pcc=0;
                    for tr = 1:szR
                        pc = floor(tr*100/szR);
                        if pcc<pc
                            pcc=pc;
                            disp(['writing resid volumes: ' num2str(pc) '%'])
                        end
                        trn = sprintf( '%06d', tr );
                        V.fname = fullfile(D(d).model(i).resid_vol_dir,['resid' trn '.nii']);
                        if ndim_resid ==4
                            spm_write_vol(V,resid_vol_z(:,:,:,tr));
                        elseif ndim_resid ==3
                            spm_write_vol(V,resid_vol_z(:,:,tr));
                        end
                    end
                end
                clear resid_vol_z
                fnames = dir(fullfile(D(d).model(i).resid_vol_dir,'*.nii'));
                resid_vol_z = fullfile(D(d).model(i).resid_vol_dir,{fnames.name});
                
                %- Filename of mapped mask image
                if ~isempty(S.MCC.mask_img)
                    VM    = S.MCC.mask_img;
                else
                    VM    = D(d).model(i).mask_img_file;
                end
                
                % smoothness estimation
                disp('MCC: smoothness estimation')
                
%                 % run one or multiple estimations?
%                 mask_dims = 1:ndims(mask_img);
%                 rep_dim = mask_dims(~ismember(mask_dims,S.MCC.dim));
%                 if isempty(rep_dim)
                    [L,Ve,ssq,Dx,Ix,Iy,Iz] = resid_smoothness(resid_vol_z,VM);
%                 else
%                     nrep = size(mask_img,rep_dim);
%                     for nr = 1:nrep
%                         [L,Ve,ssq,Dx,Ix,Iy,Iz] = resid_smoothness(resid_vol_z,VM);
%                     end
%                 end
            end

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
                    p_img = spm_read_vols(spm_vol(D(d).model(i).con(c).p_img_file));
                    szR = length(resid_vol_z);
                    if S.MCC.use_pTFCE || strcmp(S.MCC.method,'FWE_RF')
                        [FWHM,VRpv,R] = smoothness_df_correct(L,Ve,ssq,Dx,Ix,Iy,Iz,VM,[szR DF(2)]);
                        D(d).model(i).con(c).MCC.FWHM = FWHM;
                        D(d).model(i).con(c).MCC.VRpv = VRpv;
                        D(d).model(i).con(c).MCC.R = R;
                    end
                    
                    
                    % convert it to Z-score (from T or F)
                    % imgZ: Z-score image to enhance
                    if isfield(D(d).model(i).con(c),'T')
                        T_img = spm_read_vols(spm_vol(D(d).model(i).con(c).T_img_file));
                        Z_img = img_t2z(T_img, DF(2), 1);
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
                        nv = sum(mask_img(:)); % Number of voxels in mask (e.g. SPM.xVol.S)
                        % R:  Resel count, e.g. SPM.xVol.R(4).
                        minR = max(1,R(4));
                        [pTFCE_Z, pTFCE_p] = pTFCE(Z_img,mask_img,minR,nv);
                        pTFCE_Z(mask_img==0)=NaN;
                        pTFCE_p(mask_img==0)=NaN;
%                         figure; scatter(1:numel(pTFCE_Z),sort(pTFCE_Z(:)),'r'); hold on; scatter(1:numel(Z_img),sort(Z_img(:)),'b'); ylabel('Z scores');
%                         figure; scatter(Z_img(:),pTFCE_Z(:)); xlabel('original Z');ylabel('pTFCE Z');
                        Z_img=pTFCE_Z;
                        p_img=pTFCE_p;
                    end
                    zsize=size(Z_img);

                    switch S.MCC.method
                        case 'FWE_RF'
                            % FWE corrected critical height threshold at specified significance level
                            u = spm_uc_RF(S.MCC.thresh,DF,'Z',R(4),1);
                            switch length(zsize)
                                case 2
                                    MCCmask_img = reshape(Z_img>=u,zsize(1),zsize(2));
                                case 3
                                    MCCmask_img = reshape(Z_img>=u,zsize(1),zsize(2),zsize(3));
                            end

                        case 'FDR'
                            [fdr_p] = FDR(p_img(:),S.MCC.thresh);
                            psize=size(p_img);
                            switch length(psize)
                                case 1
                                    MCCmask_img = p_img<=fdr_p;
                                case 2
                                    MCCmask_img = reshape(p_img<=fdr_p,psize(1),psize(2));
                                case 3
                                    MCCmask_img = reshape(p_img<=fdr_p,psize(1),psize(2),psize(3));
                            end  
                    end

                    if ~isempty(S.MCC.mask_img_post)
                        temp = spm_read_vols(spm_vol(S.MCC.mask_img_post))>0;
                        MCCmask_img=MCCmask_img.*temp;
                    end
                    mask_p_img = p_img.*MCCmask_img;
                    mask_Z_img = Z_img.*MCCmask_img;

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
                plotdatY(:,c) = mask_Z_img(:);
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
                LR_img_thresh = LR_img>=S.MCC.LR_thresh | LR_img<=(1/S.MCC.LR_thresh);
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
                mask_Z_img = Z_img.*MCCmask_img;
                mask_p_img = p_img.*MCCmask_img;
                mask_LR_img = LR_img.*MCCmask_img;

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

function V = makeV
V.fname = '';
V.dim = [1 1 1];
V.dt = [16,0];
V.pinfo = [1;0;352];
V.mat = [-4.25000000000000,0,0,68;0,5.37500000000000,0,-100;0,0,1,-201;0,0,0,1];
V.n = [1,1];
V.descrip = '';
% V.private