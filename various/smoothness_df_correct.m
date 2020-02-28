function [FWHM,VRpv,R] = smoothness_df_correct(L,Ve,ssq,D,Ix,Iy,Iz,VM,ndf)

VM      = spm_vol(VM);
n_full  = ndf(1);
edf     = ndf(2);

%-Scale sum into an average (and account for DF)
% The standard result uses normalised residuals e/sqrt(RSS) and
%
%  \hat\Lambda = grad(e/sqrt(RSS))' * grad(e/sqrt(RSS))
%
% In terms of standardized residuals e/sqrt(RMS) this is
%
%  \hat\Lambda = (1/DF) * grad(e/sqrt(RMS))' * grad(e/sqrt(RMS))
%
% but both of these expressions assume that the RSS or RMS correspond to
% the full set of residuals considered.  However, spm_spm only considers
% upto MAXRES residual images.  To adapt, re-write the above as a scaled
% average over n scans
%
%  \hat\Lambda = (n/DF) * ( (1/n) * grad(e/sqrt(RMS))' * grad(e/sqrt(RMS)) )
%
% I.e. the roughness estimate \hat\Lambda is an average of outer products
% of standardized residuals (where the average is over scans), scaled by
% n/DF. Hence, we can use only a subset of scans simply by replacing this 
% last average term with an average over the subset.
%
% See Hayasaka et al, p. 678, for more on estimating roughness with
% standardized residuals (e/sqrt(RMS)) instead of normalised residuals
% (e/sqrt(RSS)). Note that the names arise from the fact that
% sqrt(RSS) = sqrt(r'*r) is norm(r), while sqrt(RMS) = sqrt(r'*r/edf)
% is the unbiased (ReML) estimate of the standard deviation.
%--------------------------------------------------------------------------
L  = L/Ve;      % Average
L  = L*(n_full/edf);  % Scale

%-Initialise RESELS per voxel image
%--------------------------------------------------------------------------
VRpv  = struct( 'fname',['RPV' spm_file_ext],...
    'dim',      VM.dim(1:3),...
    'dt',       [spm_type('float64') spm_platform('bigend')],...
    'mat',      VM.mat,...
    'pinfo',    [1 0 0]',...
    'descrip',  'spm_spm: resels per voxel');
VRpv   = spm_create_vol(VRpv);

%-Evaluate determinant (and xyz components for FWHM)
%--------------------------------------------------------------------------
if D == 1
    resel_xyz = L;
    resel_img = L;
end
if D == 2
    resel_xyz = [L(:,1,1) L(:,2,2)];
    resel_img = L(:,1,1).*L(:,2,2) - ...
                L(:,1,2).*L(:,1,2);
end
if D == 3
    resel_xyz = [L(:,1,1) L(:,2,2)  L(:,3,3)];
    resel_img = L(:,1,1).*L(:,2,2).*L(:,3,3) + ...
                L(:,1,2).*L(:,2,3).*L(:,1,3)*2 - ...
                L(:,1,1).*L(:,2,3).*L(:,2,3) - ...
                L(:,1,2).*L(:,1,2).*L(:,3,3) - ...
                L(:,1,3).*L(:,2,2).*L(:,1,3);
end    
resel_img(resel_img<0) = 0;
% Convert det(Lambda) and diag(Lambda) to units of resels
resel_img = sqrt(resel_img/(4*log(2))^D);
resel_xyz = sqrt(resel_xyz/(4*log(2)));


%-Optional mask-weighted smoothing of RPV image
%--------------------------------------------------------------------------
if spm_get_defaults('stats.rft.nonstat')
    fwhm_vox = 3;
else
    fwhm_vox = 0;
end
if any(fwhm_vox)
    if length(fwhm_vox) == 1, fwhm_vox = fwhm_vox([1 1 1]);  end
    
    % Convert resel_img at in-mask voxels to resel volume
    mask = spm_read_vols(VM) > 0;
    RPV = zeros(size(mask));
    RPV(mask) = resel_img;
    
    % Remove invalid mask voxels, i.e. edge voxels with missing derivatives
    smask = double(mask & isfinite(RPV)); % leaves mask for resel_img below
    
    % Smooth RPV volume (note that NaNs are treated as zeros in spm_smooth)
    spm_smooth(RPV, RPV, fwhm_vox);
    
    % Smooth mask and decide how far to trust smoothing-based extrapolation
    spm_smooth(smask, smask, fwhm_vox);
    infer = smask > 1e-3; % require sum of voxel's in-mask weights > 1e-3
    
    % Normalise smoothed RPV by smoothed mask
    RPV( infer) = RPV(infer) ./ smask(infer);
    RPV(~infer) = NaN; % spm_list handles remaining (unlikely) in-mask NaNs
    
    % Get data at in-mask voxels; smoothed resel_img conforms with original
    resel_img = RPV(mask);
end

%-Write Resels Per Voxel image
%--------------------------------------------------------------------------

for i = 1:VM.dim(3)
    d = NaN(VM.dim(1:2));
    I = find(Iz == i);
    if ~isempty(I)
        d(sub2ind(VM.dim(1:2), Ix(I), Iy(I))) = resel_img(I);
    end
    VRpv = spm_write_plane(VRpv, d, i);
end


%-(unbiased) RESEL estimator and Global equivalent FWHM
% where we desire FWHM with components proportional to 1./mean(resel_xyz),
% but scaled so prod(1./FWHM) agrees with (the unbiased) mean(resel_img).
%--------------------------------------------------------------------------
i     = isnan(ssq) | ssq < sqrt(eps);
resel_img = mean(resel_img(~i,:));
resel_xyz = mean(resel_xyz(~i,:));

RESEL = resel_img^(1/D)*(resel_xyz/prod(resel_xyz)^(1/D));
FWHM  = full(sparse(1,1:D,1./RESEL,1,3));
FWHM(isnan(FWHM)) = 0;
FWHM(~FWHM) = 1;

%-resel counts for search volume (defined by mask)
%--------------------------------------------------------------------------
% R0   = spm_resels_vol(VM,[1 1 1])';
% R    = R0.*(resel.^([0:3]/3));
% OR
R      = spm_resels_vol(VM,FWHM)';


