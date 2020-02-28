function [L,Ve,ssq,D,Ix,Iy,Iz] = resid_smoothness(V,VM)
% CAB: modified to allow V to be an image/matrix rather than filename

% Estimation of smoothness based on [residual] images
% FORMAT [FWHM,VRpv,R] = spm_est_smoothness(V,VM,[ndf])
%
% V     - Filenames of mapped standardized residual images
% VM    - Filename of mapped mask image
% ndf   - A 2-vector, [n df], the original n & dof of the linear model
%
% FWHM  - estimated FWHM in all image directions
% VRpv  - handle of Resels per Voxel image
% R     - vector of resel counts
%__________________________________________________________________________
%
% spm_est_smoothness returns a spatial smoothness estimator based on the
% variances of the normalized spatial derivatives as described in K.
% Worsley, (1996). Inputs are a mask image and a number of standardized
% residual images, or any set of mean zero, unit variance images. Output
% is a global estimate of the smoothness expressed as the FWHM of an
% equivalent Gaussian point spread function. An estimate of resels per
% voxels (see spm_spm) is written as an image file ('RPV.<ext>') to the
% current directory.
%
% To improve the accuracy of the smoothness estimation the error degrees
% of freedom can be supplied.  Since it is not assumed that all residual
% images are passed to this function, the full, original sample size n
% must be supplied as well.
%
% The mask image specifies voxels, used in smoothness estimation, by
% assigning them non-zero values. The dimensions, voxel sizes, orientation
% of all images must be the same. The dimensions of the images can be of
% dimensions 0, 1, 2 and 3.
%
% Note that 1-dim images (lines) must exist in the 1st dimension and
% 2-dim images (slices) in the first two dimensions. The estimated fwhm
% for any non-existing dimension is infinity.
%__________________________________________________________________________
%
% Refs:
%
% K.J. Worsley (1996). An unbiased estimator for the roughness of a
% multivariate Gaussian random field. Technical Report, Department of
% Mathematics and Statistics, McGill University
%
% S.J. Kiebel, J.B. Poline, K.J. Friston, A.P. Holmes, and K.J. Worsley.
% Robust Smoothness Estimation in Statistical Parametric Maps Using
% Standardized Residuals from the General Linear Model. NeuroImage,
% 10:756-766, 1999.
%
% S. Hayasaka, K. Phan, I. Liberzon, K.J. Worsley, T.E. Nichols (2004).
% Nonstationary cluster-size inference with random field and permutation
% methods. NeuroImage, 22:676-687, 2004.
%__________________________________________________________________________
% Copyright (C) 2002-2015 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel, Tom Nichols
% $Id: spm_est_smoothness.m 6894 2016-09-30 16:48:46Z spm $


%-Assign input arguments
%--------------------------------------------------------------------------
if nargin < 1
    V   = spm_select(Inf,'image','Select residual images',{},pwd,'^ResI.*.{3}$');
end
if nargin < 2
    VM  = spm_select(1,'image','Select mask image',{},pwd,'^mask\..{3}$');
end
if nargin < 3, ndf = [NaN NaN]; end
if numel(ndf) ~= 2
    error('ndf argument must be of length 2 ([n df]).')
end

%-Initialise
%--------------------------------------------------------------------------
if isnumeric(V)
    sz=size(V);
    Ve=sz(end);
  
else
    Ve=numel(V);
    V       = spm_vol(V);
   
end
VM      = spm_vol(VM);

%-Dimensionality of image
%--------------------------------------------------------------------------
D        = 3 - sum(VM.dim(1:3) == 1);
if D == 0
    FWHM = [Inf Inf Inf];
    R    = [0 0 0];
    return;
end

%-Find voxels within mask
%--------------------------------------------------------------------------
d          = spm_read_vols(VM);
d = d>0; % cab - in case of NaNs
[Ix,Iy,Iz] = ndgrid(1:VM.dim(1),1:VM.dim(2),1:VM.dim(3));
Ix = Ix(d~=0); Iy = Iy(d~=0); Iz = Iz(d~=0);

%-Compute covariance of derivatives
%--------------------------------------------------------------------------
str   = 'Spatial non-sphericity (over scans)';
fprintf('%-40s: %30s',str,'...estimating derivatives');                 %-#
spm_progress_bar('Init',100,'smoothness estimation','');

L     = zeros(size(Ix,1),D,D);
ssq   = zeros(size(Ix,1),1);
for i = 1:Ve
    %disp(['smoothness ' num2str(i) '/' num2str(Ve)])
    if isnumeric(V)
        d=V(:,:,:,i);
        [dx,dy,dz] = gradient(d); 
        d=d(Ix);dx=dx(Ix);dy=dy(Ix);dz=dz(Ix);
        dx(isnan(dx))=0;dy(isnan(dy))=0;dz(isnan(dz))=0;
    else
        [d,dx,dy,dz] = spm_sample_vol(V(i),Ix,Iy,Iz,1);
    end
    
    % sum of squares
    %----------------------------------------------------------------------
    ssq  = ssq + d.^2;
    
    % covariance of finite differences
    %----------------------------------------------------------------------
    if D >= 1
        L(:,1,1) = L(:,1,1) + dx.*dx;
    end
    if D >= 2
        L(:,1,2) = L(:,1,2) + dx.*dy;
        L(:,2,2) = L(:,2,2) + dy.*dy;
    end
    if D >= 3
        L(:,1,3) = L(:,1,3) + dx.*dz;
        L(:,2,3) = L(:,2,3) + dy.*dz;
        L(:,3,3) = L(:,3,3) + dz.*dz;
    end
    
    
end
