function D=eegstats_threshimages(S,varargin)
% Apply thresholds to any statistical image (beta, F, T, etc.) and save output images.
% Options: F images from con; beta images from coeff; logl/waic/r2 from model.

dbstop if error

%% find D saved from createimages or MCC functions
if isempty(varargin)
   % try loading the data if not supplied in varargin
   try
       load(fullfile(S.thresh.path.inputs,'D.mat'),'D') 
   catch
       error('Supply a D structure or path/file of the data')
   end
else
    D= varargin{1};
end

% set paths
S.path.code = {
    };
set_paths(S)

if length(D)>1
    save_pref = 'subject_';
else
    save_pref = 'group_';
end

%% Loop
for d=1:length(D)
    
    % get model and contrast indices
    if isempty(S.thresh.model.index)
        S.thresh.model.index = 1:length(D(d).model);
    end
    if isempty(S.thresh.model.coeff) && isfield(D(d).model(1),'coeff')
        for i = S.thresh.model.index
            S.thresh.model.coeff{i} = 1:length(D(d).model(i).coeff);
        end
    end
    if isempty(S.thresh.model.contrast) && isfield(D(d).model(1),'con')
        for i = S.thresh.model.index
            S.thresh.model.contrast{i} = 1:length(D(d).model(i).con);
        end
    end
    if isempty(S.thresh.model_comp.index)
        S.thresh.model_comp.index = 1:length(D(d).model_comp);
    end

    if S.thresh.model.index
        for i = S.thresh.model.index
            disp(['extracting clusters for ' save_pref num2str(d) ', model ' num2str(i)])
            
            % Model quantities
            if ismember(S.thresh.image,{'logl','waic','r2'})
                try
                    img_file = [S.thresh.image '_file'];
                catch
                    img_file = [S.thresh.image '_img_file'];
                end
                img_file = strrep(img_file,'.nii','');
                base_fname=D(d).model(i).(img_file);
                img = spm_read_vols(spm_vol(base_fname));

                % calculate images with a range of thresholds
                threshimgs = calc_thresh_images(S,img);

                % save thresholded images
                if S.thresh.save_thresh_images
                    save_thresh_images(threshimgs,base_fname)
                end

                % update D with file name
                threshimgs_file = strrep(img_file,'_file',['_thresh' num2str(ti) '_file']);
                D(d).model(i).(threshimgs_file) = img_file;
            end
            
            % Coefficients
            if ismember(S.thresh.image,{'b'})
                for c = S.thresh.model.coeff{i}
                    try
                        img_file = [S.thresh.image '_file'];
                    catch
                        img_file = [S.thresh.image '_img_file'];
                    end
                    img_file = strrep(img_file,'.nii','');
                    base_fname=D(d).model(i).coeff(c).(img_file);
                    img = spm_read_vols(spm_vol(base_fname));

                    % calculate images with a range of thresholds
                    threshimgs = calc_thresh_images(S,img);

                    % save thresholded images
                    if S.thresh.save_thresh_images
                        save_thresh_images(threshimgs,base_fname)
                    end

                    % update D with file name
                    threshimgs_file = strrep(img_file,'_file',['_thresh' num2str(ti) '_file']);
                    D(d).model(i).coeff(c).(threshimgs_file) = img_file;

                end
            end
            
            % Contrasts
            if ismember(S.thresh.image,{'F','T','p'})
                for c = S.thresh.model.contrast{i}
                    try
                        img_file = [S.thresh.image '_file'];
                        img_file = strrep(img_file,'.nii','');
                        try
                            base_fname=D(d).model(i).con(c).(img_file);
                            isMCC=0;
                        catch
                            base_fname=D(d).model(i).con(c).MCC.(img_file);
                            isMCC=1;
                        end
                    catch
                        img_file = [S.thresh.image '_img_file'];
                        img_file = strrep(img_file,'.nii','');
                        try
                            base_fname=D(d).model(i).con(c).(img_file);
                            isMCC=0;
                        catch
                            base_fname=D(d).model(i).con(c).MCC.(img_file);
                            isMCC=1;
                        end
                    end
                    img = spm_read_vols(spm_vol(base_fname));

                    % calculate images with a range of thresholds
                    threshimgs = calc_thresh_images(S,img);

                    % save thresholded images
                    if S.thresh.save_thresh_images
                        save_thresh_images(threshimgs,base_fname)
                    end

                    % update D with file name
                    threshimgs_file = strrep(img_file,'_file',['_thresh' num2str(ti) '_file']);
                    if isMCC
                        D(d).model(i).con(c).MCC.(threshimgs_file) = img_file;
                    else
                        D(d).model(i).con(c).(threshimgs_file) = img_file;
                    end

                end
            end
            
        end
        for i = S.thresh.model_comp.index
            disp(['extracting clusters for ' save_pref num2str(d) ', model comparison ' num2str(i)])
            % Model comparison
            if ismember(S.thresh.image,{'LR','LRStat'})
                try
                    img_file = [S.thresh.image '_file'];
                catch
                    img_file = [S.thresh.image '_img_file'];
                end
                img_file = strrep(img_file,'.nii','');
                base_fname=D(d).model_comp(i).(img_file);
                img = spm_read_vols(spm_vol(base_fname));

                % calculate images with a range of thresholds
                threshimgs = calc_thresh_images(S,img);

                % save thresholded images
                if S.thresh.save_thresh_images
                    save_thresh_images(threshimgs,base_fname)
                end

                % update D with file name
                threshimgs_file = strrep(img_file,'_file',['_thresh' num2str(ti) '_file']);
                D(d).model_comp(i).(threshimgs_file) = img_file;
            end
        end
    end
end
save(fullfile(S.thresh.path.outputs, 'D.mat'),'D');


function set_paths(S)
for p = 1:length(S.path.code)
    if S.path.code{p,1}
        addpath(genpath(S.path.code{p,2}));
    else
        addpath(S.path.code{p,2});
    end
end

function threshimgs = calc_thresh_images(S,img)
% calculate images with a one or a range of thresholds
switch S.thresh.Thresh_method
    case 'range'
        Nthresh=S.thresh.Thresh;
        incr=max(img(:))/Nthresh; % increment
        vrange = 0:incr:max(img(:))-incr; % value range

    case 'values'
        Nthresh= length(S.thresh.Thresh);
        vrange = S.thresh.Thresh;
end
if any(vrange)
    for vr = 1:Nthresh
        threshimgs{vr} = img;
        threshimgs{vr}(threshimgs{vr}<vrange(vr))=0;
        threshimgs{vr}(isnan(threshimgs{vr}))=0;
    end
end

function save_thresh_images(threshimgs,base_fname)
% save thresholded images
V=spm_vol(base_fname);
for ti = 1:length(threshimgs)
    threshimgs_fname = strrep(base_fname,'.nii',['_thresh' num2str(ti) '.nii']);
    V.fname = threshimgs_fname;
    spm_write_vol(V,threshimgs{ti})
end