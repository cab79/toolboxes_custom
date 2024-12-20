
% create 3D topo-over-time image and mask from chan x time 2D array
% code from spm_eeg_convert2images
% Y = chan x time x trials/contrasts (3D), chan x time (2D), or chan (1D) in rows.
% Or, Y can be a structure of directory/filenames, in which case output is
% also filenames after saving them to disk.
% S is a structure than contains settings, e.g. 
    % S.img.file.coord = 'path/filename'; % path and filename of channel coordinates
    % S.img.file.coord_struct = 'c2d'; % name of structure within file that has coordinates
    % S.img.chansel = [120,1,125,121,122,123,124,117,118, 111, 110, 112, 106,31,7,129, 55, 80, 79, 86, 92, 87, 93, 98, 97, 101, 102, 103, 104, 105, 108, 109,116, 115, 114, 96, 100,  99, 107, 113]; % which channels to include
    % S.img.chanswap = {[129,55]}; % swap over pairs of channels if needed, or leave empty
    % S.img.bbox  = [0.48 0.5; 0.59 0.70; 0.70 0.70; 0.73 0.79; 0.82 0.79; 0.88 0.7; 0.94 0.56; 0.94 0.45; 0.86 0.25; 0.79 0.13; 0.74 0.28; 0.68 0.25; 0.63 0.35; 0.48 0.44; 0.48 0.5];
    % S.img.bbox_inner   = [0.77 0.42; 0.77 0.47; 0.78 0.5; 0.8 0.52; 0.88 0.5; 0.88 0.38;0.80 0.36; 0.78 0.38; 0.77 0.42];
    % S.img.imgsize = 32;
    % S.img.interp_method = 'gdatav4';

function [YI, YI_mask] = topotime_3D_v2(Y, S)
    % Handle input type (Y can be a structure with filenames or a matrix)
    if isstruct(Y)
        try
            files = Y;
            Y = load(fullfile(files(1).folder, files(1).name));
            fieldname = fieldnames(Y);
            Y = Y.(fieldname{1});
        catch
            error('Input to topotime_3D must be either a double array or a cell array of filenames');
        end
    end
    
    % Assign img structure if not already present
    try S.img; 
    catch
        S.img = S.MCC;
    end

    % Data sizes
    chanind = 1:size(Y,1);  % index of channels
    timeind = 1:size(Y,2);  % index of timepoints
    trialind = 1:size(Y,3); % index of trials
    n = S.img.imgsize;      % dimension (32 is default in SPM)

    % Get channel locations
    if isfield(S.img.file, 'coord')
        temp = load(S.img.file.coord);
        c2d = temp.(S.img.file.coord_struct);
        
        % Handle channel swaps if needed
        if isfield(S.img, 'chanswap') && ~isempty(S.img.chanswap)
            for i = 1:length(S.img.chanswap)
                c2d(1, S.img.chanswap{i}(1)) = c2d(1, S.img.chanswap{i}(2));
            end
        end

        % Create Cel from selected channels and scale coordinates
        sel = c2d(:, S.img.chansel);
        sel=sel';
        Cel = scale_coor(sel', n); % Scale coordinates, 'n' is the grid size
        
        % Create grid
        [x, y] = meshgrid(1:n, 1:n);
        ch = convhull(Cel(:, 1), Cel(:, 2));

        %% NEW - Bounding box handling
        if isfield(S.img, 'bbox') && isfield(S.img, 'bbox_inner')
            MINS = min(S.img.bbox);
            MAXS = max(S.img.bbox);
            xi = linspace(MINS(1,1), MAXS(1,1), n);
            yi = linspace(MINS(1,2), MAXS(1,2), n);
            [xi, yi] = meshgrid(xi, yi);
            
            % indices to keep
            Ik = find(inpolygon(xi, yi, S.img.bbox(:,1), S.img.bbox(:,2)));
            % indices to make NaN
            Ir = find(inpolygon(xi, yi, S.img.bbox_inner(:,1), S.img.bbox_inner(:,2)));
            Ic = setdiff(Ik, Ir);
            xi = xi(Ic); yi = yi(Ic);
        else
            % Default case: indices to keep based on the convex hull
            Ic = find(inpolygon(x, y, Cel(ch, 1), Cel(ch, 2)));
        end
        x = x(Ic); y = y(Ic);

    elseif isfield(S.img.file, 'SPMimg')
        SPMdata = spm_eeg_load(fullfile(S.img.path.SPMdata, S.img.file.SPMdata)); 
        [Cel, x, y] = spm_eeg_locate_channels(SPMdata, n, chanind);
    end
    
    % Preallocate memory for YI (output image)
    dimi = [n, n, length(timeind)];
    if length(trialind) > 1
        dim = [dimi, length(trialind)];
    else
        dim = [dimi];
    end
    YI = NaN(dim);
    
    % Select interpolation method based on S.img.interp_method
    if ~isfield(S.img,'interp_method')
        S.img.interp_method = 'griddata';
    end
    switch S.img.interp_method
        case 'gdatav4'
            for i = trialind
                for j = timeind
                    YY = NaN(n, n);  % Reset YY for each timepoint
                    if length(dim) == 3
                        YY(sub2ind([n n], x, y)) = gdatav4(Cel(:,1), Cel(:,2), Y(:,j), xi, yi);
                        YI(:,:,j) = YY';
                    elseif length(dim) == 4
                        YY(sub2ind([n n], x, y)) = gdatav4(Cel(:,1), Cel(:,2), Y(:,j,i), xi, yi);
                        YI(:,:,j,i) = YY';
                    end
                end
            end
            
        case 'scatteredInterpolant'
            % Create interpolant once with coordinates Cel(:,1), Cel(:,2)
            F = scatteredInterpolant(Cel(:, 1), Cel(:, 2), zeros(size(Cel, 1), 1), 'linear', 'none');
            
            for i = trialind
                for j = timeind
                    YY = NaN(n, n);  % Reset YY for each timepoint
                    
                    % Update the interpolant's values with current Y data
                    if length(dim) == 3
                        F.Values = Y(:, j);
                    elseif length(dim) == 4
                        F.Values = Y(:, j, i);
                    end
                    
                    % Perform the interpolation and fill YY
                    YY(sub2ind([n n], x, y)) = F(x, y);

                    % Store the result in the output matrix YI
                    if length(dim) == 3
                        YI(:, :, j) = YY';
                    else
                        YI(:, :, j, i) = YY';
                    end
                end
            end
            
        otherwise  % Default to griddata if no specific method is set
            for i = trialind
                for j = timeind
                    YY = NaN(n, n);  % Reset YY for each timepoint
                    if length(dim) == 3
                        YY(sub2ind([n n], x, y)) = griddata(Cel(:,1), Cel(:,2), double(Y(:,j)), x, y, 'linear');
                        YI(:,:,j) = YY';
                    elseif length(dim) == 4
                        YY(sub2ind([n n], x, y)) = griddata(Cel(:,1), Cel(:,2), double(Y(:,j,i)), x, y, 'linear');
                        YI(:,:,j,i) = YY';
                    end
                end
            end
    end
    
    % If input was a structure (files), save the results as images
    if exist('files', 'var')
        V = spm_vol(S.img.file.SPMimg);
        for f = 1:length(files)
            loadname = fullfile(files(f).folder, files(f).name);
            Y = load(loadname);
            fieldname = fieldnames(Y);
            Y = Y.(fieldname{1});

            for i = trialind
                for j = timeind
                    YY = NaN(n, n);  % Reset YY for each timepoint
                    if strcmp(S.img.interp_method, 'scatteredInterpolant')
                        F.Values = Y(:, j); % Update values for current data
                        YY(sub2ind([n n], x, y)) = F(x, y);
                    else
                        YY(sub2ind([n n], x, y)) = gdatav4(Cel(:, 1), Cel(:, 2), Y(:, j), xi, yi);
                    end
                    YI(:,:,j) = YY';
                end
            end

            % Save the results as images
            outYI(f).name = fullfile(files(f).folder, strrep(files(f).name, '.mat', '.nii'));
            V.fname = outYI(f).name;
            V.dim = size(YI);
            spm_write_vol(V, YI);
        end
        YI = outYI;
    end

    % Create the mask
    if length(dim) == 3
        YI_mask = ~isnan(YI) & (YI ~= 0);
    elseif length(dim) == 4
        YI_mask = ~isnan(YI(:,:,:,1)) & (YI(:,:,:,1) ~= 0);
    end
end

%==========================================================================
% scale_coor - SPM subfunction from spm_eeg_locate_channels.m
%==========================================================================
function Cel = scale_coor(Cel, n)

% check limits and stretch, if possible
dx = max(Cel(1,:)) - min(Cel(1,:));
dy = max(Cel(2,:)) - min(Cel(2,:));

if dx > 1 || dy > 1
    error('Coordinates not between 0 and 1');
end

scale = (1 - 10^(-6))/max(dx, dy);
Cel(1,:) = n*scale*(Cel(1,:) - min(Cel(1,:)) + eps) + 0.5;
Cel(2,:) = n*scale*(Cel(2,:) - min(Cel(2,:)) + eps) + 0.5;

% shift to middle
dx = n+0.5 -n*eps - max(Cel(1,:));
dy = n+0.5 -n*eps - max(Cel(2,:));
Cel(1,:) = Cel(1,:) + dx/2;
Cel(2,:) = Cel(2,:) + dy/2;

% 2D coordinates in voxel-space (incl. badchannels)
Cel = round(Cel)';

end
