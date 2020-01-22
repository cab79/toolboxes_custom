
function [YI,YI_mask]=topotime_3D(Y,S)
% create 3D topo-over-time image and mask from chan x time 2D array
% code from spm_eeg_convert2images
% Y = chan x time x trials/contrasts (3D), chan x time (2D), or chan (1D) in rows.
% Or, Y can be a structure of directory/filenames, in which case output is
% also filenames after saving them to disk.

if isstruct(Y)
    try
        files = Y;
        Y = load(fullfile(files(1).folder,files(1).name));
        fieldname = fieldnames(Y);
        Y = Y.(fieldname{1});
    catch
        error('input to topotime_3D must be either a double array or a cell array of filenames')
    end
end

% inputs
try S.img;
catch
    S.img = S.MCC;
end
SPMdata = spm_eeg_load(fullfile(S.img.path.SPMdata,S.img.file.SPMdata)); % example, for chan 2D coords only
chanind = 1:size(Y,1); % index of chans
timeind = 1:size(Y,2); % index of timepoints
trialind = 1:size(Y,3); % index of trials
n = S.img.imgsize; % dimension (default from SPM)

% example data to check topos are correct
% eD = SPMdata.fttimelock.avg;

% get chan locations
[Cel, x, y] = spm_eeg_locate_channels(SPMdata, n, chanind);

% create N.dat
%dat = file_array('', [n n length(timeind)], 'FLOAT32-LE'); 
%cdat       = dat;
dimi = [n n length(timeind)];
dim   = [dimi ones(1, 3-length(dimi)) length(trialind)];
        
% convert
if exist('files','var')
    % create and save images
    V=spm_vol(S.img.file.SPMimg);
    for f = 1:length(files)
        loadname = fullfile(files(f).folder,files(f).name);
        Y = load(loadname);
        fieldname = fieldnames(Y);
        Y = Y.(fieldname{1});
        YI=[];
        %str = ['topotime: creating image ' loadname]; 
        for i = trialind
            %fprintf('%-s : %i / %i',str,i,trialind(end));
            for j = timeind % for each timepoint
                YY = NaN(n,n);

                switch length(dim)
                    case 2
                        YY(sub2ind([n n], x, y)) = griddata(Cel(:,1),Cel(:,2),...
                            double(Y(:)), x, y,'linear');
                        YI(:,:)     = YY;
                    case 3
                        YY(sub2ind([n n], x, y)) = griddata(Cel(:,1),Cel(:,2),...
                            double(Y(:,j)), x, y,'linear');
                        YI(:,:, j)  = YY;
                    case 4
                        YY(sub2ind([n n], x, y)) = griddata(Cel(:,1),Cel(:,2),...
                            double(Y(:,j,i)), x, y,'linear');
                        YI(:,:,j,i) = YY;
                    otherwise
                        error('Invalid output file');
                end
            end
        end
        % save images
        %fprintf('%-s : %s',str,'...saving images');
        outYI(f).name = fullfile(files(f).folder,strrep(files(f).name,'.mat','.nii'));
        V.fname = outYI(f).name;
        V.dim = size(YI);
        spm_write_vol(V,YI);
        %fprintf('%-s : %s \n',str,'...done');
    end
    YI=outYI;
else
    % just create image but don't save
    YI=[];
    %str = ['topotime: creating image']; 
    for i = trialind
        %fprintf('%-s : %i / %i',str,i,trialind(end));
        for j = timeind % for each timepoint
            YY = NaN(n,n);

            switch length(dim)
                case 2
                    YY(sub2ind([n n], x, y)) = griddata(Cel(:,1),Cel(:,2),...
                        double(Y(:)), x, y,'linear');
                    YI(:,:)     = YY;
                case 3
                    YY(sub2ind([n n], x, y)) = griddata(Cel(:,1),Cel(:,2),...
                        double(Y(:,j)), x, y,'linear');
                    YI(:,:, j)  = YY;
                case 4
                    YY(sub2ind([n n], x, y)) = griddata(Cel(:,1),Cel(:,2),...
                        double(Y(:,j,i)), x, y,'linear');
                    YI(:,:,j,i) = YY;
                otherwise
                    error('Invalid output file');
            end
        end
    end
    %fprintf('%-s : %s \n',str,'...done');
    
    % create mask
    switch length(dim)
        case {2,3}
            YI_mask = ~isnan(YI) & (YI~=0);
        case 4
            YI_mask = ~isnan(YI(:,:,:,1)) & (YI(:,:,:,1)~=0);
    end
    
%         % save example eD
%         V=spm_vol(S.file.SPMimg);
%         outYI.name = fullfile(pwd,'test.nii');
%         V.fname = outYI.name;
%         V.dim = size(YI);
%         spm_write_vol(V,YI);
end