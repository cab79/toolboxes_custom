function D=eegstats_extractclusters(S,varargin)
% Save mean or median of the EEG data for each statistic cluster (e.g. thresholded F images from MCC)

dbstop if error

%% find D saved from createimages or MCC functions
if isempty(varargin)
   % try loading the data if not supplied in varargin
   try
       load(fullfile(S.clus.path.inputs,'D.mat'),'D') 
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
    if isempty(S.clus.model.index)
        S.clus.model.index = 1:length(D(d).model);
    end
    if isempty(S.clus.model.contrast) && isfield(D(d).model(1),'con')
        for i = S.clus.model.index
            S.clus.model.contrast{i} = 1:length(D(d).model(i).con);
        end
    end
    if isempty(S.clus.model_comp.index)
        S.clus.model_comp.index = 1:length(D(d).model_comp);
    end

    if S.clus.model.index
        for i = S.clus.model.index
            disp(['extracting clusters for ' save_pref num2str(d) ', model ' num2str(i)])
            for c = S.clus.model.contrast{i}
                
                %% connected components analysis
                % find optimal threshold for image by maximising number of
                % independent clusters.
                C = S.clus.connected_clusters.connectivity;
                
                try
                    img_file = [S.thresh.image '_file'];
                    img_file = strrep(img_file,'.nii','');
                    try
                        base_fname=D(d).model(i).con(c).(img_file);
                    catch
                        base_fname=D(d).model(i).con(c).MCC.(img_file);
                    end
                catch
                    img_file = [S.thresh.image '_img_file'];
                    img_file = strrep(img_file,'.nii','');
                    try
                        base_fname=D(d).model(i).con(c).(img_file);
                    catch
                        base_fname=D(d).model(i).con(c).MCC.(img_file);
                    end
                end
                img = spm_read_vols(spm_vol(base_fname));

%                 if any(vrange)
%                     clear cc ncc
%                     for vr = 1:Nthresh
%                         cc(vr) = bwconncomp(threshimgs{vr},C);
%                         ncc(vr) = length(cc(vr).PixelIdxList); %number of cc
%                         scc{vr} = cellfun(@length,cc(vr).PixelIdxList); % size of each cc
%                         sscc(vr) = sum(scc{vr}); % sum of sizes
%                         pscc(vr) = prod(scc{vr}); % product of sizes
%                         npscc(vr) = pscc(vr)*ncc(vr); % product of cc sizes x number of cc
%                     end
% 
%                     % select a threshold
%                     switch S.clus.connected_clusters.maximise
%                         case 'num'
%                             [~,mi]=max(ncc);
%                         case 'prodsize'
%                             [~,mi]=max(pscc);
%                         case 'num_prodsize'
%                             [~,mi]=max(npscc);
%                     end
% 
%                     % calculate percentage of cluster voxels remaining
%                     pc_vox = sscc(mi)/sscc(1);
%                     % select cc
%                     cc = cc(mi);
                    
                    % cc
                    cc = bwconncomp(img,C);
                    % re-order with largest cluster first
                    [~,si]=sort(cellfun(@length,cc.PixelIdxList),'descend');
                    cc.PixelIdxList = cc.PixelIdxList(si);
                    % remove small clusters
                    cc.PixelIdxList(cellfun(@length,cc.PixelIdxList)<S.clus.connected_clusters.ClusExtentThresh)=[];
                    D(d).model(i).con(c).vox = cc.PixelIdxList(1:min(length(cc.PixelIdxList),S.clus.connected_clusters.ClusMaxNum));
                    D(d).model(i).con(c).pc_vox = pc_vox;
%                 end
                
            end

            %% summarise EEG data within clusters
            types = S.clus.summary_types;
            if ~isempty(types)
                % create input volume
                D(d).model(i).input_vol_file = strrep(D(d).model(i).input_file,'.mat','_vol.mat');
                if exist(D(d).model(i).input_vol_file,'file')
                    disp([save_pref num2str(d) ', model ' num2str(i) ', loading input file'])
                    load(D(d).model(i).input_vol_file);
                else
                    disp([save_pref num2str(d) ', model ' num2str(i) ', loading input file'])
                    temp=load(D(d).model(i).input_file);
                    if isfield(temp,'input_vol')
                        input_vol=temp.input_vol;
                    else
                        disp([save_pref num2str(d) ', model ' num2str(i) ', creating input image'])
                        input_vol = topotime_3D(temp.input,S);
                        if S.clus.save_vols
                            disp([save_pref num2str(d) ', model ' num2str(i) ', saving input image'])
                            save(D(d).model(i).input_vol_file,'input_vol','-v7.3');
                            %delete(D(d).model(i).input_file)
                        end
                    end
                    clear temp
                end
                input_vol=reshape(input_vol,[],size(input_vol,4));
                for c = S.clus.model.contrast{i}
                    nc=numel(D(d).model(i).con(c).vox);
                    for ci=1:nc
                        cii = D(d).model(i).con(c).vox{ci};
                        if any(strcmp(types,'median'))
                            % median value from each row (voxel), for each
                            % observation. I.e. voxel can differ depending on the
                            % observation (subject, trial, etc.)
                            D(d).model(i).con(c).clus(ci).input_median=squeeze(median(input_vol(cii,:),1));
                        end

                        if any(strcmp(types,'mean'))
                            D(d).model(i).con(c).clus(ci).input_mean=squeeze(mean(input_vol(cii,:),1));
                        end
                    end
                end
                clear input_vol
            end
        end
    end
end
save(fullfile(S.clus.path.outputs, 'D.mat'),'D');


function set_paths(S)
for p = 1:length(S.path.code)
    if S.path.code{p,1}
        addpath(genpath(S.path.code{p,2}));
    else
        addpath(S.path.code{p,2});
    end
end
