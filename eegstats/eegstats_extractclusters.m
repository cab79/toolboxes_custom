function D=eegstats_extractclusters(S,varargin)
% Save mean or median of the EEG data for each statistic cluster (e.g. thresholded F images from MCC)

dbstop if error
close all

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
   1, 'Q:\MATLAB\toolboxes_external\cbrewer'
    };
set_paths(S)

if length(D)>1
    save_pref = 'subject_';
else
    save_pref = 'group_';
end


cmp=cbrewer('qual', 'Set3', S.clus.connected_clusters.ClusMaxNum, 'pchip');
cmp=[0 0 0;cmp]; % add black as zero

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
        try
            S.clus.model_comp.index = 1:length(D(d).model_comp);
        catch
            S.clus.model_comp.index = 0;
        end
    end

    C = S.clus.connected_clusters.connectivity;
    if S.clus.model.index
        for i = S.clus.model.index
            disp(['extracting clusters for ' save_pref num2str(d) ', model ' num2str(i)])
            for c = S.clus.model.contrast{i}
                
                try
                    img_file = [S.clus.image '_file'];
                    img_file = strrep(img_file,'.nii','');
                    try
                        base_fname=D(d).model(i).con(c).(img_file);
                    catch
                        base_fname=D(d).model(i).con(c).MCC.(img_file);
                    end
                catch
                    img_file = [S.clus.image '_img_file'];
                    img_file = strrep(img_file,'.nii','');
                    try
                        base_fname=D(d).model(i).con(c).(img_file);
                    catch
                        base_fname=D(d).model(i).con(c).MCC.(img_file);
                    end
                end
                img{c} = spm_read_vols(spm_vol(base_fname));
                img{c}(isnan(img{c}))=0;
                
                % 
                testval=img{c}(img{c}>0);
                BW = imextendedmax(img{c},prctile(testval,S.clus.minpercheight_maxima));
                img{c}=img{c}.*BW;
            end
            
            imgUC={};
            cimg = img(S.clus.model.contrast{i});
            temp = {D.model(i).con(:).term};
            cterms = temp(S.clus.model.contrast{i});
            if S.clus.unique_common_clusters
                % obtains, for example with 3 contrasts:
                % 1.	Common to all 3
                % 2.	Common to all sub-pairs, not including common to all 3
                % 3.	Unique to each, not including the above.

                % common to all
                ii=1;
                imgUC{ii} = all(cat(4,cimg{:}),4);
                nameUC{ii} = strjoin(cterms,'_');
                excl_img = imgUC{ii};
                
                % common to sub-combinations or unique
                lc = 1:length(cimg);
                temp_img = excl_img;
                for ci = length(cimg)-1:-1:1
                    comb = nchoosek(lc,ci);
                    for cb = 1:size(comb,1)
                        ii=ii+1; 
                        imgUC{ii} = all(cat(4,cimg{comb(cb,:)}),4);
                        imgUC{ii}(excl_img==1) = 0; % remove common image from previous level
                        temp_img = cat(4,temp_img,imgUC{ii});
                        nameUC{ii} = strjoin(cterms(comb(cb,:)),'_');
                    end
                    excl_img = any(temp_img,4);
                end
                
%                 % plot
%                 figure
%                 pi=0;
%                 for p = 1:length(imgUC)
%                     pi=pi+1;
%                     subplot(length(imgUC),2,pi)
%                     imagesc(squeeze(any(imgUC{p},[1,2]))')
%                     ylabel(nameUC{p})
%                     set(get(gca,'YLabel'),'Rotation',45)
%                     
%                     pi=pi+1;
%                     subplot(length(imgUC),2,pi)
%                     imagesc(reshape(imgUC{p},[],size(imgUC{p},3)))
%                     ylabel(nameUC{p})
%                     set(get(gca,'YLabel'),'Rotation',45)
%                 end
            end
            
            %% run on con
            figure('name','contrast clusters')
            pi=0;
            for c = S.clus.model.contrast{i}
                cc = ccon(img{c},C,S);
                D(d).model(i).con(c).vox = cc.PixelIdxList(1:min(length(cc.PixelIdxList),S.clus.connected_clusters.ClusMaxNum));
                
                % make image
                pimg = zeros(size(img{c}));
                nclus = length(D(d).model(i).con(c).vox);
                for v=1:nclus
                    pimg(D(d).model(i).con(c).vox{v}) = v;
                end
                
                % plot
                pi=pi+1;
                subplot(length(S.clus.model.contrast{i}),1,pi)
                imagesc(reshape(pimg,[],size(pimg,3)))
                colormap(cmp)
                ylabel(cterms{pi})
                set(get(gca,'YLabel'),'Rotation',45)
            end
            
            %% run on conUC
            figure('name','unique and common')
%             pi=0;
            for c = 1:length(imgUC)
                
%                 % remove minimally connected voxels
%                 min_conn=round(S.clus.connected_clusters.connectivity/2);
%                 while 1
%                     x=0;
%                     vval = find(imgUC{c}>0)';
%                     for ii = vval
%                         NeighboursInd = findNeighbours(ii,size(imgUC{c}),S.clus.connected_clusters.connectivity);
%                         if sum(imgUC{c}(NeighboursInd)>0) < min_conn
%                             imgUC{c}(ii)=0;
%                             x=1;
%                         end
%                     end
%                     if x==0
%                         break
%                     end
%                 end
                
                % connected components
                cc = ccon(imgUC{c},C,S);
                D(d).model(i).conUC(c).vox = cc.PixelIdxList(1:min(length(cc.PixelIdxList),S.clus.connected_clusters.ClusMaxNum));
                D(d).model(i).conUC(c).term = nameUC{c};
                
                % make image
                pimg = zeros(size(imgUC{c}));
                nclus = length(D(d).model(i).conUC(c).vox);
                for v=1:nclus
                    pimg(D(d).model(i).conUC(c).vox{v}) = v;
                end
                
                % plot
                subplot(length(imgUC),1,c)
                imagesc(reshape(pimg,[],size(pimg,3)))
                colormap(cmp)
                ylabel(nameUC{c})
                set(get(gca,'YLabel'),'Rotation',45)
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
                            D(d).model(i).con(c).clus(ci).input_median=squeeze(nanmedian(input_vol(cii,:),1));
                        end

                        if any(strcmp(types,'mean'))
                            D(d).model(i).con(c).clus(ci).input_mean=squeeze(nanmean(input_vol(cii,:),1));
                        end
                        
                        % calculate covariance explained by fixed and
                        % random factors - useful for comparing models -
                        % TBC
%                         (sum(diag(lme.CoefficientCovariance),'all')+sum(diag(r{1}),'all'))/(sum(diag(lme.CoefficientCovariance),'all')+lme.MSE+sum(diag(r{1}),'all'));
                    end
                end
                for c = 1:length(imgUC)
                    nc=numel(D(d).model(i).conUC(c).vox);
                    for ci=1:nc
                        cii = D(d).model(i).conUC(c).vox{ci};
                        if any(strcmp(types,'median'))
                            % median value from each row (voxel), for each
                            % observation. I.e. voxel can differ depending on the
                            % observation (subject, trial, etc.)
                            D(d).model(i).conUC(c).clus(ci).input_median=squeeze(nanmedian(input_vol(cii,:),1));
                        end

                        if any(strcmp(types,'mean'))
                            D(d).model(i).conUC(c).clus(ci).input_mean=squeeze(nanmean(input_vol(cii,:),1));
                        end
                        
                        % calculate covariance explained by fixed and
                        % random factors - useful for comparing models -
                        % TBC
%                         (sum(diag(lme.CoefficientCovariance),'all')+sum(diag(r{1}),'all'))/(sum(diag(lme.CoefficientCovariance),'all')+lme.MSE+sum(diag(r{1}),'all'));
                    end
                end
                clear input_vol
            end
        end
    end
end
save(fullfile(S.clus.path.outputs, 'D.mat'),'D');


function set_paths(S)
for p = 1:size(S.path.code,1)
    if S.path.code{p,1}
        addpath(genpath(S.path.code{p,2}));
    else
        addpath(S.path.code{p,2});
    end
end

function cc = ccon(img,C,S)
%% connected components analysis

%                 % find optimal threshold for image by maximising number of
%                 % independent clusters.
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

cc = bwconncomp(img,C);
% re-order with largest cluster first
[~,si]=sort(cellfun(@length,cc.PixelIdxList),'descend');
cc.PixelIdxList = cc.PixelIdxList(si);
% remove small clusters
cc.PixelIdxList(cellfun(@length,cc.PixelIdxList)<S.clus.connected_clusters.ClusExtentThresh)=[];

