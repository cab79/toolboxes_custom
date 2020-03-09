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
                if S.clus.minpercheight_maxima
                    testval=img{c}(img{c}>0);
                    if ~isempty(testval)
                        BW = imextendedmax(img{c},prctile(testval,S.clus.minpercheight_maxima));
                        img{c}=img{c}.*BW;
                    end
                end
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
            
            
            %% create/load input volume
            types = S.clus.summary_types;
            if ~isempty(types)
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
            end
            input_avg=reshape(mean(input_vol,2),size(img{1},1),size(img{1},2),size(img{1},3));
            
            %% run on con
            pi=0; nclus={}; pimg={};
            for c = S.clus.model.contrast{i}
                
                if ~any(img{c}(:))
                    nclus{c} = 0;
                    continue
                end
                
                D(d).model(i).con(c).vox = [];
                D(d).model(i).con(c).clus = [];
                D(d).model(i).con(c).pca = [];
                
                % find clusters using connected components analysis
                % break up into positive and negative polarity clusters
                posc = img{c}.*input_avg>0;
                cc = ccon(posc,C,S);
                D(d).model(i).con(c).vox = cc.PixelIdxList(1:min(length(cc.PixelIdxList),S.clus.connected_clusters.ClusMaxNum));
                negc = img{c}.*input_avg<0;
                cc = ccon(negc,C,S);
                D(d).model(i).con(c).vox = [D(d).model(i).con(c).vox, cc.PixelIdxList(1:min(length(cc.PixelIdxList),S.clus.connected_clusters.ClusMaxNum))];

                % refine and/or dissect clusters using PCA
                if S.clus.cluster_pca_thresh
                    nc=numel(D(d).model(i).con(c).vox);
                    remove = [];
                    % for each cluster
                    for ci=1:nc
                        cii = D(d).model(i).con(c).vox{ci};
                        disp(['PCA on ' save_pref num2str(d) ', model ' num2str(i)  ', contrast ' num2str(c) ', cluster ' num2str(ci)])
                        X=input_vol(cii,:)';
                        % reduce cluster size for PCA by taking random samples
                        randind=randperm(size(X,2));
                        maxind = min(size(X,2),S.clus.cluster_pca_maxclussize);
                        Xsub=X(:,randind(1:maxind));
                        
                        if S.clus.cluster_pca_thresh==1
                            % Do PCA
                            [coeff,score,~,~,explained] = pca(Xsub,'Algorithm','eig');
                            % Create sub-clusters that explain sufficient variance
                            score_sub = score(:,explained>S.clus.cluster_pca_minvarexplained | explained<-S.clus.cluster_pca_minvarexplained);
                            score_sub = score_sub(:,1:min(size(score_sub,2),round(size(Xsub,2)*S.clus.cluster_pca_fracretainfac)));
                            D(d).model(i).con(c).pca(ci).explained = explained;
                        elseif S.clus.cluster_pca_thresh==2
                            % use erp_pca
                            NUM_FAC = [];
                            FactorResults = erp_pca(Xsub,NUM_FAC);
                            randResults = erp_pca(randn(size(Xsub)),NUM_FAC);
                            randResults.scree = randResults.scree*(sum(FactorResults.scree)/sum(randResults.scree)); %scale to same size as real data
                            % estimate the number of components
                            nfac = find(FactorResults.scree<randResults.scree);
                            % apply PCA to all voxels to get coefficients
                            % and scores
                            FactorResults = erp_pca(X,nfac(1));
                            FactorResults.nfac = nfac(1);
                            
                            % keep all components explaining over e.g. 5%
                            varkeep = FactorResults.facVar>S.clus.cluster_pca_min_var_explained;
                            % keep any others required to reach overall
                            % explained of e.g. 50%
                            varthresh = find(cumsum(FactorResults.facVar)>S.clus.cluster_pca_frac_explained);
                            varkeep(1:varthresh(1))=1;
                            FactorResults.nfac_varkeep = sum(varkeep);
                            
                            if sum(varkeep)==0
                                continue
                            end
                            
                            % To obtain sparse voxels/components, threshold the coefficients by maximising product-of-sums of cluster size over cluster overlap 
                            % first, weight the absolute coefficients by
                            % variance explained
                            abs_coeff = bsxfun(@times,abs(FactorResults.FacCof(:,varkeep)),FactorResults.facVar(1,varkeep));
                            thresh_range = linspace(min(abs_coeff(:)),max(abs_coeff(:)),100);
                            prod_ratio=[];
                            for th = 1:length(thresh_range)
                                thresh=thresh_range(th);
                                coeff_thresh=abs_coeff>thresh;
                                % cluster size: prod of non-zero sums
                                clus_size = sum(coeff_thresh,1);
                                prod_clus_size = prod(clus_size(clus_size>0));
                                % cluster overlap: prod of non-zero sums
                                clus_over = sum(coeff_thresh,2);
                                prod_clus_over = prod(clus_over(clus_over>0));
                                % ratio
                                prod_ratio(th) = prod_clus_size/prod_clus_over;
                            end
                            pratios = find(prod_ratio==max(prod_ratio));
                            thresh = thresh_range(pratios(1));
                            coeff_new = FactorResults.FacCof(:,varkeep).*(abs_coeff>thresh);
                            
                            figure('name',['model ' num2str(i)  ', contrast ' num2str(c) ', cluster ' num2str(ci)]); 
                            subplot(2,2,1);hold on; 
                            yyaxis left
                            plot(FactorResults.scree(1:FactorResults.nfac),'b'); 
                            plot(randResults.scree(1:FactorResults.nfac),'k')
                            line([FactorResults.nfac,FactorResults.nfac],[0,max(FactorResults.scree(1:FactorResults.nfac))],'LineStyle','--','Color','r')
                            line([FactorResults.nfac_varkeep,FactorResults.nfac_varkeep],[0,max(FactorResults.scree(1:FactorResults.nfac))],'LineStyle','--','Color','b')
                            ylabel('scree')
                            yyaxis right
                            bar(1:length(FactorResults.facVar),FactorResults.facVar,'FaceColor',[0.7 0.7 0.7])
                            plot(cumsum(FactorResults.facVar),'g'); 
                            ylabel('cum var explained')
                            subplot(2,2,2);plot(thresh_range,prod_ratio);
                            subplot(2,2,3);imagesc(abs_coeff);
                            subplot(2,2,4);imagesc(abs_coeff>thresh);
                            drawnow
                            pause(1)
                            
                            % calculate new scores
                            if any(coeff_new(:))
                                score_sub = X*coeff_new;
                                remove = [remove ci];
                                ind = length(D(d).model(i).con(c).vox);
                                for ve = 1:size(score_sub,2)
                                    D(d).model(i).con(c).vox{ind+ve}=cii(coeff_new(:,ve)~=0);
                                    D(d).model(i).con(c).clus(ind+ve).eig=score_sub(:,ve)';
                                end
                            end
                            
                            %OLD
%                             score_sub = FactorResults.FacScr(:,FactorResults.facVar>(S.clus.cluster_pca_minvarexplained/100) | FactorResults.facVar<-(S.clus.cluster_pca_minvarexplained/100));
%                             score_sub = score_sub(:,1:min(size(score_sub,2),round(FactorResults.nfac*S.clus.cluster_pca_fracretainfac)));
                            
                            FactorResults.FacCofNew=coeff_new;
                            FactorResults.nFacThresh = sum(any(coeff_new,1));
                            FactorResults.FacCofThresh=thresh;
                            D(d).model(i).con(c).pca(ci).FactorResults = FactorResults;
%                         end
%                         if ~isempty(score_sub)
%                             remove = [remove ci];
%                             ind = length(D(d).model(i).con(c).vox);
%                             for ve = 1:size(score_sub,2)
%                                 r=corr(X,score_sub(:,ve));
%                                 D(d).model(i).con(c).vox{ind+ve}=cii(r>S.clus.cluster_pca_minvoxelcorr);
%                                 D(d).model(i).con(c).clus(ind+ve).eig=score_sub(:,ve)';
%                             end
                        else % store 1st PC
                            [~,score] = pca(Xsub,'Algorithm','eig','NumComponents',1);
                            D(d).model(i).con(c).clus(ci).eig=score'; 
                        end
                    end
                    if ~isempty(remove)
                        D(d).model(i).con(c).vox(remove)=[];
                        D(d).model(i).con(c).clus(remove)=[];
                    end
                    
                
                end
                % re-order with largest cluster first
                [~,si]=sort(cellfun(@length,D(d).model(i).con(c).vox),'descend');
                D(d).model(i).con(c).vox = D(d).model(i).con(c).vox(si);
                if ~isempty(D(d).model(i).con(c).clus)
                    D(d).model(i).con(c).clus = D(d).model(i).con(c).clus(si);
                end
                % remove small clusters
                rmclus = cellfun(@length,D(d).model(i).con(c).vox)<S.clus.connected_clusters.ClusExtentThresh;
                D(d).model(i).con(c).vox(rmclus)=[];
                if ~isempty(D(d).model(i).con(c).clus)
                    D(d).model(i).con(c).clus(rmclus)=[];
                end
                D(d).model(i).con(c).nRemoved = sum(rmclus);
                % reduce number of clusters
                rdclus = 1:min(length(D(d).model(i).con(c).vox),S.clus.connected_clusters.ClusMaxNum);
                D(d).model(i).con(c).nReduced = length(D(d).model(i).con(c).vox) - max(rdclus);
                D(d).model(i).con(c).vox = D(d).model(i).con(c).vox(rdclus);
                if ~isempty(D(d).model(i).con(c).clus)
                    D(d).model(i).con(c).clus = D(d).model(i).con(c).clus(rdclus);
                end
                
                % make image
                pimg{c} = zeros(size(img{c}));
                nclus{c} = length(D(d).model(i).con(c).vox);
                if nclus{c}>0
                    for v=1:nclus{c}
                        pimg{c}(D(d).model(i).con(c).vox{v}) = v;
                    end
                end

            end
            
            if ~isempty(pimg)
                figure('name','contrast clusters')
                rmclus = cellfun(@isempty,nclus);
                nclus(rmclus)=[];
                pimg(rmclus)=[];
                plotind = find(cellfun(@(x) x>0,nclus));
                plotimg = pimg(plotind);
                pterms = cterms(plotind);
                for pi = 1:length(plotimg)
                    subplot(length(plotimg),1,pi)
                    imagesc(reshape(plotimg{pi},[],size(plotimg{pi},3)))
                    colormap(cmp)
                    ylabel(pterms{pi})
                    set(get(gca,'YLabel'),'Rotation',45)
                    set(gca,'xtick',[])
                end
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
                disp(['summary stats for ' save_pref num2str(d) ', model ' num2str(i)])    
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
                            % disp(['means for ' save_pref num2str(d) ', model ' num2str(i)  ', contrast ' num2str(c) ', cluster ' num2str(ci)])
                            D(d).model(i).con(c).clus(ci).input_mean=squeeze(nanmean(input_vol(cii,:),1));
                        end

                        if any(strcmp(types,'eig'))
                            disp(['eigenvariate for ' save_pref num2str(d) ', model ' num2str(i)  ', contrast ' num2str(c) ', cluster ' num2str(ci)])
                            X=input_vol(cii,:)';
                            randind=randperm(size(X,2));
                            maxind = min(size(X,2),S.clus.cluster_pca_maxclussize);
                            Xsub=X(:,randind(1:maxind));
                            [~,score] = pca(Xsub,'Algorithm','eig','NumComponents',1);
                            D(d).model(i).con(c).clus(ci).eig=score';
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

