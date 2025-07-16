function plot_eegstats_LMM_topo_2D(S)


%% issues
% we have a number of models i, a number of contrasts (within function to prevent loading D multiple times), and a number of
% clusters. Best to instead:
% Select a single model from calling function - if necessary include a loop
% there.
% Select a single contrast, despite being slower - include loop in calling
% function.
% Clusters have to be plotted individually for topos and regression ERPs.
% Ends up with same number of plots as before, unless:
% Identify 4D clusters by re-clustering; get maximal time & frequency from
% each cluster; obtain the topo and plot the mean TF values from within the
% topo cluster.


dbstop if error
% close all
% set paths and file names 
S.path.code = {
   1, 'E:\Q_backup\MATLAB\toolboxes_external\cbrewer'
   0, 'E:\Q_backup\MATLAB\toolboxes_custom\eegstats'
   1, 'E:\Q_backup\MATLAB\toolboxes_external\spm12'
%    0, 'Q:\MATLAB\toolboxes_external' % for suplabel.m
    };

% set paths
set_paths(S)

% set L1 model index
if ~isfield(S,'L1_modeln')
    S.L1_modeln=1;
end

% frequencies
load(S.path.erp_load,'freq_pnts');
freq_npnts = cellfun(@length, freq_pnts);
freq_pnts_select={};
freq_indices={};
for f=1:length(S.run_freq)
    freq = find(ismember(S.freq_bands, S.run_freq{f}));
    cumPnts = cumsum(freq_npnts);
    if freq == 1
        startIndex = 1;
    else
        startIndex = cumPnts(freq-1) + 1;
    end
    endIndex = cumPnts(freq);
    freq_indices{f} = startIndex:endIndex;
    freq_pnts_select{f} = freq_pnts{freq};
end

% plotting times - use highest sampled freq
freq_npnts_select = freq_npnts(ismember(S.freq_bands, S.run_freq));
[~,mi] = max(freq_npnts_select);
S.time = freq_pnts_select{mi};

% L2 data: load, upsample, and combine freq data into 2D array
coeff_4D = [];
mask_4D = [];
cimg_4D = [];
cimg_4D_clusidx=[];
clus_count=0;
for f=1:length(S.run_freq)

    % L2 model to obtain coeffs and clusters
    S.path.stats_load = [S.path.L1 S.modelname '\' S.model.coeff_term_L1 '\' S.run_freq{f} '\eegstats_LMM']; % topo 
    % if ~exist(S.path.stats_load,'dir')
    %     continue
    % end
    load(fullfile(S.path.stats_load,'D.mat'),'D');

    % gather all D structs for every freq, for use getting clusters later -
    % not needed?
    % Df(f)=D;

    % L1 model indices to extract coeffs from this freq
    % sample_idx = freq_indices{f};

    % check contrast exists
    if ~any(strcmp({D.model(S.model.index).con(:).term},S.model.contrast_term))
        error([S.model.contrast_term ' is not a constrast term in this model'])
    end

    % get contrast index
    c = find(strcmp({D.model(S.model.index).con(:).term},S.model.contrast_term));

    % get cluster images
    disp('loading image file')
    if isfield(S.img.path,'swap') && ~isempty(S.img.path.swap)
        newpath = strrep(D.model(S.model.index).con(c).MCC.([S.clus_image '_img_file']),S.img.path.swap{1},S.img.path.swap{2});
        cimg_all = spm_read_vols(spm_vol(newpath));
    else
        cimg_all = spm_read_vols(spm_vol(D.model(S.model.index).con(c).MCC.([S.clus_image '_img_file'])));
    end

    % reduce to only those clusters within D.model(S.model.index).con(c).vox
    c_img = zeros(size(cimg_all));
    cimg_clusidx = zeros(size(cimg_all));

    if ~isempty(D.model(S.model.index).con(c).vox)
        idx = sort(vertcat(D.model(S.model.index).con(c).vox{intersect(S.model.clus,1:length(D.model(S.model.index).con(c).vox))}));
        c_img(idx) = cimg_all(idx);
    
        % to get indices of each cluster
        clus_vox_all = D.model(S.model.index).con(c).vox(intersect(S.model.clus,1:length(D.model(S.model.index).con(c).vox)));
        for cl = 1:length(clus_vox_all)
            clus_count=clus_count+1;
            cimg_clusidx(clus_vox_all{cl})=clus_count;
        end
    end
    
    % get coefficient image
    coeff_idx = find(strcmp({D.model(S.model.index).coeff(:).name},S.model.coeff_term));
    if isfield(S.img.path,'swap') && ~isempty(S.img.path.swap)
        newpath = strrep(D.model(S.model.index).coeff(coeff_idx).b_img_file,S.img.path.swap{1},S.img.path.swap{2});
        coeff_img = spm_read_vols(spm_vol(newpath));
    else
        coeff_img=spm_read_vols(spm_vol(D.model(S.model.index).coeff(coeff_idx).b_img_file));
    end

    % get mask
    if ~exist('maskimg','var')
        if isfield(S.img.path,'swap') && ~isempty(S.img.path.swap)
            newpath = strrep(D.model(S.model.index).mask_img_file,S.img.path.swap{1},S.img.path.swap{2});
            mask_img = spm_read_vols(spm_vol(newpath));
        else
            mask_img = spm_read_vols(spm_vol(D.model(S.model.index).mask_img_file));
        end
        mask_img(mask_img==0)=nan;
    end
    
    % upsample if needed
    if size(mask_img,3)~=length(S.time)
        mask_img = permute(interp1(freq_pnts_select{f}, permute(mask_img,[3 1 2]), S.time, 'nearest', 'extrap'),[2 3 1]);
        c_img = permute(interp1(freq_pnts_select{f}, permute(c_img,[3 1 2]), S.time, 'nearest', 'extrap'),[2 3 1]);
        cimg_clusidx = permute(interp1(freq_pnts_select{f}, permute(cimg_clusidx,[3 1 2]), S.time, 'nearest', 'extrap'),[2 3 1]);
        coeff_img = permute(interp1(freq_pnts_select{f}, permute(coeff_img,[3 1 2]), S.time, 'nearest', 'extrap'),[2 3 1]);
    end

    % combine into 4D TF-topo arrays
    coeff_4D = cat(4,coeff_4D,coeff_img);
    mask_4D = cat(4,mask_4D,mask_img);
    cimg_4D = cat(4,cimg_4D,c_img);
    cimg_4D_clusidx = cat(4,cimg_4D_clusidx,cimg_clusidx);

end
%figure;histogram(cimg_4D_clusidx);

% regression ERP for plotting
regress_erp = coeff_4D;

% % Make this optional later
% % re-cluster the data
% cc = bwconncomp(cimg_4D,conndef(4,"minimal"));
% % cc = bwconncomp(cimg_4D,26); % seems to ignore freq dimension since results in 10 clusters (1 per freq)
% cc.PixelIdxList
% 
% % PCA on clusters?



% L1 data
if S.plot_regressionerp_L1 && isfield(S.path,'inputs_L1') && ~isempty(S.path.inputs_L1)

    % load L1 coeffs
    temp=load(S.path.inputs_L1,'D');
    D_L1 = temp.D;

    % grand avg of L1 coeffs (regression ERP)
    ci = find(strcmp({D_L1(1).model(S.L1_modeln).coeff(:).name}, S.model.coeff_term_L1));
    allL1=[];
    for d = 1:length(D_L1)
        allL1 = cat(3,allL1, D_L1(d).model(S.L1_modeln).coeff(ci).b);
    end
    L1_avg = mean(allL1,3);
    L1_avg = topotime_3D(L1_avg,S);


    L1_avg_4D = [];
    for f=1:length(S.run_freq)

         % upsample if needed
        if length(freq_indices{f})~=length(S.time)
            L1_avg_freq = permute(interp1(freq_pnts_select{f}, permute(L1_avg(:,:,freq_indices{f}),[3 1 2]), S.time, 'nearest', 'extrap'),[2 3 1]);
        else
            L1_avg_freq = L1_avg(:,:,freq_indices{f});
        end
    
        % combine into 4D TF-topo arrays
        L1_avg_4D = cat(4,L1_avg_4D,L1_avg_freq);
    end

end

if S.plot_erp || S.plot_input 

    input_vol_avg_exist=0;
    try
        nempty = find(~cellfun(@isempty,{D.model(:).input_vol_avg_file}));
        input_vol_avg_exist=1;
    catch
        input_vol_avg_exist=0;
    end

    if input_vol_avg_exist && exist(D.model(nempty(1)).input_vol_avg_file,'file')
        load(D.model(nempty(1)).input_vol_avg_file);
        erp=input_avg;

    elseif S.plot_erp && ~S.plot_input && exist(S.path.erp_load,'file')
        load(S.path.erp_load);

    else
        disp('loading data table')
        Ddtab = load(S.path.dtab_inputs);
        
        
        % build input
        if isfield(Ddtab.D(1).prep,'Y')
            dtab = Ddtab.D.prep.dtab;
            if ~exist('freq_pnts','var') && isfield(Ddtab,'freq_pnts')
                freq_pnts = Ddtab.freq_pnts;
            end
            datadim = Ddtab.D.prep.dim;
            samples=struct;
            for s = 1:length(Ddtab.D.prep.Y)
                input=Ddtab.D.prep.Y(s).dtab.data;
                if S.plot_input
                    samples(s).input=input;
                end
                if S.plot_erp
                    if isfield(Ddtab.D.prep.Y(s),'data_std')
                        input_scaled = bsxfun(@plus,bsxfun(@times,input,Ddtab.D.prep.Y(s).data_std),Ddtab.D.prep.Y(s).data_mean); % re-scale
                    else
                        input_scaled = input;
                    end
                    samples(s).input_scaled=input_scaled;
                end
            end
            if S.plot_input
                input_mat = reshape(vertcat([samples(:).input]'),datadim(1),datadim(2),[]);
                input_chanavg = mean(input_mat,3);
                input_avg = topotime_3D(input_chanavg,S);
                clear input input_mat input_chanavg
            end
            if S.plot_erp
                input_scaled_mat = reshape(vertcat([samples(:).input_scaled]'),datadim(1),datadim(2),[]);
                input_scaled_chanavg = mean(input_scaled_mat,3);
                erp = topotime_3D(input_scaled_chanavg,S);
                if exist('freq_pnts','var')
                    save(S.path.erp_load,'erp','freq_pnts');
                else
                    save(S.path.erp_load,'erp');
                end
                clear input_scaled input_scaled_mat input_scaled_chanavg
            end
            clear samples Ddtab dtab
        else
            alldata=[];
            for d = 1:length(Ddtab.D)
                alldata = [alldata; Ddtab.D(d).prep.data];
            end
            datadim = Ddtab.D(1).prep.dim;
            input_mat = reshape(alldata',datadim(1),datadim(2),[]);
            input_chanavg = mean(input_mat,3);
            input_avg = topotime_3D(input_chanavg,S);
            erp = input_avg;
            save(S.path.erp_load,'erp');
            clear input input_mat input_chanavg Ddtab dtab
        end
    end

    % % renoved this part
    % if isfield(S,'freq_indices')
    %     erp = erp(:,:,S.freq_indices);
    % end

    % reshape and resample input into 4D
    ERP_4D = [];
    for f=1:length(S.run_freq)

         % upsample if needed
        if length(freq_indices{f})~=length(S.time)
            ERP_freq = permute(interp1(freq_pnts_select{f}, permute(erp(:,:,freq_indices{f}),[3 1 2]), S.time, 'nearest', 'extrap'),[2 3 1]);
        else
            ERP_freq = erp(:,:,freq_indices{f});
        end
    
        % combine into 4D TF-topo arrays
        ERP_4D = cat(4,ERP_4D,ERP_freq);
    end


end

% % use model comp image as mask? - removed this code for now
% if ~isempty(S.model_comp_mask.index)
%     % find con
%     i = S.model_comp_mask.index;
%     try D.model(S.model.index);
%     catch
%         error(['this model comp number does not exist'])
%     end
% 
%     % get image
%     disp('loading mask file from model comp')
%     maskimg = spm_read_vols(spm_vol(D.model_comp(i).MCC.([S.model_comp_mask.img '_img_file'])));
% end

% specify colormaps
cmp=cbrewer('seq', 'Purples', 100, 'pchip');
cmp2=cbrewer('div', 'Spectral', 100, 'pchip');
cmp2 = flipud(rescale(exp(cmp2),0,0.9));
cmp3=cbrewer('div', 'Spectral', 100, 'pchip');
cmp3 = flipud(rescale(exp(cmp3),0,0.9));

% if S.plot_fitted %load once per model only
%     disp('loading fitted file')
%     load(D.model(S.model.index).fitted_file);
%     fitted_chanavg = nanmean(fitted,3);
%     fitted_avg{i} = topotime_3D(fitted_chanavg,S);
%     clear fitted fitted_chanavg
% end
% if S.plot_resid
%     disp('loading resid file')
%     load(D.model(S.model.index).resid_file);
%     resid_chanavg = nanmean(resid,3);
%     resid_avg{i} = topotime_3D(resid_chanavg,S);
%     clear resid resid_chanavg
% end

% plot clusters per contrast
nclus = length(unique(cimg_4D_clusidx))-1;

% subplots per cluster
for cl = 1:nclus


    % cluster images: create masks to specific cluster of interest
    clusvox = find(cimg_4D_clusidx==cl);
    clusmask = nan(size(cimg_4D));
    clusmask(clusvox)=1;
    if exist('mask_4D','var')
        clusmask = clusmask .* mask_4D>0; % if using mask
    end
    clusmask0 = zeros(size(clusmask));
    clusmask0(clusvox)=1;
    if exist('maskimg','var')
        clusmask0 = clusmask0 .* mask_4D>0; % if using mask
    end
    % % only include overlaps?
    % if S.model.contrast_overlaps
    %      clusmask =  clusmask .* all_cimg;
    %      clusmask0 =  clusmask0 .* all_cimg;
    % end

    if all(isnan(clusmask(:)) | (clusmask(:)==0))
        disp(['cluster ' num2str(cl) ' is not within the mask - moving on...'])
        continue
    end

    % do the masking
    cimgmask0 = cimg_4D.*clusmask0;
    cimgmask = cimg_4D.*clusmask;

    % reshape for plotting
    clusmask_re = reshape(clusmask,[],size(clusmask,3),size(clusmask,4));
    cimgmask_re = reshape(cimgmask,[],size(cimgmask,3),size(cimgmask,4));
    cimg_re = reshape(cimg_4D,[],size(cimg_4D,3),size(cimg_4D,4));

    % Average over topo voxels to obtain time/frequency array for the
    % cluster.
    % first replace zeros with nans to ensure zeros are not included in
    % average, then add the zeros back to the averages.
    cimg_re(cimg_re==0) = nan;
    cimgmask_re(cimgmask_re==0) = nan;
    switch S.summarytype
        case 'mean'
            cimg_avg = squeeze(nanmean(cimg_re,1))';
            cimg_avgmask = squeeze(nanmean(cimgmask_re,1))';
        case 'median'
            cimg_avg = squeeze(nanmedian(cimg_re,1))';
            cimg_avgmask = squeeze(nanmedian(cimgmask_re,1))';
    end
    cimg_avg(isnan(cimg_avg)) = 0;
    % cimg_avgmask(isnan(cimg_avgmask)) = 0;
    
    %% Create figure
    % plot img values over time/frequency
    figure('Name',[S.model.coeff_term_L1 ', model ' num2str(S.model.index) ', ' S.model.contrast_term ', cluster ' num2str(cl)], 'units', 'normalized', 'outerposition', [0 0 0.75 1]);
    
    clear ax1 ax2
    
    % Create the first axes (ax1)
    ax1 = subplot(5, 6, [1:3]);
    im1 = imagesc(ax1, S.time, 1:length(S.run_freq), flipud(cimg_avgmask)); 
    cb = colorbar(ax1, 'Location', 'eastoutside', 'FontSize', S.FontSize);
    cb.Label.String = 'z-value';
    set(ax1, 'Color', 'none');
    im1.AlphaData = ~isnan(flipud(cimg_avgmask));
    hold all
    
    % Create the second axes (ax2) if the mask is disabled
    if ~S.clus_image_mask
        ax2 = axes('Position', get(ax1, 'Position'), 'Color', 'none'); % Align with ax1
        im2 = imagesc(ax2, S.time, 1:length(S.run_freq), flipud(cimg_avg));
        % im2.AlphaData = 0.5; % Set transparency for ax2
        colormap(ax2, flipud(gray)); % Apply gray colormap to ax2
        caxis(ax2, S.clus_image_range); % Adjust color axis for ax2 if specified
        
        % Link axes for consistent interaction
        linkaxes([ax1, ax2]);
        ax2.Visible = 'off'; % Hide ax2 ticks and labels
        ax2.XTick = [];
        ax2.YTick = [];
        
        % Use `uistack` to push ax2 behind ax1
        uistack(ax2, 'bottom'); % Ensure ax2 is behind ax1
    end
    
    % Set the colormap for ax1 (main axes)
    colormap(ax1, cmp);
    
    % Adjust y-ticks and labels for ax1
    yticks(ax1, 1:length(S.run_freq));
    yticklabels(ax1, fliplr(S.run_freq));
    
    % Adjust color axis for ax1 if specified
    if isfield(S, 'clus_image_range') && ~isempty(S.clus_image_range)
        caxis(ax1, S.clus_image_range);
    end
    
    % Add additional plot elements to ax1
    plot(ax1, [0 0], [0.5, length(S.run_freq) + 0.5], 'k');
    title(ax1, 'Cluster Z-values over time (ms)');

    % add peak highlights
    % get cluster extent over time
    cimg_extent_ms{cl} = S.time([find(~isnan(median(cimg_avgmask,1,'omitnan')),1,'first'),find(~isnan(median(cimg_avgmask,1,'omitnan')),1,'last')]);
    locs = findPeakOrPlateau(median(cimg_avgmask,1,'omitnan'),S.time,30);
    % [~,locs]=findpeaks(median(cimg_avgmask,1,'omitnan'),S.time,'MinPeakDistance',30);
    % if isempty(locs)
    %     % use middle latency
    %     locs = S.time(dsearchn(S.time',mean(cimg_extent_ms{cl})));
    % 
    %     % use regression ERP? - ONLY WORKS FOR ERP NOT TF
    % end
    
    if S.order_locs_by_extent
        % sort by cluster extent
        ext=[];
        for pk=1:length(locs)
            ext(pk) = sum(mask_4D(:,:,S.time==locs(pk),:).*cimgmask0(:,:,S.time==locs(pk),:)>0,'all');
        end
        [~,si] = sort(ext,'descend');
        locs = locs(si);
    end
    
    % save locations of peaks
    D.model(S.model.index).con(c).clus(cl).peaks = locs;
    D.model(S.model.index).con(c).clus(cl).extent = cimg_extent_ms{cl};
    
    hold on
    if length(locs)==1
        plot(ax1,[locs(1) locs(1)],[0.5,length(S.run_freq) + 0.5],'k--', 'LineWidth', 1.5)
%                     ti=[min(find(cimg_avgmask))-1,max(find(cimg_avgmask))+1];
%                     ri=S.time([max(ti(1),1),min(ti(end),length(cimg_avgmask))]);
%                     rectangle('Position', [ri(1), 0.5, ri(end)-ri(1), 1], 'EdgeColor', [0.1 0.1 0.1],'LineWidth',2);
    elseif length(locs)>1
        for pk = 1:min([3,S.Nlocs,length(locs)])
            plot(ax1,[locs(pk) locs(pk)],[0.5,length(S.run_freq) + 0.5],'k--', 'LineWidth', 1.5)
%                         if pk<4 % only plot first 3 as rectangles
%                             ti=[locs(pk)-S.MinPeakDistance/2,locs(pk)+S.MinPeakDistance/2];
%                             ri=[max(ti(1),min(S.time)),min(ti(end),max(S.time))];
%                             rectangle('Position', [ri(1), 0.5, ri(end)-ri(1), 1], 'EdgeColor', [0.1 0.1 0.1],'LineWidth',2);
%                         end
        end
    end
    hold off

    % plot image topo
    spi=[4:6];
    for pk=1:min([3,S.Nlocs,length(locs)])
        ax3(spi(pk))=subplot(5,6,spi(pk)); 
        plotimg = max(mask_4D(:,:,S.time==locs(pk),:).*cimgmask0(:,:,S.time==locs(pk),:),[],4);
        pcolor(plotimg), shading interp; axis off
        colormap(ax3(spi(pk)),cmp); 
        if isfield(S,'clus_image_range') && ~isempty(S.clus_image_range)
            caxis(S.clus_image_range)
        end
        title([num2str(locs(pk)) ' ms'])
    end
    
    % plot regression ERP
    if S.plot_regressionerp

        % plot
        ax1(2)=subplot(5,6,[7:9]);
        tf = plot_TF(regress_erp,clusmask_re,S,cmp2,ax1(2));
        title('Cluster TF (coefficients)')

        % set limits
        xlim([S.time(1) S.time(end)])
        limc=[-max(abs(tf(:))),max(abs(tf(:)))];
        colorbar('Location','eastoutside','FontSize',S.FontSize);
        caxis(limc);

        % Adjust y-ticks and labels for ax1
        ylim([0.5, length(S.run_freq)+0.5])
        yticks(1:length(S.run_freq));
        yticklabels(S.run_freq);
        
        % Add additional plot elements to ax1
        plot([0 0], [0.5, length(S.run_freq) + 0.5], 'k');

        % lines
        for pk = 1:min([3,S.Nlocs,length(locs)])
            plot([locs(pk) locs(pk)],[0.5,length(S.run_freq) + 0.5],'k--', 'LineWidth', 1.5)
        end

        % MIP mask for topoplots
        clusmask_2D = reshape(clusmask_re,size(clusmask_re,1),[]);
        mip_mask = repmat(any(clusmask_2D,1),size(clusmask_2D,1),1);
        mip_mask = reshape(mip_mask,size(clusmask_re));
        mip_mask = reshape(mip_mask,size(clusmask));
        mip_mask = mip_mask.*mask_4D;
        mip_mask(mip_mask==0)=nan;

        % topoplots
        spi=[10:12];
        for pk=1:min([3,S.Nlocs,length(locs)])
            ax3(spi(pk))=subplot(5,6,spi(pk)); 
            plotimg = mean(mip_mask(:,:,S.time==locs(pk),:).*regress_erp(:,:,S.time==locs(pk),:),4,'omitnan');
            pcolor(plotimg), shading interp; axis off
            colormap(ax3(spi(pk)),cmp2); 
            caxis(limc);
            title([num2str(locs(pk)) ' ms'])
        end

    end

    % plot grand average and cluster-over-time
    if S.plot_regressionerp_L1 

        % plot
        ax1(3)=subplot(5,6,[13:15]);
        tf = plot_TF(L1_avg_4D,clusmask_re,S,cmp3,ax1(3));
        title('Cluster L1 TF (coefficients)')

        % set limits
        xlim([S.time(1) S.time(end)])
        limc=[-max(abs(tf(:))),max(abs(tf(:)))];
        colorbar('Location','eastoutside','FontSize',S.FontSize);
        caxis(limc);

        % Adjust y-ticks and labels for ax1
        ylim([0.5, length(S.run_freq)+0.5])
        yticks(1:length(S.run_freq));
        yticklabels(S.run_freq);
        
        % Add additional plot elements to ax1
        plot([0 0], [0.5, length(S.run_freq) + 0.5], 'k');

        % lines
        for pk = 1:min([3,S.Nlocs,length(locs)])
            plot([locs(pk) locs(pk)],[0.5,length(S.run_freq) + 0.5],'k--', 'LineWidth', 1.5)
        end

        % MIP mask for topoplots
        clusmask_2D = reshape(clusmask_re,size(clusmask_re,1),[]);
        mip_mask = repmat(any(clusmask_2D,1),size(clusmask_2D,1),1);
        mip_mask = reshape(mip_mask,size(clusmask_re));
        mip_mask = reshape(mip_mask,size(clusmask));
        mip_mask = mip_mask.*mask_4D;
        mip_mask(mip_mask==0)=nan;

        % topoplots
        spi=[16:18];
        for pk=1:min([3,S.Nlocs,length(locs)])
            ax3(spi(pk))=subplot(5,6,spi(pk)); 
            plotimg = mean(mip_mask(:,:,S.time==locs(pk),:).*L1_avg_4D(:,:,S.time==locs(pk),:),4,'omitnan');
            pcolor(plotimg), shading interp; axis off
            colormap(ax3(spi(pk)),cmp2); 
            caxis(limc);
            title([num2str(locs(pk)) ' ms'])
        end

    end

    
    if S.plot_erp

        % plot
        ax1(4)=subplot(5,6,[19:21]);
        tf = plot_TF(ERP_4D,clusmask_re,S,cmp3,ax1(4));
        title('Cluster TF MIP (log power)')

        % set limits
        xlim([S.time(1) S.time(end)])
        limc=[-max(abs(tf(:))),max(abs(tf(:)))];
        colorbar('Location','eastoutside','FontSize',S.FontSize);
        caxis(limc);

        % Adjust y-ticks and labels for ax1
        ylim([0.5, length(S.run_freq)+0.5])
        yticks(1:length(S.run_freq));
        yticklabels(S.run_freq);
        
        % Add additional plot elements to ax1
        plot([0 0], [0.5, length(S.run_freq) + 0.5], 'k');

        % lines
        for pk = 1:min([3,S.Nlocs,length(locs)])
            plot([locs(pk) locs(pk)],[0.5,length(S.run_freq) + 0.5],'k--', 'LineWidth', 1.5)
        end

        % MIP mask for topoplots
        clusmask_2D = reshape(clusmask_re,size(clusmask_re,1),[]);
        mip_mask = repmat(any(clusmask_2D,1),size(clusmask_2D,1),1);
        mip_mask = reshape(mip_mask,size(clusmask_re));
        mip_mask = reshape(mip_mask,size(clusmask));
        mip_mask = mip_mask.*mask_4D;
        mip_mask(mip_mask==0)=nan;

        % topoplots
        spi=[22:24];
        for pk=1:min([3,S.Nlocs,length(locs)])
            ax3(spi(pk))=subplot(5,6,spi(pk)); 
            plotimg = mean(mip_mask(:,:,S.time==locs(pk),:).*ERP_4D(:,:,S.time==locs(pk),:),4,'omitnan');
            pcolor(plotimg), shading interp; axis off
            colormap(ax3(spi(pk)),cmp2); 
            caxis(limc);
            title([num2str(locs(pk)) ' ms'])
        end

    end

    % if S.plot_fitted
    % 
    %     ax1(5)=subplot(5,6,[19:21]);
    %     [wf,mn,mx] = plot_clusterovertime(fitted_avg{i},clusmask,S,'b');
    %     title('fitted: MIP ERP (blue) and MIP rERP (green)')
    %     plot([0 0],[mn,mx],'k')
    %     plot([S.time(1) S.time(end)],[0 0],'k-')
    %     for pk = 1:min([3,S.Nlocs,length(locs)])
    %         plot([locs(pk) locs(pk)],[mn,mx],'k--', 'LineWidth', 1.5)
    %     end
    % 
    %     spi=[22:24];
    %     for pk=1:min([3,S.Nlocs,length(locs)])
    %         ax2(spi(pk))=subplot(5,6,spi(pk)); 
    %         pcolor(fitted_avg{i}(:,:,S.time==locs(pk))), shading interp; axis off
    %         caxis(limy);
    %         title([num2str(locs(pk)) ' ms'])
    %     end
    % end
    % if S.plot_resid
    % 
    %     ax1(6)=subplot(5,6,[25:27]);
    %     [wf,mn,mx] = plot_clusterovertime(resid_avg{i},clusmask,S,'b');
    %     title('resid: MIP ERP (blue) and MIP rERP (green)')
    %     plot([0 0],[mn,mx],'k')
    %     plot([S.time(1) S.time(end)],[0 0],'k-')
    %     for pk = 1:min([3,S.Nlocs,length(locs)])
    %         plot([locs(pk) locs(pk)],[mn,mx],'k--', 'LineWidth', 1.5)
    %     end
    % 
    %     spi=[28:30];
    %     for pk=1:min([3,S.Nlocs,length(locs)])
    %         ax2(spi(pk))=subplot(5,6,spi(pk)); 
    %         pcolor(resid_avg{i}(:,:,S.time==locs(pk))), shading interp; axis off
    %         caxis(limy);
    %         title([num2str(locs(pk)) ' ms'])
    %     end
    % end
    
    % add spacing
    for a=1:length(ax1)
        ax1(a).Position(2) = ax1(a).Position(2)-S.vertical_spacing*a;  
        ax1(a).FontSize = S.FontSize; 
    end
    if exist('ax2','var')
        ax2.Position(2) = ax2.Position(2)-S.vertical_spacing;  
        ax2.FontSize = S.FontSize; 
    end
    if exist('ax3','var')
        for a=1:length(ax3)
            
            if ismember(a,4:6) mult = 1;
            elseif ismember(a,10:12) mult = 2;
            elseif ismember(a,16:18) mult = 3;
            elseif ismember(a,22:24) mult = 4;
            elseif ismember(a,28:30) mult = 5;
            end
            
            try
                ax3(a).Position(2) = ax3(a).Position(2)-S.vertical_spacing*mult;  
                ax3(a).FontSize = S.FontSize; 
            end
        end
    end
    
    % match the x axes
    pos1=get(ax1(1),'Position');
    for a = 2:length(ax1)
        try
            pos=get(ax1(a),'Position');
            pos(3)=pos1(3);
            set(ax1(a),'Position',pos);
        end
    end
    if exist('ax2','var')
        pos=get(ax1(1),'Position');
        set(ax2,'Position',pos);
    end
    disp(['model ' num2str(S.model.index) ', contrast ' num2str(ci) ', coeff ' num2str(coeff_idx) ', cluster ' num2str(cl) ' extent: ' num2str(cimg_extent_ms{cl}) ' ms'])
    disp(['model ' num2str(S.model.index) ', contrast ' num2str(ci) ', coeff ' num2str(coeff_idx) ', cluster ' num2str(cl) ' peaks: ' num2str(locs) ' ms'])
    drawnow

    sgtitle([S.model.coeff_term_L1 ', cluster ' num2str(cl)])


    pause(1)

    exportgraphics(gcf, fullfile(S.path.stats_load, [S.model.coeff_term_L1 ', model ' num2str(S.model.index) '_' S.model.contrast_term '_cluster' num2str(cl) '.png']), 'Resolution', 300); % Adjust resolution as needed
end    

    

save(fullfile(S.path.stats_load,'D.mat'),'D','-v7.3');

function set_paths(S)
for p = 1:size(S.path.code,1)
    if S.path.code{p,1}
        addpath(genpath(S.path.code{p,2}));
    else
        addpath(S.path.code{p,2});
    end
end

function [mipmask,mn,mx] = plot_TF(avg,clusmask,S,col,ax)
avg_re = reshape(avg,[],size(avg,3),size(avg,4));

% mask for TF values in cluster only
avgmask_re = avg_re.*clusmask;

% mask for maximum intensity projection across time and frequency
clusmask_2D = reshape(clusmask,size(clusmask,1),[]);
mip_mask = repmat(any(clusmask_2D,2),1,size(clusmask_2D,2));
mip_mask = reshape(mip_mask,size(avg_re));
mip = avg_re.*mip_mask;

avg_re(avg_re==0) = nan;
avgmask_re(avgmask_re==0) = nan;
mip(mip==0) = nan;

gfp = squeeze(std(avg_re,[],1,'omitnan'));
switch S.summarytype
    case 'mean'
        avgmask = squeeze(mean(avgmask_re,1,'omitnan'))';
        mipmask = squeeze(mean(mip,1,'omitnan'))';
    case 'median'
        avgmask = squeeze(median(avgmask_re,1,'omitnan'))';
        mipmask = squeeze(median(mip,1,'omitnan'))';
end
gfp(isnan(gfp)) = 0;
avgmask(isnan(avgmask)) = 0;
mipmask(isnan(mipmask)) = 0;

% baseline correct
if ~isempty(S.base)
    bp = dsearchn(S.time',S.base');
    mipmask = mipmask-repmat(mean(mipmask(bp,:)),size(mipmask,1),1);
end

mn = min([mipmask,avgmask]);
mx = max([mipmask,avgmask]);
        
% hold on
% h=plot(S.time,mipmask,'LineWidth', 3);
% xlim([S.time(1) S.time(end)])

% plot values over time/frequency
hold all
h = imagesc(ax, S.time, 1:length(S.run_freq),mipmask); 
    
if ischar(col) || size(col,1)==1
    h.color = col;
else
    % cmap
    colormap(ax,col); 
    % maxval = max(abs(mipmask(:)));
    % cm = interp1(linspace(-maxval,maxval,length(cm)),cm,mipmask); % map color to y values
    % cm = uint8(cm'*255); % need a 4xN uint8 array
    % cm(4,:) = 255; % last column is transparency
    % drawnow
    % set(h.Edge,'ColorBinding','interpolated','ColorData',cm)
end

function peakIndex = findPeakOrPlateau(data,times,MinPeakDistance)
    
    % Check if the data is constant (excluding NaNs)
    validData = data(~isnan(data)); % Remove NaN values for checking
    if isempty(validData)
        % If all values are NaN, return NaN (or handle as you see fit)
        peakIndex = NaN;
        return;
    end
    if all(validData == validData(1))
        % If all valid values are the same, return the middle index of the full array
        peakIndex = times(round(length(data) / 2));
        return;
    end
    
    % Use findpeaks to detect peaks
    [peaks, locs] = findpeaks(data,times,'MinPeakDistance',MinPeakDistance);
    
    if isempty(peaks)
        % If no peaks are found, treat it as a plateau
        % Find the start and end of the plateau
        plateauStart = find(data == max(data), 1, 'first');
        plateauEnd = find(data == max(data), 1, 'last');
        
        % Return the middle index of the plateau
        peakIndex = times(round((plateauStart + plateauEnd) / 2));
    else
        % If there are peaks, return the index of the highest peak
        [~, maxPeakIdx] = sort(peaks,'descend');
        peakIndex = locs(maxPeakIdx);
    end

