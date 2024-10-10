function plot_eegstats_LMM_topo(S)
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

% load cluster data
load(fullfile(S.path.stats_load,'D.mat'),'D');
cd(S.path.stats_load)

if S.plot_erp || S.plot_input 

    nempty = find(~cellfun(@isempty,{D.model(:).input_vol_avg_file}));
    if exist(D.model(nempty(1)).input_vol_avg_file,'file')
        load(D.model(nempty(1)).input_vol_avg_file);
        erp=input_avg;

    else
        disp('loading data table')
        Ddtab = load(S.path.dtab_inputs);
        
        dtab = Ddtab.D.prep.dtab;
        datadim = Ddtab.D.prep.dim;
        
        % build input
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
            save(fullfile(S.path.stats_load,'grand_avg_ERP_rescaled.mat'),'erp');
            clear input_scaled input_scaled_mat input_scaled_chanavg
        end
        clear samples Ddtab dtab
    end
end

% use model comp image as mask?
if ~isempty(S.model_comp_mask.index)
    % find con
    i = S.model_comp_mask.index;
    try D.model(i);
    catch
        error(['this model comp number does not exist'])
    end
    
    % get image
    disp('loading mask file from model comp')
    maskimg = spm_read_vols(spm_vol(D.model_comp(i).MCC.([S.model_comp_mask.img '_img_file'])));
end

cmp=cbrewer('seq', 'Purples', 100, 'pchip');
% colmod = repmat(linspace(0.9,1,100),3,1)';
% cmp=cmp.*colmod; % make darker
% cmp=1-cmp; % make darker
cmp2=cbrewer('div', 'Spectral', 100, 'pchip');
% cmp2 = 1-cmp2;
cmp2 = flipud(rescale(exp(cmp2),0,0.9));
% cmp2 = 1-cmp2;
% cmp2=cmp2.*colmod; % make darker
cmp3=cbrewer('div', 'Spectral', 100, 'pchip');
% cmp3 = 1-cmp3;
cmp3 = flipud(rescale(exp(cmp3),0,0.9));
% cmp3=cmp3.*colmod; % make darker
% add in white for -1 values
% cmp = [1 1 1; cmp];

if ~isempty(S.model.index)
    % find con
    for i = S.model.index
        try D.model(S.model.index)
        catch
            error(['this model number does not exist'])
        end
        if S.plot_fitted %load once per model only
            disp('loading fitted file')
            load(D.model(i).fitted_file);
            fitted_chanavg = nanmean(fitted,3);
            fitted_avg{i} = topotime_3D(fitted_chanavg,S);
            clear fitted fitted_chanavg
        end
        if S.plot_resid
            disp('loading resid file')
            load(D.model(i).resid_file);
            resid_chanavg = nanmean(resid,3);
            resid_avg{i} = topotime_3D(resid_chanavg,S);
            clear resid resid_chanavg
        end
        % get cluster images
        for ci = 1:length(S.model.contrast_term)
            if ~any(strcmp({D.model(i).con(:).term},S.model.contrast_term{ci}))
                error([S.model.contrast_term{ci} ' is not a constrast term in this model'])
            end
            c = find(strcmp({D.model(i).con(:).term},S.model.contrast_term{ci}));
            disp('loading image file')
            if isfield(S.img.path,'swap') && ~isempty(S.img.path.swap)
                newpath = strrep(D.model(i).con(c).MCC.([S.clus_image '_img_file']),S.img.path.swap{1},S.img.path.swap{2});
                cimg_all = spm_read_vols(spm_vol(newpath));
            else
                cimg_all = spm_read_vols(spm_vol(D.model(i).con(c).MCC.([S.clus_image '_img_file'])));
            end
            
            % reduce to only those clusters within D.model(i).con(c).vox
            cimg{ci} = zeros(size(cimg_all));
            idx = sort(vertcat(D.model(i).con(c).vox{intersect(S.model.clus,1:length(D.model(i).con(c).vox))}));
            cimg{ci}(idx) = cimg_all(idx);
        end
        % only include overlaps?
        if S.model.contrast_overlaps
            all_cimg = cat(4,cimg{:});
            all_cimg = all(all_cimg,4);
        end
        % coefficient image
        coeff_term = S.model.coeff_term;
        coeff_idx = find(strcmp({D.model(i).coeff(:).name},coeff_term));
        if isfield(S.img.path,'swap') && ~isempty(S.img.path.swap)
            newpath = strrep(D.model(i).coeff(coeff_idx).b_img_file,S.img.path.swap{1},S.img.path.swap{2});
            coeff_img = spm_read_vols(spm_vol(newpath));
        else
            coeff_img=spm_read_vols(spm_vol(D.model(i).coeff(coeff_idx).b_img_file));
        end
        
        regress_erp = coeff_img;
        % mask
        if ~exist('maskimg','var')
            if isfield(S.img.path,'swap') && ~isempty(S.img.path.swap)
                newpath = strrep(D.model(i).mask_img_file,S.img.path.swap{1},S.img.path.swap{2});
                maskimg = spm_read_vols(spm_vol(newpath));
            else
                maskimg = spm_read_vols(spm_vol(D.model(i).mask_img_file));
            end
            maskimg(maskimg==0)=nan;
        end
        % plot clusters per contrast
        for ci = 1:length(S.model.contrast_term)
            c = find(strcmp({D.model(i).con(:).term},S.model.contrast_term{ci}));
            nclus = S.model.clus;
            nclus(nclus>length(D.model(i).con(c).clus))=[];

            % subplots per cluster
            for cl = nclus

                figure('Name',['model ' num2str(i) ', ' S.model.contrast_term{ci} ', cluster ' num2str(cl)],'units','normalized','outerposition',[0 0 0.75 1]);
                %set(gcf, 'renderer', 'opengl');
                clear ax1 ax2
                % find clusters
                clusvox = D.model(i).con(c).vox{cl};
                clusmask = nan(size(cimg{ci}));
                clusmask(clusvox)=1;
                if exist('maskimg','var')
                    clusmask = clusmask .* maskimg>0; % if using mask
                end
                clusmask0 = zeros(size(cimg{ci}));
                clusmask0(clusvox)=1;
                if exist('maskimg','var')
                    clusmask0 = clusmask0 .* maskimg>0; % if using mask
                end
                % only include overlaps?
                if S.model.contrast_overlaps
                     clusmask =  clusmask .* all_cimg;
                     clusmask0 =  clusmask0 .* all_cimg;
                end
        
                if all(isnan(clusmask(:)) | (clusmask(:)==0))
                    disp(['cluster ' num2str(cl) ' is not within the mask - moving on...'])
                    continue
                end
                
                
                % do the masking
                cimgmask0 = cimg{ci}.*clusmask0;
                cimgmask = cimg{ci}.*clusmask;
%                 coeff_imgmask = coeff_img.*clusmask;

                % reshape for plotting
                clusmask = reshape(clusmask,[],size(clusmask,3));
                cimgmask_re = reshape(cimgmask,[],size(cimgmask,3));
                cimg_re = reshape(cimg{ci},[],size(cimg{ci},3));
%                 coeff_imgmask_re = reshape(coeff_imgmask,[],size(coeff_imgmask,3));

                % plot img values over time
                cimg_re(cimg_re==0) = nan;
                cimgmask_re(cimgmask_re==0) = nan;
%                 coeff_imgmask_re(coeff_imgmask_re==0) = nan;
                switch S.summarytype
                    case 'mean'
                        cimg_avg = nanmean(cimg_re,1);
                        cimg_avgmask = nanmean(cimgmask_re,1);
%                         coeff_avgmask = nanmean(coeff_imgmask_re,1);
                    case 'median'
                        cimg_avg = nanmedian(cimg_re,1);
                        cimg_avgmask = nanmedian(cimgmask_re,1);
%                         coeff_avgmask = nanmedian(coeff_imgmask_re,1);
                end
                cimg_avg(isnan(cimg_avg)) = 0;
                cimg_avgmask(isnan(cimg_avgmask)) = 0;
                cimg_extent_ms{cl} = S.time([find(cimg_avgmask,1,'first'),find(cimg_avgmask,1,'last')]);
%                 coeff_avgmask(isnan(coeff_avgmask)) = 0;
                ax1(1)=subplot(5,6,[1:3]);
                hold on
                if S.clus_image_mask % mask the image to show only one cluster per plot 
                    p1=imagesc(ax1(1),S.time,[],cimg_avgmask); 
                else
                    p1=imagesc(ax1(1),S.time,[],cimg_avg); 
                end
                set(gca, 'YTick', []);
                xlim([S.time(1) S.time(end)])
                if isfield(S,'clus_image_range') && ~isempty(S.clus_image_range)
                    caxis(S.clus_image_range)
                end
        %         ax(2)=subplot(5,6,[1:3]); 
        %         p2=imagesc(ax(2),S.time,[],cimg_avgmask,'AlphaData',double(cimg_avgmask>0));
        %         linkaxes([ax(1) ax(2)]) 
        %         p2=pcolor(S.time,[1:2],repmat(cimg_avgmask,2,1));shading flat;
        %         colormap(ax(2),cbrewer('seq', 'Reds', 100, 'pchip')); 
                colormap(ax1(1),cmp);
                cb=colorbar('Location','eastoutside','FontSize',S.FontSize);
                cb.Label.String='z-value';
%                 cb.Label.Position=[0 S.clus_image_range(2)/2];
%                 cb.Label.HorizontalAlignment='center';
%                 cb.Label.VerticalAlignment='bottom';
%                 cb.Label.FontSize = 8;
                plot([0 0],[0.5,1.5],'k')
                title('Cluster Z-values over time (ms)')
                % add peak highlights
                [~,locs]=findpeaks(cimg_avgmask,S.time,'MinPeakDistance',S.MinPeakDistance);
                
                if S.order_locs_by_extent
                    % sort by cluster extent
                    ext=[];
                    for pk=1:length(locs)
                        ext(pk) = sum(maskimg(:,:,S.time==locs(pk)).*cimgmask0(:,:,S.time==locs(pk))>0,'all');
                    end
                    [~,si] = sort(ext,'descend');
                    locs = locs(si);
                end
                
                D.model(i).con(c).clus(cl).peaks = locs;
                D.model(i).con(c).clus(cl).extent = cimg_extent_ms{cl};
                
                if length(locs)==1
                    plot([locs(1) locs(1)],[0.5,1.5],'k--', 'LineWidth', 1.5)
%                     ti=[min(find(cimg_avgmask))-1,max(find(cimg_avgmask))+1];
%                     ri=S.time([max(ti(1),1),min(ti(end),length(cimg_avgmask))]);
%                     rectangle('Position', [ri(1), 0.5, ri(end)-ri(1), 1], 'EdgeColor', [0.1 0.1 0.1],'LineWidth',2);
                elseif length(locs)>1
                    for pk = 1:min([3,S.Nlocs,length(locs)])
                        plot([locs(pk) locs(pk)],[0.5,1.5],'k--', 'LineWidth', 1.5)
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
                for pk=1:min([2,S.Nlocs,length(locs)])
                    ax2(spi(pk))=subplot(5,6,spi(pk)); 
                    plotimg = maskimg(:,:,S.time==locs(pk)).*cimgmask0(:,:,S.time==locs(pk));
                    pcolor(plotimg), shading interp; axis off
                    colormap(ax2(spi(pk)),cmp); 
                    if isfield(S,'clus_image_range') && ~isempty(S.clus_image_range)
                        caxis(S.clus_image_range)
                    end
                    title([num2str(locs(pk)) ' ms'])
                end
                
                % plot regression ERP
                if S.plot_regressionerp

                    ax1(2)=subplot(5,6,[7:9]);
                    wf = plot_clusterovertime(regress_erp,clusmask,S,cmp2,ax1(2));
                    D.model(i).con(c).clus_rerp_mip{cl}=wf;
                    title('Cluster regression ERP (coefficients)')
                    %%align zero for left and right
                    ylimr = get(gca,'Ylim'); 
                    limy=[-max(abs(ylimr)),max(abs(ylimr))];
                    ylim(limy)
                    plot([0 0],limy,'k')
                    plot([S.time(1) S.time(end)],[0 0],'k-')
                    % lines
                    for pk = 1:min([3,S.Nlocs,length(locs)])
                        plot([locs(pk) locs(pk)],limy,'k--', 'LineWidth', 1.5)
                    end


                    spi=[10:12];
                    for pk=1:min([3,S.Nlocs,length(locs)])
                        ax2(spi(pk))=subplot(5,6,spi(pk)); 
                        plotimg = regress_erp(:,:,S.time==locs(pk));
                        pcolor(plotimg), shading interp; axis off
                        colormap(ax2(spi(pk)),cmp2); 
                        title([num2str(locs(pk)) ' ms'])
                    end

                end

                % plot grand average and cluster-over-time
                if S.plot_erp

                    ax1(3)=subplot(5,6,[13:15]);
                    [wf,mn,mx] = plot_clusterovertime(erp,clusmask,S,cmp3,ax1(3));
                    D.model(i).con(c).clus_erp_mip{cl}=wf;
                    title('Cluster maximum intensity projection')
                    %%align zero for left and right
                    yliml = get(gca,'Ylim'); 
                    limy=[-max(abs(yliml)),max(abs(yliml))];
                    ylim(limy)
                    plot([0 0],limy,'k')
                    plot([S.time(1) S.time(end)],[0 0],'k-')
                    % lines
                    for pk = 1:min([3,S.Nlocs,length(locs)])
                        plot([locs(pk) locs(pk)],limy,'k--', 'LineWidth', 1.5)
                    end

                    spi=[16:18];
                    for pk=1:min([3,S.Nlocs,length(locs)])
                        ax2(spi(pk))=subplot(5,6,spi(pk)); 
                        pcolor(erp(:,:,S.time==locs(pk))), shading interp; axis off
                        colormap(ax2(spi(pk)),cmp3); 
                        title([num2str(locs(pk)) ' ms'])
                    end

                end
%                 if S.plot_input
% 
%                     ax1(3)=subplot(5,6,[13:15]);
%                     [wf,mn,mx] = plot_clusterovertime(input_avg,clusmask,S,'b');
%                     title('input: MIP ERP (blue) and MIP rERP (green)')
%                     plot([0 0],[mn,mx],'k')
%                     plot([S.time(1) S.time(end)],[0 0],'k-')
%                     for pk = 1:min([3,S.Nlocs,length(locs)])
%                         plot([locs(pk) locs(pk)],[mn,mx],'k--', 'LineWidth', 1.5)
%                     end
% 
%                     spi=[16:18];
%                     for pk=1:min([3,S.Nlocs,length(locs)])
%                         ax2(spi(pk))=subplot(5,6,spi(pk)); 
%                         pcolor(input_avg(:,:,S.time==locs(pk))), shading interp; axis off
%                         title([num2str(locs(pk)) ' ms'])
%                     end
% 
%                 end
                if S.plot_fitted

                    ax1(4)=subplot(5,6,[19:21]);
                    [wf,mn,mx] = plot_clusterovertime(fitted_avg{i},clusmask,S,'b');
                    title('fitted: MIP ERP (blue) and MIP rERP (green)')
                    plot([0 0],[mn,mx],'k')
                    plot([S.time(1) S.time(end)],[0 0],'k-')
                    for pk = 1:min([3,S.Nlocs,length(locs)])
                        plot([locs(pk) locs(pk)],[mn,mx],'k--', 'LineWidth', 1.5)
                    end

                    spi=[22:24];
                    for pk=1:min([3,S.Nlocs,length(locs)])
                        ax2(spi(pk))=subplot(5,6,spi(pk)); 
                        pcolor(fitted_avg{i}(:,:,S.time==locs(pk))), shading interp; axis off
                        title([num2str(locs(pk)) ' ms'])
                    end
                end
                if S.plot_resid

                    ax1(5)=subplot(5,6,[25:27]);
                    [wf,mn,mx] = plot_clusterovertime(resid_avg{i},clusmask,S,'b');
                    title('resid: MIP ERP (blue) and MIP rERP (green)')
                    plot([0 0],[mn,mx],'k')
                    plot([S.time(1) S.time(end)],[0 0],'k-')
                    for pk = 1:min([3,S.Nlocs,length(locs)])
                        plot([locs(pk) locs(pk)],[mn,mx],'k--', 'LineWidth', 1.5)
                    end

                    spi=[28:30];
                    for pk=1:min([3,S.Nlocs,length(locs)])
                        ax2(spi(pk))=subplot(5,6,spi(pk)); 
                        pcolor(resid_avg{i}(:,:,S.time==locs(pk))), shading interp; axis off
                        title([num2str(locs(pk)) ' ms'])
                    end
                end
                
                % add spacing
                for a=1:length(ax1)
                    ax1(a).Position(2) = ax1(a).Position(2)-S.vertical_spacing*a;  
                    ax1(a).FontSize = S.FontSize; 
                end
                if exist('ax2','var')
                    for a=1:length(ax2)
                        
                        if ismember(a,4:6) mult = 1;
                        elseif ismember(a,10:12) mult = 2;
                        elseif ismember(a,16:18) mult = 3;
                        elseif ismember(a,22:24) mult = 4;
                        elseif ismember(a,28:30) mult = 5;
                        end
                        
                        try
                            ax2(a).Position(2) = ax2(a).Position(2)-S.vertical_spacing*mult;  
                            ax2(a).FontSize = S.FontSize; 
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
                disp(['model ' num2str(i) ', contrast ' num2str(ci) ', coeff ' num2str(coeff_idx) ', cluster ' num2str(cl) ' extent: ' num2str(cimg_extent_ms{cl}) ' ms'])
                disp(['model ' num2str(i) ', contrast ' num2str(ci) ', coeff ' num2str(coeff_idx) ', cluster ' num2str(cl) ' peaks: ' num2str(locs) ' ms'])
                drawnow
                pause(1)
            end
            
        end
    end
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

function [mipmask,mn,mx] = plot_clusterovertime(avg,clusmask,S,col,ax)
avg_re = reshape(avg,[],size(avg,3));
avgmask_re = avg_re.*clusmask;
mip_mask = repmat(any(clusmask,2),1,size(clusmask,2));
mip = avg_re.*mip_mask;

avg_re(avg_re==0) = nan;
avgmask_re(avgmask_re==0) = nan;
mip(mip==0) = nan;

gfp = nanstd(avg_re,[],1);
switch S.summarytype
case 'mean'
    avgmask = nanmean(avgmask_re,1);
    mipmask = nanmean(mip,1);
case 'median'
    avgmask = nanmedian(avgmask_re,1);
    mipmask = nanmedian(mip,1);
end
gfp(isnan(gfp)) = 0;
avgmask(isnan(avgmask)) = 0;
mipmask(isnan(mipmask)) = 0;

% baseline correct
if ~isempty(S.base)
    bp = dsearchn(S.time',S.base');
    mipmask = mipmask-mean(mipmask(bp));
end
        
hold on
h=plot(S.time,mipmask,'LineWidth', 3);
xlim([S.time(1) S.time(end)])

if ischar(col) || size(col,1)==1
    h.color = col;
else
    % cmap
    cm = colormap(ax,col); 
    maxval = max(abs(mipmask));
    cm = interp1(linspace(-maxval,maxval,length(cm)),cm,mipmask); % map color to y values
    cm = uint8(cm'*255); % need a 4xN uint8 array
    cm(4,:) = 255; % last column is transparency
    drawnow
    set(h.Edge,'ColorBinding','interpolated','ColorData',cm)
end

mn = min([mipmask,avgmask]);
mx = max([mipmask,avgmask]);