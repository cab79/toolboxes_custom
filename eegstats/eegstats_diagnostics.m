function D=eegstats_diagnostics(S,varargin)
dbstop if error

%% find D saved from createimages or MCC functions
if isempty(varargin)
   % try loading the data if not supplied in varargin
   try
       load(fullfile(S.diag.path.inputs,'D.mat'),'D') 
   catch
       error('Supply a D structure or path/file of the data')
   end
else
    D= varargin{1};
end

disp('loading data table')
Ddtab = load(S.diag.path.dtab_inputs);
dtab = Ddtab.D.prep.dtab;

% set paths
S.path.code = {
    };
set_paths(S)

if length(D)>1
    save_pref = 'subject_';
else
    save_pref = 'group_';
end

for d=1:length(D)
    
    types = S.diag.summary_types;
    
    % get model and contrast indices
    if isempty(S.diag.model.index)
        S.diag.model.index = 1:length(D(d).model);
    end
    if isempty(S.diag.model.contrast) && isfield(D(d).model(1),'con')
        for i = S.diag.model.index
            S.diag.model.contrast{i} = 1:length(D(d).model(i).con);
        end
    end
    if isempty(S.diag.model_comp.index)
        try
            S.diag.model_comp.index = 1:length(D(d).model_comp);
        catch
            S.diag.model_comp.index = 0;
        end
    end

    if S.diag.model.index
        for i = S.diag.model.index
            
            %residuals
            if isfield(D(d).model(i),'resid_vol')
                resid_vol=D(d).model(i).resid_vol;
                D(d).model(i).resid_vol=[]; % prevent duplication in memory
            else
                D(d).model(i).resid_vol_file = strrep(D(d).model(i).resid_file,'.mat','_vol.mat');
                if exist(D(d).model(i).resid_vol_file,'file')
                    disp([save_pref num2str(d) ', model ' num2str(i) ', loading residm file'])
                    load(D(d).model(i).resid_vol_file);
                else
                    disp([save_pref num2str(d) ', model ' num2str(i) ', loading residm file'])
                    r=load(D(d).model(i).resid_file);
                    residm=r.resid;
                    clear r;
                    if ~exist('resid_vol','var')
                        disp([save_pref num2str(d) ', model ' num2str(i) ', creating residm image'])
                        resid_vol = topotime_3D(residm,S);
                        disp([save_pref num2str(d) ', model ' num2str(i) ', saving residm image'])
                        save(D(d).model(i).resid_vol_file,'resid_vol','-v7.3');
                        clear residm
                    end
                end
            end
            resid_vol=reshape(resid_vol,[],size(resid_vol,4)); % space-time (voxels) x observations
            for c = S.diag.model.contrast{i}
                nc=numel(D(d).model(i).con(c).vox);
                for ci=1:nc
                    cii = D(d).model(i).con(c).vox{ci};
                    if any(strcmp(types,'median'))
                        % median value from each row (voxel), for each
                        % observation. I.e. voxel can differ depending on the
                        % observation (subject, trial, etc.)
                        D(d).model(i).con(c).clus(ci).resid_median=squeeze(nanmedian(resid_vol(cii,:),1));
                    end

                    if any(strcmp(types,'mean'))
                        D(d).model(i).con(c).clus(ci).resid_mean=squeeze(nanmean(resid_vol(cii,:),1));
                    end
                end
            end
            clear resid_vol

            %fitted
            D(d).model(i).fitted_vol_file = strrep(D(d).model(i).fitted_file,'.mat','_vol.mat');
            if exist(D(d).model(i).fitted_vol_file,'file')
                disp([save_pref num2str(d) ', model ' num2str(i) ', loading fitted file'])
                load(D(d).model(i).fitted_vol_file);
            else
                disp([save_pref num2str(d) ', model ' num2str(i) ', loading fitted file'])
                load(D(d).model(i).fitted_file);
                if ~exist('fitted_vol','var')
                    disp([save_pref num2str(d) ', model ' num2str(i) ', creating fitted image'])
                    fitted_vol = topotime_3D(fitted,S);
                    disp([save_pref num2str(d) ', model ' num2str(i) ', saving fitted image'])
                    save(D(d).model(i).fitted_vol_file,'fitted_vol','-v7.3');
                    clear fitted
                end
            end
            fitted_vol=reshape(fitted_vol,[],size(fitted_vol,4));
            for c = S.diag.model.contrast{i}
                nc=numel(D(d).model(i).con(c).vox);
                for ci=1:nc
                    cii = D(d).model(i).con(c).vox{ci};
                    if any(strcmp(types,'median'))
                        % median value from each row (voxel), for each
                        % observation. I.e. voxel can differ depending on the
                        % observation (subject, trial, etc.)
                        D(d).model(i).con(c).clus(ci).fitted_median=squeeze(nanmedian(fitted_vol(cii,:),1));
                    end

                    if any(strcmp(types,'mean'))
                        D(d).model(i).con(c).clus(ci).fitted_mean=squeeze(nanmean(fitted_vol(cii,:),1));
                    end
                end
            end
            clear fitted_vol

            % plots
            for c = S.diag.model.contrast{i}
                nc=numel(D(d).model(i).con(c).vox);
                fsize=max(5,min(15,100/nc));
                for ci=1:nc
                    close all
                    f=figure('units','normalized','outerposition',[0 0 1 1]);
                    for tp = 1:length(types)
                        if any(~isnan(D(d).model(i).con(c).clus(ci).resid_median))
                            subplot(length(types),3,3*(tp-1)+1)
                            histfit(D(d).model(i).con(c).clus(ci).(['resid_' types{tp}])')
                            subplot(length(types),3,3*(tp-1)+2)
                            boxplot(D(d).model(i).con(c).clus(ci).(['resid_' types{tp}])')
                            subplot(length(types),3,3*(tp-1)+3)
                            qqplot(D(d).model(i).con(c).clus(ci).(['resid_' types{tp}])')
                        end
                    end
                    sname = fullfile(S.diag.path.outputs, [save_pref num2str(d) ', model_' num2str(i) '_con_' num2str(c) '_clus_' num2str(ci) '_residnorm.png']);
                    saveas(f,sname);
                end
                close all
                for tp = 1:length(types)
                    pi(1) = floor(sqrt(nc));
                    pi(2) = ceil(nc/pi(1));

                    % fitted vs response for all clusters
                    f=figure('units','normalized','outerposition',[0 0 1 1]);
                    for ci = 1:nc
                        subplot(pi(1),pi(2),ci)
                        X=D(d).model(i).con(c).clus(ci).(['input_' types{tp}])';
                        Y=D(d).model(i).con(c).clus(ci).(['fitted_' types{tp}])';
                        rho=corr(X,Y,'type','Spearman');
                        scatter(X,Y,20,(0.5+rho/2)*ones(1,3),'filled');
                        rho=num2str(rho);rho=rho(1:4);
                        set(gca,'fontsize',fsize)
                        title(['m' num2str(i) ',sz' num2str(length(D(d).model(i).con(c).vox{ci})) ',r' rho])
                        xlabel('response')
                        ylabel('fitted')
                    end
                    sname=fullfile(S.diag.path.outputs, [save_pref num2str(d) ', model_' num2str(i) '_con_' num2str(c) '_' types{tp} '_fitresp.png']);
                    saveas(f,sname);

                    % fitted vs residuals for all clusters
                    f=figure('units','normalized','outerposition',[0 0 1 1]);
                    for ci = 1:nc
                        subplot(pi(1),pi(2),ci)
                        X=D(d).model(i).con(c).clus(ci).(['fitted_' types{tp}])';
                        Y=D(d).model(i).con(c).clus(ci).(['resid_' types{tp}])';
                        hold on
                        scatter(X,Y,'bx');
                        set(gca,'fontsize',fsize)
                        title(['m' num2str(i) ',sz' num2str(length(D(d).model(i).con(c).vox{ci}))])
                        ylabel('residuals')
                        xlabel('fitted')
                        line([min(X), max(X)],[0 0],'Color','k')
                        if 0
                            yy = smooth(X,Y,0.75,'loess');
                            plot(X,yy,'r--')
                        end
                        hold off
                    end
                    sname=fullfile(S.diag.path.outputs, [save_pref num2str(d) ', model_' num2str(i) '_con_' num2str(c) '_' types{tp} '_residfit.png']);
                    saveas(f,sname);

                    for pd = 1:length(S.diag.pred)
                        X=dtab.(S.diag.pred{pd});
                        % fitted vs residuals for all clusters
                        f=figure('units','normalized','outerposition',[0 0 1 1]);
                        for ci = 1:nc
                            Y=D(d).model(i).con(c).clus(ci).(['resid_' types{tp}])';
                            % max variance ratio
                            uX=unique(X,'stable');
                            vr=[];
                            for x=1:numel(uX)
                                if isnumeric(uX(x))
                                    vr(x) = var(Y(X==uX(x)));
                                elseif iscategorical(uX(x))
                                    X=double(X);uX=double(uX);
                                    vr(x) = var(Y(X==uX(x)));
                                elseif iscell(uX(x))
                                    vr(x) = var(Y(strcmp(X,uX{x})));
                                end
                            end
                            vr_ratio = num2str(max(vr)/min(vr));
                            % boxplot
                            subplot(pi(1),pi(2),ci)
                            boxplot(Y,X);
                            set(gca,'fontsize',fsize)
                            title({['m' num2str(i) ',sz' num2str(length(D(d).model(i).con(c).vox{ci}))],['vRatio ' vr_ratio(1:min(3,length(vr_ratio)))]})
                            ylabel('residuals')
                            xlabel(S.diag.pred{pd})
                        end
                        sname=fullfile(S.diag.path.outputs, [save_pref num2str(d) ', model_' num2str(i) '_con_' num2str(c) '_' types{tp} '_resid' S.diag.pred{pd} '.png']);
                        saveas(f,sname);
                    end
                end
            end

%             delete(D(d).model(i).input_vol_file)
%             try
%                 delete(D(d).model(i).resid_vol_file)
%                 delete(D(d).model(i).fitted_vol_file)
%             end
        end
    end
end
save(fullfile(S.diag.path.outputs, 'D.mat'),'D');


function set_paths(S)
for p = 1:length(S.path.code)
    if S.path.code{p,1}
        addpath(genpath(S.path.code{p,2}));
    else
        addpath(S.path.code{p,2});
    end
end

