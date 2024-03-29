function D=eegstats_diagnostics_channel(S,varargin)
% This version operates on single EEG channels or CCA components

dbstop if error

%% find D 
if isempty(varargin)
   % try loading the data if not supplied in varargin
   try
       load(S.diag.path.inputs,'D') 
   catch
       error('Supply a D structure or path/file of the data')
   end
else
    D= varargin{1};
end

disp('loading data table')
Ddtab = load(S.diag.path.dtab_inputs);

if length(D)>1
    save_pref = 'subject_';
else
    save_pref = 'group_';
end

if isfield(S.diag,'subjects')
    di = S.diag.subjects;
else
    di = 1:length(D);
end

for d=di
    if length(D)>1
        dtab = Ddtab.D(d).prep.dtab;
    else
        dtab = Ddtab.D(1).prep.dtab;
    end
    
%     types = S.diag.summary_types;
    
    % get model and contrast indices
    if isempty(S.diag.model.index)
        S.diag.model.index = 1:length(D(d).model);
    end

    if S.diag.model.index
        for i = S.diag.model.index
            
            disp([save_pref num2str(d) ', model ' num2str(i) ', loading resid file'])
            load(D(d).model(i).resid_file);
            resid=squeeze(resid);

            disp([save_pref num2str(d) ', model ' num2str(i) ', loading fitted file'])
            load(D(d).model(i).fitted_file);
            fitted=squeeze(fitted);

            disp([save_pref num2str(d) ', model ' num2str(i) ', loading input file'])
            load(D(d).model(i).input_file);
            input=squeeze(input);

            % N channels
            if ~isempty(S.diag.channels)
                chan = S.diag.channels;
            else
                chan = 1:size(input,1);
            end
            nc = length(chan);
            fsize=max(5,min(15,100/nc));

            % plots
            f=figure('units','normalized','outerposition',[0 0 1 1]);
            tiledlayout(nc,3)
            for ci=1:nc
%                 close all
                nexttile
                histfit(resid(ci,:))
                ylabel(['channel ' num2str(ci)])
                nexttile
                boxplot(resid(ci,:))
                nexttile
                qqplot(resid(ci,:))
%                 sname = fullfile(S.diag.path.outputs, [save_pref num2str(d) ', model_' num2str(i) '_con_' num2str(c) '_clus_' num2str(ci) '_residnorm.png']);
%                 saveas(f,sname);
            end
%             close all
%             for tp = 1:length(types)
%                 pi(1) = floor(sqrt(nc));
%                 pi(2) = ceil(nc/pi(1));

                % fitted vs response for all channels
                f=figure('units','normalized','outerposition',[0 0 1 1]);
                tiledlayout(nc,1)
                for ci = 1:nc
                    nexttile
                    X=input(ci,:)';
                    Y=fitted(ci,:)';
                    rho=corr(X,Y,'type','Spearman');
                    scatter(X,Y,20,(0.5+rho/2)*ones(1,3),'filled');
                    rho=num2str(rho);rho=rho(1:4);
                    set(gca,'fontsize',fsize)
                    title(['m' num2str(i) ', channel ' num2str(ci) ', rho=' rho])
                    xlabel('response')
                    ylabel('fitted')
                end
%                 sname=fullfile(S.diag.path.outputs, [save_pref num2str(d) ', model_' num2str(i) '_con_' num2str(c) '_' types{tp} '_fitresp.png']);
%                 saveas(f,sname);

                % fitted vs residuals for all channels
                f=figure('units','normalized','outerposition',[0 0 1 1]);
                tiledlayout(nc,1)
                for ci = 1:nc
                    nexttile
                    X=fitted(ci,:)';
                    Y=resid(ci,:)';
                    hold on
                    scatter(X,Y,'bx');
                    set(gca,'fontsize',fsize)
                    title(['m' num2str(i) ', channel ' num2str(ci)])
                    ylabel('residuals')
                    xlabel('fitted')
                    line([min(X), max(X)],[0 0],'Color','k')
                    if S.diag.loess

                        % take a random sample for speed
                        sampleSize = min(length(X),3000);  % Number of points to include in the random sample
                        rng('default');  % For reproducibility
                        randomIndices = randperm(length(X), sampleSize);
                        xSample = X(randomIndices);
                        ySample = Y(randomIndices);
                        
                        % Perform smoothing on the random sample
                        disp(['loess for channel ' num2str(ci)])
                        yy = smooth(xSample, ySample,100,'loess');
                        plot(xSample,yy,'r--')
                    end
                    hold off
                end
%                 sname=fullfile(S.diag.path.outputs, [save_pref num2str(d) ', model_' num2str(i) '_con_' num2str(c) '_' types{tp} '_residfit.png']);
%                 saveas(f,sname);

                for pd = 1:length(S.diag.pred)
                    X=dtab.(S.diag.pred{pd});
                    % fitted vs residuals for each channel
                    f=figure('units','normalized','outerposition',[0 0 1 1]);
                    tiledlayout(nc,1)
                    for ci = 1:nc
                        Y=resid(ci,:)';
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
                        nexttile
                        boxplot(Y,X);
                        set(gca,'fontsize',fsize)
                        title({['m' num2str(i) ', channel ' num2str(ci)],['vRatio ' vr_ratio(1:min(3,length(vr_ratio)))]})
                        ylabel('residuals')
                        xlabel(S.diag.pred{pd})
                    end
%                     sname=fullfile(S.diag.path.outputs, [save_pref num2str(d) ', model_' num2str(i) '_con_' num2str(c) '_' types{tp} '_resid' S.diag.pred{pd} '.png']);
%                     saveas(f,sname);
                end
%             end

        end
    end
end
% save(fullfile(S.diag.path.outputs, 'D.mat'),'D');
