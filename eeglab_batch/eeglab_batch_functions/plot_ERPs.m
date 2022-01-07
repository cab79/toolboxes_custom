function S=plot_ERPs(S,varargin)

dbstop if error

if nargin>1
    type = varargin{1};
else
    type = 'subject';
end

S.func = 'ploterp';
switch S.(S.func).select.datatype
    case 'ERP'
        S.path.file = S.path.erp;
    case 'TF'
        S.path.file = S.path.tf;
    case 'Freq'
        S.path.file = S.path.freq;
end

% load directory
if ~isfield(S.(S.func),'loaddir')
    S.(S.func).loaddir = fullfile(S.path.file,S.(S.func).load.suffix{:});
end
S.path.file = S.(S.func).loaddir;

% select channels
S=select_chans(S);

switch type
    case 'subject'

    % GET FILE LIST
    S = getfilelist(S);

    % show data quality from metrics
    pdata = readtable(S.path.datfile);
%     pdata(all(cellfun(@(x) any(isnan(x)),pdata),2),:) = []; % remove NaN rows
%     pdata(:,all(cellfun(@(x) any(isnan(x)),pdata),1)) = []; % remove NaN rowscols

    % get subjects and include
%     subs = pdata.Subject;
%     include = pdata.Include;
%     subs=subs(include);

%     quality_col = strfind(pdata(1,:),'EEG_quality_');
%     quality_col = find(~cellfun(@isempty,quality_col));
%     metric={};
%     if ~isempty(quality_col)
%         for i = 1:length(quality_col)
%             metric{i,1}=pdata{1,quality_col(i)};
%             metric{i,1}=strrep(metric{i,1},'EEG_quality_','');
%             metric{i,2}=[pdata{2:end,quality_col(i)}];
%             metric{i,2} = metric{i,2}(include);
%         end
%     end

    for i = S.(S.func).startfile : length(S.(S.func).filelist)

%         if exist(fullfile(S.path.file,'data_quality.mat'),'file')
%             load(fullfile(S.path.file,'data_quality'),'qual')
%         end
        file = S.(S.func).filelist{i};
        filesubs=S.(S.func).designtab.subjects;
        load(file)
        
        if exist('gadata','var')
            tldata{1}=gadata;
        else
            % select events
            if iscell(S.ploterp.select.events)
                for ev = 1:length(S.ploterp.select.events)
                    newdata{ev} = tldata{1};
                    selectdata = tldata(S.ploterp.select.events{ev});
                    datstruct = cell2mat(selectdata(~cellfun(@isempty,selectdata)));
                    datmat = cat(1,datstruct(:).trial);
                    newdata{ev}.trial = datmat;
                    newdata{ev}.avg = squeeze(mean(newdata{ev}.trial,1));
                end
                tldata=newdata;
            else
                temp=intersect(S.ploterp.select.events,find(~cellfun(@isempty,tldata)));
                tldata = tldata(temp);
            end
        end
        
        % calculate weighted mean
        noemp = find(~cellfun(@isempty,tldata));
        tldata = tldata(noemp);
        avgdata = tldata{1};
        datstruct = cell2mat(tldata);
        datmat = cat(3,datstruct(:).avg);
        weights = cat(3,datstruct(:).ntrials);
        weighted_data = bsxfun(@times,datmat,weights/mean(weights));
        avgdata.avg = mean(weighted_data,3);
        %avgdata.avg = sum(datmat*weights,2)./sum(datmat,2);

        [f,f2]=plotmulti(S,datmat,tldata,file,avgdata)
        waitfor(f)
        close(f2)
%         % show outlier-ness
%         if exist('outlist','var');
%             max([outlist{:,2}])
%             %outval = round(100*outlist{strcmp(outlist(:,1),subs{i}),2}/max([outlist{:,2}]));
%             outval = round(invprctile([outlist{:,2}],outlist{strcmp(outlist(:,1),subs{i}),2}));
%         else
%             outval = 0;
%         end

%         % good, so-so, bad data?
%         hm=figure('Name','Z-scores of each metric',...
%         'units','normalized','outerposition',[0.75 0.60 0.25 0.40]);
%         hold on
%         col='bgrmy'; 
%         ii=find(strcmp(subs,filesubs{i}));
%         qstr = '';
%         for m = 1:length(metric)
%             if all(ismember(metric{m,2},[0 1]))
%                 %qstr = [qstr metric{m,1} ': ' num2str(metric{m,2}(ii)) ', '];
%             else
%                 %qstr = [qstr metric{m,1} ': ' num2str(round(invprctile(metric{m,2},metric{m,2}(ii)))) '%, '];
%                 sortm=sort(metric{m,2});
%                 p(m)=plot(sortm,col(m)); 
%                 scatter(find(sortm==metric{m,2}(ii)),metric{m,2}(ii),col(m),'filled');
%             end
%         end
%         %thresholds
%         qual_metric = metric{strcmp(metric(:,1),S.(S.func).qual_metric),2};
%         stat(1).data = qual_metric(strcmp(qual(:,2),'Good'));
%         stat(2).data = qual_metric(strcmp(qual(:,2),'So-so'));
%         stat(3).data = qual_metric(strcmp(qual(:,2),'Bad'));
%         if length(stat(1).data)>=S.ploterp.minN_thresh && length(stat(2).data)>=S.ploterp.minN_thresh && length(stat(3).data)>=S.ploterp.minN_thresh
%             for st = 1:length(stat)
%                 stat(st).mean=mean(stat(st).data);
%                 stat(st).std=std(stat(st).data);
%             end
%             % create thresholds based on means weighted by std
%             std_norm = 2*[stat(2).std,stat(1).std]/(stat(1).std+stat(2).std);
%             thresh(1)=mean(std_norm.*[stat(1).mean,stat(2).mean]);
%             perc(1)=dsearchn(sort(qual_metric)',thresh(1));
%             std_norm = 2*[stat(3).std,stat(2).std]/(stat(2).std+stat(3).std);
%             thresh(2)=mean(std_norm.*[stat(2).mean,stat(3).mean]);
%             perc(2)=dsearchn(sort(qual_metric)',thresh(2));
%             plot([0 length(subs)],[thresh(1) thresh(1)],'--k'); % plot thresholds
%             plot([0 length(subs)],[thresh(2) thresh(2)],'--k'); % plot thresholds
%             yl=ylim;
%             plot([perc(1) perc(1)],[yl(1) yl(2)],'--k'); % plot thresholds
%             plot([perc(2) perc(2)],[yl(1) yl(2)],'--k'); % plot thresholds
%         end
%         legend(p,metric(:,1),'Location','southeast');
%         xlabel('subjects')
%         ylabel('z-scores')
%         hold off
%         qual{i,1} = file;
%         if ~isempty(qual{i,2})
%             dft=qual{i,2};
%         else
%             dft='So-so';
%         end
%         qual{i,2} = MFquestdlg([0.9 0.3],qstr,'Data quality','Good','So-so','Bad',dft);
%         if isempty(qual{i,2})
%             return
%         end
%         close(f); close(f2); close(hm);
%         qual=qual;
%         autoind(qual_metric<thresh(1))=1;
%         autoind(qual_metric>=thresh(1) & qual_metric<thresh(2))=2;
%         autoind(qual_metric>=thresh(2))=3;
%         autoqual = {'Good','So-so','Bad'}; autoqual = autoqual(autoind);
%         save(fullfile(S.path.file,'data_quality'),'qual','thresh','perc','autoqual')
% 
%         if strcmp(qual{i,2},'')
%             return
%         end

    end
    
    case 'grandavg'

        fname = dir(fullfile(S.path.file,[S.ploterp.select.sessions{:} '_' S.ploterp.select.blocks{:} '_' S.ploterp.select.conds{:} '_' S.ploterp.load.suffix{:} '.' S.ploterp.fname.ext{:}]));
        fname = fname.name;
        load(fullfile(S.path.file,fname));
        avgdata = gadata.gavg;
        tldata = gadata.events;
        datstruct = cell2mat(tldata);
        datmat = cat(3,datstruct(:).avg);

        if S.ploterp.allsubjects
            title = 'grand average';
            [f,f2]=plotmulti(S,datmat,tldata,title,avgdata)
        end
        
        if S.ploterp.allsubjects_var
            [f3,f4]=plotvar(S,datmat,tldata,title,avgdata)
        end
        
        if S.ploterp.robustsubjects
            tldata = gadata.events_acc;
            title = 'grand average robust';
            [f5,f6]=plotmulti(S,datmat,tldata,title,avgdata)
        end

        if S.ploterp.outliersubjects
            tldata = gadata.events_rej;
            title = 'grand average rejected';
            [f7,f8]=plotmulti(S,datmat,tldata,title,avgdata)
        end
 
end

function [f,f2] = plotmulti(S,datmat,tldata,file,avgdata)
% check chan inputs
if length(S.(S.func).inclchan) ~= length(avgdata.label)
    disp(['no. ERP channels = ' num2str(length(avgdata.label))])
    disp(['no. selected channels = ' num2str(length(S.(S.func).inclchan))])
    error('number of selected channels do not match the number of channels in ERP data')
end
f=figure('units','normalized','outerposition',[0 0 1 1]);
cfg = [];
cfg.layout = S.ploterp.layout;
cfg.ylim = [prctile(datmat(:),0.1),prctile(datmat(:),99.9)];%S.ploterp.ylim;
if ~isfield(S.ploterp,'event_labels') || isempty(S.ploterp.event_labels{:})
    temp=1:length(tldata);
    labels = cellstr(num2str(temp(:)))';
else
    labels = S.ploterp.event_labels;
end
if isfield(tldata{1},'trial')
    for d = 1:length(tldata)
        cfg.dataname{d} = ['event: ' labels{d} ', num trials: ' num2str(size(tldata{d}.trial,1))];
    end
else
    for d = 1:length(tldata)
        cfg.dataname{d} = ['event: ' labels{d}];
    end
end
cfg.graphcolor=S.ploterp.select.eventcolors;
ft_multiplotER_cab(cfg, tldata);
title(file)

f2 = figure('units','normalized','outerposition',[0 1-(0.2*length(S.ploterp.times)) 0.1 0.2*length(S.ploterp.times)]);
for t = 1:length(S.ploterp.times)
    subplot(length(S.ploterp.times),1,t);
    tim = dsearchn(avgdata.time',S.ploterp.times{t}');
    topoplot(mean(avgdata.avg(:,tim(1):tim(2)),2),S.(S.func).chanlocs(S.(S.func).inclchan));
    title([num2str(S.ploterp.times{t}(1)) ':' num2str(S.ploterp.times{t}(2))])
end
%cfg = [];                            
%cfg.xlim = [0.3 0.5];                
%cfg.zlim = [0 6e-14];                
%cfg.layout = S.ploterp.layout;            
%cfg.parameter = 'avg'; % the default 'avg' is not present in the data
%figure; ft_topoplotER(cfg,avgdata); 


function [f3,f4] = plotvar(S,datmat,tldata,file,avgdata)

% settings
lat_ms = [268]; %ms
pos_neg_peak = [1,-1]; %1 or -1
num_ele = 3;
plotidx=[1 2];
col=colormap(jet(length(tldata)));

for i = 1:length(S.ploterp.times)
    tim{i} = dsearchn(avgdata.time',S.ploterp.times{i}');
end

for p = 1:length(pos_neg_peak)

    for i = 1:length(S.ploterp.times)
        
        % identify peak electrodes
        if pos_neg_peak(p)==1
            [~,tps] = sort(mean(avgdata.avg(:,tim{i}(1):tim{i}(2)),2),'descend');
            tp{i} = tps(1:num_ele);
        elseif pos_neg_peak(p)==-1
            [~,tps] = sort(mean(avgdata.avg(:,tim{i}(1):tim{i}(2)),2),'ascend');
            tp{i} = tps(1:num_ele);
        end
        
        f3(p,i)=figure('units','normalized','outerposition',[0 0 1 1]);
        ts=avgdata.time;
        ptimes = 1:length(avgdata.time);
        hold all
        st=avgdata.time(tim{i}(1):tim{i}(2));
        fill([st, fliplr(st)], [ones(1,length(st))*100, fliplr(ones(1,length(st))*-100)], ...
        [0.5 0.5 0.5], 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
        %plot([lat_ms lat_ms],[-10 10],'k');
        lowest=Inf;
        highest=-Inf;
        for ii = 1:length(tldata)
            peakdata = tldata{ii}.avg;
            vardata = tldata{ii}.var;
            ERP = mean(peakdata(tp{i},:),1);
            VAR = mean(sqrt(vardata(tp{i},:)),1);
            %nsub = length(tldata(ii,:));
            %SEM = VAR/sqrt(nsub);               % Standard Error
            %tscore = -tinv(0.025,nsub-1);      % T-Score
            %CI = tscore*SEM;                      % Confidence Intervals
            upper = ERP(ptimes)+VAR(ptimes);
            lower = ERP(ptimes)-VAR(ptimes);
            lowest = min(min(lower),lowest);
            highest = max(max(upper),highest);
            fill([ts, fliplr(ts)], [(upper), fliplr((lower))], ...
            col(ii,:), 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
            plot(ts,ERP(ptimes),'color',col(ii,:)); 
        end
        ylabel('Amplitude, uV') % label for y axis
        xlabel('Time (ms)') % label for x axis
        set(gca,'FontSize',15);
        set(findall(gcf,'type','text'),'FontSize',15);
        ylim([lowest highest]);
        hold off
    end

    f4(p) = figure('units','normalized','outerposition',[0 0 0.1*length(tldata) 0.2*length(S.ploterp.times)]);
    plotchans=1:length(S.(S.func).chanlocs);
    iii=0;
    for i = 1:length(S.ploterp.times)
        for ii = 1:length(tldata);
            iii = iii+1;
            peakdata = tldata{ii}.avg;
            [~,markchans] = intersect(plotchans,tp{i});
            subplot(length(tldata),length(S.ploterp.times),iii); 
            topoplot(mean(peakdata(:,tim{i}(1):tim{i}(2)),2), S.(S.func).chanlocs(S.(S.func).inclchan),'maplimits',[-4 4],'electrodes','on','plotchans',plotchans,'emarker2',{markchans,'o','w',3,1}); colorbar
        end
    end
end
