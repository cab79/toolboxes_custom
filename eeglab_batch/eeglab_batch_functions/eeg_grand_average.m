function S=eeg_grand_average(S)

dbstop if error

S.func = 'ga';
switch S.(S.func).select.datatype
    case 'ERP'
        S.path.file = S.path.erp;
    case 'TF'
        S.path.file = S.path.tf;
    case 'Freq'
        S.path.file = S.path.freq;
end

% GET FILE LIST
S = getfilelist(S);

% select channels
S=select_chans(S);

% unique indices of S.filelist to include in separate grand averages
col_ind = find(ismember(S.(S.func).designmat(1,:),S.(S.func).grand_avg.parts));
designtab = cell2table(S.(S.func).designmat(2:end,col_ind)); % convert to table because unique with rows does not work on cell arrays!
[~,first_ind,file_ind]=unique(designtab,'rows','stable');
uni_ind = unique(file_ind);

% file loop
data_all = {}; % empty cell array for all subjects' tldata
S.(S.func).gadata = {}; % empty cell array for grand average data
for ga = 1:length(uni_ind)
    
    % FIND THE FILES
    S.(S.func).gafiles{ga} = S.(S.func).filelist(file_ind==uni_ind(ga));
    
    % savename
    tabcell = table2cell(designtab(first_ind(ga),:));
    sname = datestr(now,30);
    if S.ga.grand_avg.outliers 
        S.(S.func).ganame{ga} = [sname '_' strjoin(tabcell,'') 'grandavg_outliers'];
    else
        S.(S.func).ganame{ga} = [sname '_' strjoin(tabcell,'') 'grandavg'];
    end
    switch S.(S.func).select.datatype
        case {'TF','Freq'}
            fr_name = '_freq';
            for fr = 1:length(S.(S.func).select.freqs)
                fr_name = [fr_name '_' num2str(S.(S.func).select.freqs(fr))];
            end
            S.(S.func).ganame{ga} = [S.(S.func).ganame{ga} fr_name];
    end
    
    % function to select events and frequencies
    S.(S.func).files = S.(S.func).gafiles{ga};
    S = select_eventsfreqs(S);
    data_all{ga} = S.(S.func).data;
    
    method = {'across','within'}; if S.(S.func).grand_avg.weighted; method = method{2}; else method = method{1};end
    
    % get index of files each condition came from
    Nconds = cellfun(@size,data_all{ga},'UniformOutput',0);
    fileidx=[];
    for n = 1:length(Nconds)
        fileidx(n,1:max(Nconds{n})) = n;
    end
    %fileidx = reshape(fileidx',1,[]);
    temp=data_all{ga}';
    data_all{ga}=vertcat(temp{:});
    nev = size(data_all{ga},2); % number of events
    clear temp
    %data_all{ga}=horzcat(data_all{ga}{:});
    
    % remove empty cells
    noemp = ~cellfun(@isempty,data_all{ga});
    %data_all{ga} = reshape(data_all{ga}(noemp),[],nev);
    S.(S.func).fileidx = fileidx(noemp);
    
    % multivariate outliers
    % calculate multivariate outliers (SLOW): 1=all events, 2=each event, 3=mean of events
    if S.ga.grand_avg.outliers==1
        if ~isempty(S.ga.grand_avg.out_win)
            cfg.latency        = S.ga.grand_avg.out_chan; 
        else
            cfg.channel        = 'all';
        end
        if ~isempty(S.ga.grand_avg.out_win)
            cfg.latency        = S.ga.grand_avg.out_win;
        else
            cfg.latency        = 'all';
        end
        cfg.keepindividual = 'no';
        cfg.normalizevar   = 'N-1';
        cfg.method         = method;
        cfg.parameter      = 'avg';
        for i = 1:numel(data_all{ga})
            if noemp(i)
                disp(['collecting data ' num2str(i) '/' num2str(numel(data_all{ga}))])
                temp{ga}{i} = ft_selectdata(cfg, data_all{ga}{i});
            end
        end
        S=MultiOutliers(S,temp{ga}(noemp));
        outdata = S.(S.func).multout;
        outlist = S.(S.func).multoutlist;
        clear temp
    elseif S.ga.grand_avg.outliers==2
        for n = 1:nev
            S=MultiOutliers(S,data_all{ga}(:,n));
            outdata{n} = S.(S.func).multout;
            outlist{n} = S.(S.func).multoutlist;
        end
    elseif S.ga.grand_avg.outliers==3
        if ~isempty(S.ga.grand_avg.out_win)
            cfg.latency        = S.ga.grand_avg.out_chan; 
        else
            cfg.channel        = 'all';
        end
        if ~isempty(S.ga.grand_avg.out_win)
            cfg.latency        = S.ga.grand_avg.out_win;
        else
            cfg.latency        = 'all';
        end
        cfg.keepindividual = 'no';
        cfg.normalizevar   = 'N-1';
        cfg.method         = method;
        cfg.parameter      = 'avg';
        for i = 1:size(data_all{ga},1)
            temp{ga}{i,1} = ft_timelockgrandaverage_cab(cfg, data_all{ga}(i,noemp(i,:)));
        end
        S.(S.func).fileidx = 1:length(temp{ga});
        S=MultiOutliers(S,temp{ga});
        outdata = S.(S.func).multout;
        outlist = S.(S.func).multoutlist;
        nev=1;
        % reduce data_all to a single event
        if ~isempty(S.ga.grand_avg.out_win) || ~isempty(S.ga.grand_avg.out_chan)
            cfg.channel        = 'all';
            cfg.latency        = 'all';
            for i = 1:size(data_all{ga},1)
                temp{ga}{i,1} = ft_timelockgrandaverage_cab(cfg, data_all{ga}(i,noemp(i,:)));
            end
        end
        data_all{ga}=temp{ga};
        noemp=any(noemp,2);
            
    end
    
    if S.ga.grand_avg.outliers 
        rej_thresh = 1.25;%log10(outdata.chi_crt(1,3));
        rej_thresh = {num2str(rej_thresh)};%inputdlg({'Select RD threshold from chart'},'',1,{num2str(rej_thresh)});
        rejsub = outlist((cell2mat(outlist(:,2))>=str2double(rej_thresh{:})),1);
        subs = S.(S.func).designmat(2:end,find(strcmp(S.(S.func).designmat(1,:),'subjects')));
        [unisubs,~,subsidx] = unique(subs,'stable');
        data_all_rej{ga} = data_all{ga}(ismember(subs,rejsub),:);
        data_all_acc{ga} = data_all{ga}(~ismember(subs,rejsub),:);
        noemp_noout = noemp(ismember(subs,rejsub),:);
        noemp_out = noemp(~ismember(subs,rejsub),:);
        save(fullfile(S.path.file,[S.(S.func).ganame{ga} '_outdata.mat']),'outdata','outlist','rejsub','rej_thresh');
    end
    
    % for each event
    for n = 1:nev
        switch S.(S.func).select.datatype
            case {'ERP'}
                % grand average using Fieldtrip
                cfg.channel        = data_all{ga}{1,1}.label(S.(S.func).inclchan);
                cfg.latency        = 'all';
                cfg.keepindividual = 'no';
                cfg.normalizevar   = 'N-1';
                cfg.method         = method;
                cfg.parameter      = 'avg';
                S.(S.func).gadata{ga}.events{n} = ft_timelockgrandaverage_cab(cfg, data_all{ga}(noemp(:,n),n));
                if n==1
                    S.(S.func).gadata{ga}.gavg = ft_timelockgrandaverage_cab(cfg, data_all{ga}{noemp});
                    S.(S.func).gadata{ga}.avg = S.(S.func).gadata{ga}.gavg.avg;
                end
                if S.ga.grand_avg.outliers 
                    S.(S.func).gadata{ga}.events_rej{n} = ft_timelockgrandaverage_cab(cfg, data_all_rej{ga}(noemp_noout(:,n),n));
                    S.(S.func).gadata{ga}.events_acc{n} = ft_timelockgrandaverage_cab(cfg, data_all_acc{ga}(noemp_out(:,n),n));
                end
            case {'Freq','TF'}
                cfg.keepindividual = 'no';
                cfg.foilim         = 'all'; %[fmin fmax] or 'all', to specify a subset of frequencies (default = 'all')
                cfg.toilim         = 'all'; % to specify a subset of latencies (default = 'all')
                cfg.channel        = data_all{ga}{1,1}.label(S.(S.func).inclchan);
                cfg.parameter      = 'powspctrm';
                S.(S.func).gadata{ga}.events{n} = ft_freqgrandaverage_cab(cfg, data_all{ga}(:,n));
                if n==1
                    S.(S.func).gadata{ga}.gavg = ft_freqgrandaverage_cab(cfg, data_all{ga}{noemp});
                    S.(S.func).gadata{ga}.avg = S.(S.func).gadata{ga}.gavg.avg;
                end
                if S.ga.grand_avg.outliers 
                    S.(S.func).gadata{ga}.events_rej{n} = ft_freqgrandaverage_cab(cfg, data_all_rej{ga}(noemp_noout(:,n),n));
                    S.(S.func).gadata{ga}.events_acc{n} = ft_freqgrandaverage_cab(cfg, data_all_acc{ga}(noemp_out(:,n),n));
                end
        end
    end
    
    gadata = S.(S.func).gadata{ga};
    save(fullfile(S.path.file,S.(S.func).ganame{ga}),'gadata'); 
    
    if S.ga.grand_avg.similarity_score
        % select chans/times
        if ~isempty(S.ga.grand_avg.sim_chan)
            cfg.channel        = S.ga.grand_avg.sim_chan; 
        else
            cfg.channel        = 'all';
        end
        if ~isempty(S.ga.grand_avg.sim_win)
            cfg.latency        = S.ga.grand_avg.sim_win;
        else
            cfg.latency        = 'all';
        end
        cfg.keepindividual = 'no';
        cfg.normalizevar   = 'N-1';
        cfg.method         = method;
        cfg.parameter      = 'avg';
        if ~isempty(S.ga.grand_avg.out_win) || ~isempty(S.ga.grand_avg.out_chan) || ~isempty(S.ga.grand_avg.sim_win) || ~isempty(S.ga.grand_avg.sim_chan)
            for i = 1:size(data_all{ga}(:))
                if noemp(i)
                    sim_all{ga}{i} = ft_selectdata(cfg, data_all{ga}{i});
                end
            end
        else 
            sim_all{ga}=data_all{ga};
        end
        
        % re-do GA
        if S.ga.grand_avg.sim_use_robustavg
            if ~isempty(S.ga.grand_avg.sim_chan)
                allchans=data_all{ga}{1,1}.label(S.(S.func).inclchan);
                selectchanind = find(ismember(allchans,S.ga.grand_avg.sim_chan));
            else
                selectchanind=1:size(S.(S.func).multout.mu,2);
            end
            if ~isempty(S.ga.grand_avg.sim_win) && isempty(S.ga.grand_avg.out_win)
                temp=dsearchn(data_all{ga}{1,1}.time',S.ga.grand_avg.sim_win');
                selecttimeind=temp(1):temp(2);
            else
                selecttimeind=1:size(S.(S.func).multout.mu,1);
            end
            template=S.(S.func).multout.mu(selecttimeind,selectchanind)';
        else
            if ~isempty(S.ga.grand_avg.sim_chan)
                allchans=data_all{ga}{1,1}.label(S.(S.func).inclchan);
                selectchanind = ismember(allchans,S.ga.grand_avg.sim_chan);
                cfg.channel = allchans(selectchanind); 
            else
                cfg.channel = data_all{ga}{1,1}.label(S.(S.func).inclchan); 
            end
            if ~isempty(S.ga.grand_avg.sim_win)
                cfg.latency = S.ga.grand_avg.sim_win;
            else
                cfg.latency = 'all';
            end
            for n = 1:nev
                S.(S.func).gadata{ga}.events{n} = ft_timelockgrandaverage_cab(cfg, data_all{ga}(noemp(:,n),n));
                if n==1
                    S.(S.func).gadata{ga}.gavg = ft_timelockgrandaverage_cab(cfg, data_all{ga}{noemp});
                    S.(S.func).gadata{ga}.avg = S.(S.func).gadata{ga}.gavg.avg;
                end
            end
            template=S.(S.func).gadata{ga}.avg;
        end
        
        sim_score.sim=zeros(size(data_all{ga}));
        sim_score.inv=zeros(size(data_all{ga}));
        sim_score.combined=zeros(size(data_all{ga}));
        xcorr_score=zeros(size(data_all{ga}));
        for i = 1:numel(data_all{ga})
            disp(['Calculating similarity for subject ' num2str(i) '/' num2str(numel(data_all{ga}))])
            if noemp(i)
                if ~isempty(S.ga.grand_avg.sim_chan)
                    dat=sim_all{ga}{i}.avg;
                else
                    dat=sim_all{ga}{i}.avg(S.(S.func).inclchan,:);
                end
            
                for e = 1:size(dat,1)
                    dat(e,:) = smooth(squeeze(dat(e,:)),30,'lowess');
                end
                sim = dat .* template;
                inv = (1./dat) .* (1./template);
                max_sim = template .* template;
                max_inv = (1./template) .* (1./template);
                sim_score.sim(i) = sum(sim,'all') / sum(max_sim,'all');
                sim_score.inv(i) = sum(inv,'all') / sum(max_inv,'all');
                xcorr_score(i) = sum(xcorr2(dat,template),'all');
            end
        end
        outdata.sim_score.sim=mean(sim_score.sim,2);
        outdata.sim_score.inv=mean(sim_score.inv,2);
        outdata.xcorr_score=mean(xcorr_score,2);
        save(fullfile(S.path.file,[S.(S.func).ganame{ga} '_outdata.mat']),'outdata','outlist','rejsub','rej_thresh');
    end

end

