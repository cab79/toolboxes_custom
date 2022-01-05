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

% load directory
if ~isfield(S.(S.func),'loaddir')
    S.(S.func).loaddir = fullfile(S.path.erp,S.(S.func).load.suffix{:});
end


% save directory
if ~isfield(S.(S.func),'outdir')
    S.path.outfile = fullfile(S.path.ga,S.(S.func).select.datatype);
else
    S.path.outfile = S.(S.func).outdir;
end
if ~exist(S.path.outfile,'dir')
    mkdir(S.path.outfile)
end

% GET FILE LIST
S.path.file = S.(S.func).loaddir;
S = getfilelist(S);

% select channels
S=select_chans(S);

% unique indices of S.filelist to include in separate grand averages
col_ind = find(ismember(S.(S.func).designtab.Properties.VariableNames,S.(S.func).grand_avg.parts));
designtab = S.(S.func).designtab(:,col_ind); % convert to table because unique with rows does not work on cell arrays!
[S.(S.func).designtab_temp,first_ind,file_ind]=unique(designtab,'rows','stable');
uni_ind = unique(file_ind);

% file loop
data_all = {}; % empty cell array for all subjects' tldata
S.(S.func).gadata = {}; % empty cell array for grand average data
for ga = 1:length(uni_ind)
    
    % FIND THE FILES
    S.(S.func).gafiles{ga} = S.(S.func).filelist(file_ind==uni_ind(ga));

    % initiate/clear outputs
    S.(S.func).multout = struct;
    S.(S.func).multoutlist = table;
    S.(S.func).multoutlist.subjects = S.(S.func).designtab.subjects(file_ind==uni_ind(ga));
    
    % savename
    tabcell = table2cell(designtab(first_ind(ga),:));
    sname = datestr(now,30);
    if S.ga.grand_avg.outliers 
        S.(S.func).ganame{ga} = [strjoin(tabcell,'_') '_grandavg_outliers'];
    else
        S.(S.func).ganame{ga} = [strjoin(tabcell,'_') '_grandavg'];
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

    % trim
    temp=data_all{ga}';
    data_all{ga}=vertcat(temp{:});
    %data_all{ga}=horzcat(data_all{ga}{:});
    keepind = any(~cellfun(@isempty,data_all{ga}),1);
    keepev = find(keepind);
    data_all{ga} = data_all{ga}(:,keepind);
    nev = size(data_all{ga},2); % number of events
    noemp = ~cellfun(@isempty,data_all{ga}); % some files might be missing conditions
    clear temp

    method = {'across','within'}; if S.(S.func).grand_avg.weighted; method = method{2}; else method = method{1};end
    
    % get index of files each condition came from
    Nconds = cellfun(@size,data_all{ga},'UniformOutput',0);
    fileidx=[];
    for n = 1:length(Nconds)
        fileidx(n,1:max(Nconds{n})) = n;
    end

    % multivariate outliers
    % calculate multivariate outliers (SLOW): 1=all events/subjects
    % separately, 2=each event separately over subjects, 3=concatenate events, 4=mean of events over
    % subjects
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
        if ~isempty(S.ga.grand_avg.out_chan)
            cfg.channel        = S.ga.grand_avg.out_chan; 
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

        % selecting chans and timewin
        temp = {};
        if ~isempty(S.ga.grand_avg.out_win) || ~isempty(S.ga.grand_avg.out_chan)
            for i = 1:size(data_all{ga}(:))
                if noemp(i)
                    temp{ga}{i} = ft_selectdata(cfg, data_all{ga}{i});
                end
            end
        else 
            temp{ga}=data_all{ga};
        end
        
        for n = 1:nev
            S.(S.func).fileidx = 1:length(temp{ga});
            S=MultiOutliers(S,temp{ga}(:,n)); %CHECK: NOEMP NEEDED?
            if any(isempty(fieldnames(S.(S.func).multout)))
                continue % due to not enough subjects
            end
            outdata{n} = S.(S.func).multout; % CHECK THIS
            outlist{n} = S.(S.func).multoutlist; % CHECK THIS - probably breaks the code later in "if S.ga.grand_avg.outliers..."
        end

%         % grand average, this time keeping all chans and timewin
%         cfg.channel        = 'all';
%         cfg.latency        = 'all';
%         for i = 1:size(data_all{ga},1)
%             temp{ga}{i,1} = ft_timelockgrandaverage_cab(cfg, data_all{ga}(i,noemp(i,:)));
%         end
%         data_all{ga}=temp{ga};
%         noemp=any(noemp,2);

    elseif S.ga.grand_avg.outliers==3
        if ~isempty(S.ga.grand_avg.out_chan)
            cfg.channel        = S.ga.grand_avg.out_chan; 
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

        % selecting chans and timewin
        temp = {};
        if ~isempty(S.ga.grand_avg.out_win) || ~isempty(S.ga.grand_avg.out_chan)
            for i = 1:size(data_all{ga}(:))
                if noemp(i)
                    temp{ga}{i} = ft_selectdata(cfg, data_all{ga}{i});
                end
            end
        else 
            temp{ga}=data_all{ga};
        end

        % conatenate over event types
        temp{ga} = reshape(temp{ga},[],1);
        S.(S.func).noemp_vector = reshape(noemp,[],1);

        % file index matches to concatenated events
        S.(S.func).fileidx = repmat(fileidx,size(data_all{ga},2),1);
        S.(S.func).eventidx = reshape(repmat(keepev,size(data_all{ga},1),1),{},1);
        %remove empty
        S.(S.func).fileidx = S.(S.func).fileidx(S.(S.func).noemp_vector);
        S.(S.func).eventidx = S.(S.func).eventidx(S.(S.func).noemp_vector);
        
        %S.(S.func).fileidx = 1:length(temp{ga});
        S=MultiOutliers(S,temp{ga}(S.(S.func).noemp_vector)); %CHECK: NOEMP NEEDED?
        if any(isempty(fieldnames(S.(S.func).multout)))
            continue % due to not enough subjects/events
        end
        outdata = S.(S.func).multout; 
        outdata.fileidx = S.(S.func).fileidx;
        outdata.eventidx = S.(S.func).eventidx;
        outlist = S.(S.func).multoutlist; 

        % avg trials per subject and event type, keeping all timewin and chan
        %if ~isempty(S.ga.grand_avg.out_win) || ~isempty(S.ga.grand_avg.out_chan)
        clear temp;
        cfg.channel        = 'all';
        cfg.latency        = 'all';
        for i = 1:size(data_all{ga},1)
            for n = 1:nev
                if noemp(i,n)
                    temp{ga}{i,n} = ft_timelockgrandaverage_cab(cfg, data_all{ga}{i,n});
                end
            end
        end
        %end
        data_all{ga}=temp{ga};
        
    elseif S.ga.grand_avg.outliers==4
        if ~isempty(S.ga.grand_avg.out_chan)
            cfg.channel        = S.ga.grand_avg.out_chan; 
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

        % average over all event types for each subject after selecting
        % chans and timewin
        for i = 1:size(data_all{ga},1)
            temp{ga}{i,1} = ft_timelockgrandaverage_cab(cfg, data_all{ga}(i,noemp(i,:)));
        end
        S.(S.func).fileidx = 1:length(temp{ga});
        S.(S.func).eventidx = reshape(repmat(keepev,size(data_all{ga},1),1),{},1);
        S=MultiOutliers(S,temp{ga});
        if any(isempty(fieldnames(S.(S.func).multout)))
            continue % due to not enough subjects
        end
        outdata = S.(S.func).multout;
        outdata.fileidx = S.(S.func).fileidx;
        outdata.eventidx = S.(S.func).eventidx;
        outlist = S.(S.func).multoutlist;

        % reduce data_all to a single event, this time keeping all chans
        % and timewin
        nev=1;
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
        if S.ga.grand_avg.reject_threshold
            rej_thresh = S.ga.grand_avg.reject_threshold;
        else
            rej_thresh = log10(outdata.chi_crt(1,5)); % 0.01 significance level
        end
        rej_thresh = {num2str(rej_thresh)};%inputdlg({'Select RD threshold from chart'},'',1,{num2str(rej_thresh)});
        
        rej={};
        rejsub={};
        if S.ga.grand_avg.reject_subjects
            rejsub = outlist.subjects(outlist.RDsub>=str2double(rej_thresh{:}));
            data_all_rej{ga} = data_all{ga}(ismember(outlist.subjects,rejsub),:);
            data_all_acc{ga} = data_all{ga}(~ismember(outlist.subjects,rejsub),:);
            noemp_noout = noemp(ismember(outlist.subjects,rejsub),:);
            noemp_out = noemp(~ismember(outlist.subjects,rejsub),:);
        else
            for i = 1:nev
                rej{i} = outlist.subjects(outlist.(['RD_event' num2str(keepev(i))])>=str2double(rej_thresh{:}));
                data_all_rej{ga}{i} = data_all{ga}(ismember(outlist.subjects,rej{i}),i);
                data_all_acc{ga}{i} = data_all{ga}(~ismember(outlist.subjects,rej{i}),i);
                noemp_noout{i} = noemp(ismember(outlist.subjects,rej{i}),i);
                noemp_out{i} = noemp(~ismember(outlist.subjects,rej{i}),i);
            end

        end

        save(fullfile(S.path.outfile,[S.(S.func).ganame{ga} '_outdata.mat']),'outdata','outlist','rejsub','rej','rej_thresh');
    end
    
    % for each event
    for n = 1:nev
        switch S.(S.func).select.datatype
            case {'ERP'}
                % grand average using Fieldtrip - includes all subjects

                cfg.channel        = data_all{ga}{1,1}.label(S.(S.func).inclchan);
                cfg.latency        = 'all';
                cfg.keepindividual = 'no';
                cfg.normalizevar   = 'N-1';
                cfg.method         = method;
                cfg.parameter      = 'avg';
                % includes all subjects
                S.(S.func).gadata{ga}.events{n} = ft_timelockgrandaverage_cab(cfg, data_all{ga}(noemp(:,n),n)); % for each event type
                if n==1 % over all events
                    S.(S.func).gadata{ga}.gavg = ft_timelockgrandaverage_cab(cfg, data_all{ga}{noemp});
                    S.(S.func).gadata{ga}.avg = S.(S.func).gadata{ga}.gavg.avg;
                end
                % separate GA for included and rejected subjects
                if S.ga.grand_avg.outliers 
                    if ~isempty(noemp_noout)
                        S.(S.func).gadata{ga}.events_rej{n} = ft_timelockgrandaverage_cab(cfg, data_all_rej{ga}{n}(noemp_noout{n}));
                    end
                    if ~isempty(noemp_out)
                        S.(S.func).gadata{ga}.events_acc{n} = ft_timelockgrandaverage_cab(cfg, data_all_acc{ga}{n}(noemp_out{n}));
                    end
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
                    if ~isempty(noemp_noout)
                        S.(S.func).gadata{ga}.events_rej{n} = ft_freqgrandaverage_cab(cfg, data_all_rej{ga}{n}(noemp_noout{n}));
                    end
                    if ~isempty(noemp_out)
                        S.(S.func).gadata{ga}.events_acc{n} = ft_freqgrandaverage_cab(cfg, data_all_acc{ga}{n}(noemp_out{n}));
                    end
                end
        end
    end
    
    % save
    gadata = S.(S.func).gadata{ga}; 
    save(fullfile(S.path.outfile,S.(S.func).ganame{ga}),'gadata'); 
    
    if S.ga.grand_avg.similarity_score
        % select chans/times for all data (including rejected subjects)
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
        sim_all{ga}=data_all{ga};
        if ~isempty(S.ga.grand_avg.out_win) || ~isempty(S.ga.grand_avg.out_chan) || ~isempty(S.ga.grand_avg.sim_win) || ~isempty(S.ga.grand_avg.sim_chan)
            for i = 1:size(data_all{ga}(:))
                if noemp(i)
                    sim_all{ga}{i} = ft_selectdata(cfg, data_all{ga}{i});
                end
            end
        end
        
        % re-do GA after selecting chans/timewin specifically for
        % similarity scores
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

        % per event
        for n = 1:size(sim_score.sim,2)
            outdata.sim_score.(['sim_event' num2str(keepev(n))])=sim_score.sim(:,n);
            outdata.sim_score.(['inv_event' num2str(keepev(n))])=sim_score.inv(:,n);
            outdata.(['xcorr_score_event' num2str(keepev(n))])=xcorr_score(:,n);
        end
        % means
        outdata.sim_score.sim_mean=mean(sim_score.sim,2);
        outdata.sim_score.inv_mean=mean(sim_score.inv,2);
        outdata.xcorr_score_mean=mean(xcorr_score,2);
        save(fullfile(S.path.outfile,[S.(S.func).ganame{ga} '_outdata.mat']),'outdata','outlist','rejsub','rej','rej_thresh');

        % excel table
        % per event
        for n = 1:size(sim_score.sim,2)
            outlist.(['sim_event' num2str(keepev(n))])=sim_score.sim(:,n);
            outlist.(['inv_event' num2str(keepev(n))])=sim_score.inv(:,n);
            outlist.(['xcorr_score_event' num2str(keepev(n))])=xcorr_score(:,n);
            outlist.(['outlier_event' num2str(keepev(n))]) = ismember(outlist.subjects, rej{n});  
        end
        % means
        outlist.sim_score_mean=mean(sim_score.sim,2);
        outlist.inv_sim_score_mean=mean(sim_score.inv,2);
        outlist.xcorr_score_mean=mean(xcorr_score,2);
        outlist.outlier_anyevent = ismember(outlist.subjects, vertcat(rej{:}));  
        writetable(outlist,fullfile(S.path.outfile,[S.(S.func).ganame{ga} '_summary.xlsx']));
    end

end

