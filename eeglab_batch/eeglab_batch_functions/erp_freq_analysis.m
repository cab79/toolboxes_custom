function S=erp_freq_analysis(S)
%% ERP, Frequency, Time-Frequency and Coherence analysis
% To be run after data can be cleaned 
% Before running this, check data quality in EEGLAB
% Requires Fieldtrip to be in the Matlab path

try 
    if strcmp(S.erp.select.datatype,'ERP')
        S.func = 'erp';
        S.(S.func).nf = 1;
    end
catch
    if strcmp(S.(S.func).select.datatype,'TF')
        S.func = 'tf';
        S.(S.func).nf = length(S.(S.func).select.freq);
    end
end
if strcmp(S.(S.func).select.datatype,'ERP') && S.(S.func).CSD.apply 
    load(S.(S.func).CSD.montage); % path to CSD montage (only needed if using CSD)
end 

if S.(S.func).epoch.combinemarkers
    S.(S.func).nMarkType = 1;
else
    S.(S.func).nMarkType = length(S.(S.func).epoch.markers);
end

% % load directory
% if ~isfield(S.(S.func),'loaddir')
%     S.(S.func).loaddir = fullfile(S.path.prep,S.(S.func).load.suffix{:});
% end
% 
% % previously processed files
% if isfield(S.(S.func),'filelist')
%     if S.(S.func).overwrite==0
%         prev_filelist = S.(S.func).filelist;
%         prev_dirlist = S.(S.func).dirlist;
%         prev_subj_pdat_idx = S.(S.func).subj_pdat_idx;
%         prev_designtab = S.(S.func).designtab;
%     end
%     S.(S.func) = rmfield(S.(S.func),'dirlist');
%     S.(S.func) = rmfield(S.(S.func),'subj_pdat_idx');
%     S.(S.func) = rmfield(S.(S.func),'designtab');
% end
% 
% % GET FILE LIST
% S.path.file = S.(S.func).loaddir;
% if ~isempty(S.(S.func).load.suffix{:})
%     S = getfilelist(S,S.(S.func).load.suffix);
% else
%     S = getfilelist(S);
% end
% 
% % select which to process
% if S.(S.func).overwrite==0 && exist('prev_filelist','var')
%     idx_remove = ismember(S.(S.func).filelist,prev_filelist);
%     S.(S.func).filelist(idx_remove)=[];
%     S.(S.func).dirlist(idx_remove)=[];
%     S.(S.func).subj_pdat_idx(idx_remove)=[];
%     S.(S.func).designtab(idx_remove,:)=[];
% end

S = filehandler(S,'start');

% change to the input directory
eval(sprintf('%s', ['cd(''' S.path.file ''')']));

% report if there are no such files
if isempty(S.(S.func).filelist)
    error('No files found!\n');
end

% run though all files in a loop
for f = S.(S.func).startfile:length(S.(S.func).filelist)

    file = S.(S.func).filelist{f};
    EEG = pop_loadset('filename',file,'filepath',S.path.file);

    % collect summary data
    S.(S.func).outtable.file{S.fn+f} = file;
    
    if isempty(EEG.epoch)
        continue
    end

    % SAVE CHANLOCS FOR PLOTTING LATER (needed for EEGLAB's topoplot function)
    %if f==1
    %    if ~exist(fullfile(S.path.prep,'chanlocs.mat'),'file')
    %        chanlocs = EEG.chanlocs;
    %        save(fullfile(S.path.prep,'chanlocs.mat'),'chanlocs')
    %    end
    %end
    
    % add in any missing chans, fill with nan
    try
        load(fullfile(S.path.main,'chanlocs.mat'));
    catch
        chanlocs = EEG.chanlocs;
        save(fullfile(S.path.main,'chanlocs.mat'),'chanlocs')
    end
    labs = {chanlocs.labels};
    actlabs = {EEG.chanlocs.labels};
    missing = find(~ismember(labs,actlabs));
    if ~isempty(missing)
        % first ensure chanlocs and EEG.chanlocs are similar structures
        fdnames=fieldnames(EEG.chanlocs);
        fm=~ismember(fdnames,fieldnames(chanlocs));
        EEG.chanlocs=rmfield(EEG.chanlocs,fdnames(fm));
        % then replace missing channels
        for m = 1:length(missing)
            EEG.chanlocs(missing(m)+1:end+1) = EEG.chanlocs(missing(m):end);
            EEG.chanlocs(missing(m)) = chanlocs(missing(m));
            EEG.data(missing(m)+1:end+1,:,:) = EEG.data(missing(m):end,:,:);
            EEG.data(missing(m),:,:) = nan;
        end
        EEG.nbchan = length(EEG.chanlocs);
        if EEG.nbchan~=length(labs)
            error('not enough chans')
        end
    end

    % Remove ECG if not done already
    ecgchan = strcmp({EEG.chanlocs.labels},'ECG');
    if any(ecgchan)
        EEG = pop_select(EEG,'nochannel',find(ecgchan));
    end
    
    % downsample
    if S.(S.func).epoch.downsample
        EEG = pop_resample(EEG, S.(S.func).epoch.downsample);
    end
    
    % LINEAR DETREND
    if S.(S.func).epoch.detrend
        tim = dsearchn(EEG.times',S.(S.func).epoch.detrend'); 
        for i = 1:EEG.trials, EEG.data(:,tim(1):tim(2),i) = detrend(EEG.data(:,tim(1):tim(2),i)')'; end
    end

    % remove baseline
    if S.(S.func).epoch.rmbase
        EEG = pop_rmbase( EEG, [S.(S.func).epoch.basewin(1)*1000 S.(S.func).epoch.basewin(2)*1000]);
    end
    
    dsize = size(EEG.data);
    S.(S.func).intimes=EEG.times;
    %SINGDAT = cell(1);
    
    % prepare for saving
    [pth nme ext] = fileparts(file); 
    sname_ext = [S.(S.func).select.datatype '.mat'];
    sname = strrep(file,[S.(S.func).load.suffix{:} '.' S.(S.func).fname.ext{:}],sname_ext);
    
    % This is study-specific code for data recorded from the Cambridge EGI system
    if any(find(strcmp('STIM',EEG.epoch(1).eventtype)))
        [conds, tnums, fnums, bnums] = get_markers(EEG);
    end
    
%     % This is study-specific code for MoNoPly study where there are
%     % duplicate events - only use the second event in each epoch
%     if S.(S.func).epoch.secondeventonly
%         for e = 1:length(EEG.epoch)
%             if length(EEG.epoch(e).eventtype)>1
%                 EEG.epoch(e).event = EEG.epoch(e).event(end);
%                 EEG.epoch(e).eventtype = EEG.epoch(e).eventtype(end);
%                 EEG.epoch(e).eventlatency = EEG.epoch(e).eventlatency(end);
%                 EEG.epoch(e).eventcode = EEG.epoch(e).eventcode(end);
%                 EEG.epoch(e).eventduration = EEG.epoch(e).eventduration(end);
%                 EEG.epoch(e).eventchannel = EEG.epoch(e).eventchannel(end);
%                 EEG.epoch(e).eventbvtime = EEG.epoch(e).eventbvtime(end);
%                 EEG.epoch(e).eventurevent = EEG.epoch(e).eventurevent(end);
%             end
%         end
%     end

    % flip channels right to left (e.g. CORE study)
    if isfield(S.(S.func),'flipchans') && ~isempty(S.(S.func).flipchans)
        flipidx = find(ismember(fnums,S.(S.func).flipchans));
        EEG.data(:,:,flipidx) = flipchan(EEG.data(:,:,flipidx),EEG.chanlocs);
    end
    
    % switch to selected analysis
    switch S.(S.func).select.datatype
        case 'ERP'
            EEGall=EEG;
            tldata={};
            maxlat_all=nan(1,S.(S.func).nMarkType);
            for mt = 1:S.(S.func).nMarkType
                if S.(S.func).epoch.combinemarkers
                    if any(find(strcmp('STIM',EEG.epoch(1).eventtype))) % EGI
                        % not complete
                    else
                        EEG = pop_selectevent(EEGall,'type',S.(S.func).epoch.markers,'latency',S.(S.func).epoch.latency); 
                    end
                else
                    %if any(find(strcmp('STIM',EEG.epoch(1).eventtype))) % EGI
                        % not complete
                    %else
                        try
                            EEG = pop_selectevent(EEGall,'type',S.(S.func).epoch.markers{mt},'latency',S.(S.func).epoch.latency);
                        catch
                            continue
                        end
                    %end
                end

                % This is study-specific code for MoNoPly study where there are
                % duplicate events - only use the second event in each epoch
                events = [];
                if S.(S.func).epoch.secondeventonly
                    for e = 1:length(EEG.epoch)
                        if length(EEG.epoch(e).eventtype)>1
                            EEG.epoch(e).event = EEG.epoch(e).event(end);
                            EEG.epoch(e).eventtype = EEG.epoch(e).eventtype(end);
                            EEG.epoch(e).eventlatency = EEG.epoch(e).eventlatency(end);
                            EEG.epoch(e).eventcode = EEG.epoch(e).eventcode(end);
                            EEG.epoch(e).eventduration = EEG.epoch(e).eventduration(end);
                            EEG.epoch(e).eventchannel = EEG.epoch(e).eventchannel(end);
                            if ~isempty(EEG.epoch(e).eventbvtime); EEG.epoch(e).eventbvtime = EEG.epoch(e).eventbvtime(end); end
                            EEG.epoch(e).eventurevent = EEG.epoch(e).eventurevent(end);
                        end
                        if iscell(EEG.epoch(e).eventurevent) events = [events EEG.epoch(e).eventurevent{:}];  
                        else events = [events EEG.epoch(e).eventurevent]; 
                        end
                    end
                    % update events
                    
                end
            
                % convert to Fieldtrip
                FTEEG=convertoft(EEG);
                FTEEG.event=EEG.event;

                % remove incorrect events
                if ~isempty(events)
                    FTEEG.event = FTEEG.event(ismember([FTEEG.event(:).urevent],events));
                end
                
                % ERP: timelocked (evoked) data
                if ~isfield(S.(S.func).epoch,'keeptrials')
                    S.(S.func).epoch.keeptrials = 'no';
                end
                cfg = [];
                cfg.keeptrials = S.(S.func).epoch.keeptrials;
                cfg.feedback = 'textbar';
                tldata{mt} = ft_timelockanalysis(cfg,FTEEG);
                tldata{mt}.ntrials = length(FTEEG.trial);
                tldata{mt}.median = squeeze(median(tldata{mt}.trial,1));

                % apply CSD
                if S.(S.func).CSD.apply
                    [ntrial,nchan,nsamp] = size(tldata{mt}.trial);
                    trials = reshape(permute(tldata{mt}.trial,[2 3 1]),nchan,nsamp*ntrial);
                    trials=CSD(trials,G,H); 
                    trials = permute(reshape(trials,nchan,nsamp,ntrial),[3 1 2]);
                    tldata{mt}.trial = trials;
                end
                %tldata{mt} = ft_timelockanalysis(cfg,tldata{mt});
                %erp = squeeze(double(mean(tldata.trial,1)));

                % collect summary data: number of trials
                mname = S.(S.func).epoch.markers{mt};
                mname(isspace(mname)) = [];
                S.(S.func).outtable.(['Ntrials_' mname]){S.fn+f} = tldata{mt}.ntrials;

                % collect summary data: SNR
                basewin=dsearchn(EEG.times',1000*S.erp.SNR.basewin');
                sigwin=dsearchn(EEG.times',1000*S.erp.SNR.signalwin');
                base=EEG.data(:,basewin(1):basewin(2),:);
                sig=EEG.data(:,sigwin(1):sigwin(2),:);
                S.(S.func).outtable.(['signal_' mname]){S.fn+f} = mean(rms(std(sig,0,1),2),3); % mean over trials of the RMS over time of the GFP over channels
                S.(S.func).outtable.(['noise_' mname]){S.fn+f} = mean(rms(std(base,0,1),2),3); % mean over trials of the RMS over time of the GFP over channels
                S.(S.func).outtable.(['SNR_' mname]){S.fn+f} = mean(rms(std(sig,0,1),2),3) / mean(rms(std(base,0,1),2),3); % mean over trials of the RMS over time of the GFP over channels
                
                % SME (see ERPLAB)
                tw=S.erp.SNR.timewin; % time window to average over
                [~,maxlat] = max(mean(std(sig,[],1),3));
                mlat = sigwin(1)+maxlat-1; % adjust to sig onset
                meanamp = squeeze(mean(EEG.data(:,floor(mlat-EEG.srate/(tw*1000)):ceil(mlat+EEG.srate/(tw*1000)),:),2));
                SME = std(meanamp,[],2)./sqrt(size(meanamp,2));
                % mean over chans
                S.(S.func).outtable.(['SME_' mname]){S.fn+f} = mean(SME);
                % latency
                S.(S.func).outtable.(['maxlat_' mname]){S.fn+f} = EEG.times(mlat);
                maxlat_all(mt)=EEG.times(mlat);

                % topo correlation
                topo_rho=[];
                topo_fracsig=[];
                for sm = 1:size(sig,2)
                    [rho, pval] = corr(squeeze(sig(:,sm,:)),'type','Spearman');
                    upperT = triu(ones(size(rho)));
                    topo_rho(sm)= median(rho(upperT==1));
                    topo_fracsig(sm)= sum((pval(upperT==1)<0.05))/numel(pval);
                end
                S.(S.func).outtable.(['signal_topocorr_' mname]){S.fn+f} = mean(topo_rho);
                S.(S.func).outtable.(['fracsig_topocorr_' mname]){S.fn+f} = mean(topo_fracsig);

                % collect summary data: signal difference from zero (t-test over trials)
                if size(sig,3)>1
                    sig2d = reshape(sig,[],size(sig,3));
                    hs=[];zs=[];
                    for dp = 1:size(sig2d,1)
                        if all(isnan(sig2d(dp,:)))
                            hs(dp) = nan;
                            zs(dp) = nan;
                        else
                            [~,hs(dp),stats] = signrank(sig2d(dp,:),0,'method','approximate');
                            zs(dp) = abs(stats.zval);
                        end
                    end
                    S.(S.func).outtable.(['signrank_meanabsZ_' mname]){S.fn+f} = nanmean(zs); 
                    S.(S.func).outtable.(['signrank_fracSig_' mname]){S.fn+f} = nanmean(hs); 
                else
                    S.(S.func).outtable.(['signrank_meanabsZ_' mname]){S.fn+f} = nan; 
                    S.(S.func).outtable.(['signrank_fracSig_' mname]){S.fn+f} = nan; 
                end

            end
            S.(S.func).outtable.nMarkerTypes{S.fn+f} = sum(~cellfun(@isempty,tldata));
            
            %summarise over events
            totalNtrials = 0;
            FracSigMean = nan(1,length(tldata));
            SNRmean = nan(1,length(tldata));
            Signalmean = nan(1,length(tldata));
            Noisemean = nan(1,length(tldata));
            topomean = nan(1,length(tldata));
            topofracsig = nan(1,length(tldata));
            meanSME = nan(1,length(tldata));
            for mt = 1:length(tldata)
                mname = S.(S.func).epoch.markers{mt};
                mname(isspace(mname)) = [];
                if isfield(tldata{mt},'ntrials')
                    totalNtrials = totalNtrials+tldata{mt}.ntrials;
                    FracSigMean(1,mt) = S.(S.func).outtable.(['signrank_fracSig_' mname]){S.fn+f};
                    SNRmean(1,mt) = S.(S.func).outtable.(['SNR_' mname]){S.fn+f};
                    Signalmean(1,mt) = S.(S.func).outtable.(['signal_' mname]){S.fn+f};
                    Noisemean(1,mt) = S.(S.func).outtable.(['noise_' mname]){S.fn+f};
                    topomean(1,mt) = S.(S.func).outtable.(['signal_topocorr_' mname]){S.fn+f};
                    topofracsig(1,mt) = S.(S.func).outtable.(['fracsig_topocorr_' mname]){S.fn+f};
                    meanSME(1,mt) = S.(S.func).outtable.(['SME_' mname]){S.fn+f};
                end
            end
            S.(S.func).outtable.(['Ntrials_Total']){S.fn+f} = totalNtrials;
            S.(S.func).outtable.(['SignRank_FracSigMean']){S.fn+f} = nanmean(FracSigMean);
            S.(S.func).outtable.(['SNRmean']){S.fn+f} = nanmean(SNRmean);
            S.(S.func).outtable.(['Signalmean']){S.fn+f} = nanmean(Signalmean);
            S.(S.func).outtable.(['Noisemean']){S.fn+f} = nanmean(Noisemean);
            S.(S.func).outtable.(['signal_topocorr_mean']){S.fn+f} = nanmean(topomean);
            S.(S.func).outtable.(['fracsig_topocorr_mean']){S.fn+f} = nanmean(topofracsig);
            S.(S.func).outtable.(['SME_mean']){S.fn+f} = nanmean(meanSME);
            S.(S.func).outtable.(['SME_maxlat']){S.fn+f} = nanstd(maxlat_all,[],2)./sqrt(sum(~isnan(maxlat_all)));
            
            % SAVE
            save(fullfile(S.path.erp,S.erp.save.dir{:},sname),'tldata');
            
        case {'Freq','TF','Coh'}
            
            % convert to Fieldtrip
            FTEEG=convertoft(EEG);
            FTEEG.event=EEG.event;
            
            if strcmp(S.(S.func).freq.type,'induced') 
                % timelocked (evoked) data
                cfg = [];
                cfg.keeptrials = 'yes';
                cfg.feedback = 'textbar';
                tldata = ft_timelockanalysis(cfg,FTEEG);
                % subtract out ERP for induced activity
                for t = 1:length(FTEEG)
                    FTEEG.trial{t} = FTEEG.trial{t}-tldata.avg; 
                end
            end
            
            % freq analysis
            fdata = FTfreqanalysis(FTEEG,S);
            fdata.time = FTEEG.time;
            act_freq = fdata.freq;
            disp(['actual frequencies: ' num2str(act_freq)])
            
            % moving averages over time
            if S.(S.func).op.mov_avg_trials
                nAvg = floor((size(fdata.powspctrm,1)/S.(S.func).op.mov_avg_trials_step)-(S.(S.func).op.mov_avg_trials/S.(S.func).op.mov_avg_trials_step))+1;
                for ma = 1:nAvg
                    st_ind = (ma-1)*S.(S.func).op.mov_avg_trials_step+1;
                    fdata.madata(:,:,ma) = squeeze(nanmean(fdata.powspctrm(st_ind:st_ind+S.(S.func).op.mov_avg_trials-1,:,:),1));
                end
            end

            % COHERENCE (WPLI)
            if strcmp(S.(S.func).select.datatype,'Coh')
                % coherence matrices
                nchan = length(fdata.label);
                matrix=zeros(size(S.(S.func).select.freq,1),nchan,nchan); 
                coh = zeros(nchan,nchan);
                cohboot = zeros(nchan,nchan);

                wpli = ft_connectivity_wpli(fdata.crsspctrm,'debias',true,'dojack',false);
                for c = 1:size(S.(S.func).select.freq,1)
                    fprintf('f %d',c);
                    [M, bstart] = min(abs(fdata.freq-freqlist(c,1)));
                    [M, bend] = min(abs(fdata.freq-freqlist(c,2)));
                    coh(:) = 0;
                    coh(logical(tril(ones(size(coh)),-1))) = max(wpli(:,bstart:bend),[],2);
                    coh = tril(coh,1)+tril(coh,1)';
                    matrix(c,:,:) = coh;
                end

                for nboot = 1:S.(S.func).freq.bootrep
                    fprintf('nboot %d',nboot);
                    for tr = 1:length(fdata.trial)
                        for ele = 1:size(fdata.trial{tr},1);
                            trial = fdata.trial{tr};
                            surrEEG=[];
                            surrEEG = [phaseran(squeeze(trial(tr,:)'),1); 0]';
                            %surrEEG = [phaseran(EEG.trial{1}',1); zeros(1, length(EEG.label))]';
                            trial(tr,:) = surrEEG;
                        end
                        fdata.trial{tr} = trial;
                    end

                    fdata = FTfreqanalysis(fdata,S);
                    wpli_boot = ft_connectivity_wpli(fdata.crsspctrm,'debias',true,'dojack',false);

                    for c = 1:size(freqlist,1)
                        [M, bstart] = min(abs(fdata.freq-freqlist(c,1)));
                        [M, bend] = min(abs(fdata.freq-freqlist(c,2)));

                        cohboot(:) = 0;
                        cohboot(logical(tril(ones(size(cohboot)),-1))) = max(wpli_boot(:,bstart:bend),[],2);
                        cohboot = tril(cohboot,1)+tril(cohboot,1)';
                        bootmat(c,:,:,nboot) = cohboot;
                    end
                    clear wpli_boot
                end
            end
            freqs = S.(S.func).select.freq;

            if strcmp(S.(S.func).select.datatype,'TF') && ~isempty(S.(S.func).freq.basenorm)
                % remove baseline
                cfg.baseline     = S.(S.func).epoch.basewin;
                cfg.baselinetype = S.(S.func).freq.basenorm; %'absolute', 'relative', 'relchange', 'normchange' or 'db' (default = 'absolute')
                cfg.parameter    = 'powspctrm';
                fdata = ft_freqbaseline_CAB(cfg, fdata);
            end
            
            % separate marker types
            if ~S.(S.func).epoch.combinemarkers
                % get marker occuring at zero latency
                events={};
                for m=1:length(EEG.epoch)
                    if iscell(EEG.epoch(m).eventlatency)
                        events{m} = EEG.epoch(m).eventtype{find([EEG.epoch(m).eventlatency{:}]==0)};
                    else
                        events{m} = EEG.epoch(m).eventtype;
                    end
                end
                for mt = 1:S.(S.func).nMarkType
                    cfg.trials      = find(strcmp(events,S.(S.func).epoch.markers{mt}));
                    cfg.avgoverrpt  = 'no';
                    fdatat{mt} = ft_selectdata(cfg, fdata);
                end
            else
                fdatat = {fdata};
            end
            
            if isfield(fdata,'madata')
                fdatat{1}.madata=fdata.madata;
            end
            
            fdata=fdatat;
            %include trial info
            if exist('conds','var')
                fdata{1}.conds = conds; fdata{1}.tnums = tnums; fdata{1}.fnums = fnums; fdata{1}.bnums = bnums;
            end
            clear fdatat;
            % SAVE
            switch S.(S.func).select.datatype
                case 'Freq'
                    if ~exist(S.path.freq,'dir')
                        mkdir(S.path.freq)
                    end
                    save(fullfile(S.path.freq,S.erp.save.dir{:},sname),'fdata');
                case 'TF'
                    if ~exist(S.path.tf,'dir')
                        mkdir(S.path.tf)
                    end
                    save(fullfile(S.path.tf,S.erp.save.dir{:},sname),'fdata');
            end
                    
    end

S.(S.func).file_processed=file;
S = filehandler(S,'update');
    
end


function [conds, tnums, fnums, bnums] = get_markers(EEG,type)

conds = nan(1,length(EEG.epoch));
tnums = nan(1,length(EEG.epoch));
fnums = nan(1,length(EEG.epoch));
bnums = nan(1,length(EEG.epoch));

switch type
    case 'EGI'
        for ep = 1:length(EEG.epoch)
            stimevidx = find(strcmp('STIM',EEG.epoch(ep).eventtype));
            if ep<length(EEG.epoch); stimevidx1 = find(strcmp('STIM',EEG.epoch(ep+1).eventtype));end
            if ~isempty(stimevidx)
                stimcodes = EEG.epoch(ep).eventcodes{stimevidx(end)};
                if ~any(strcmp('CNUM',stimcodes(:,1)))
                    error('change CNUM to FNUM to analyse conditions');
                end
                conds(1,ep) = stimcodes{strcmp('CNUM',stimcodes(:,1)),2};
                tnums(1,ep) = stimcodes{strcmp('TNUM',stimcodes(:,1)),2};
                fnums(1,ep) = stimcodes{strcmp('FNUM',stimcodes(:,1)),2};
                bnums(1,ep) = stimcodes{strcmp('BNUM',stimcodes(:,1)),2};

            else
                if length(stimevidx1)==2
                    stimcodes = EEG.epoch(ep+1).eventcodes{stimevidx1(1)};
                    conds(1,ep) = stimcodes{strcmp('CNUM',stimcodes(:,1)),2};
                    tnums(1,ep) = stimcodes{strcmp('TNUM',stimcodes(:,1)),2};
                    fnums(1,ep) = stimcodes{strcmp('FNUM',stimcodes(:,1)),2};
                    bnums(1,ep) = stimcodes{strcmp('BNUM',stimcodes(:,1)),2};
                else
                    error(['too many / too few STIMs on trial ' num2str(ep+1)])
                end
            end
        end
    case 'BV'
        evtypes = unique(horzcat(EEG.epoch(:).eventtype));
        for ep = 1:length(EEG.epoch)
            stimevidx = find([EEG.epoch(ep).eventlatency{:}]==0);
            if ~isempty(stimevidx)
                str=EEG.epoch(ep).eventtype{stimevidx};
                conds(1,ep) = str2double(regexp(str,'.* (\d+)','tokens','once'));
                tnums(1,ep) = EEG.epoch(ep).eventurevent{stimevidx};
            end
        end
end