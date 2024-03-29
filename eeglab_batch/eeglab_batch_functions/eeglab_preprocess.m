function S=eeglab_preprocess(S,part)
% PREPROCESSING FOR CONTINUOUS EEGLAB .SET FILES
S.func = 'prep';

switch part
    
    %% EPOCH
    case 'epoch'
   
    
    S = filehandler(S,'start');

    % change to the input directory
    eval(sprintf('%s', ['cd(''' S.path.prep ''')']));

    % report if there are no such files
    if isempty(S.prep.filelist)
        disp('No files found!\n');
    end

    % run though all files in a loop
    for f = S.prep.startfile:length(S.prep.filelist)
        file = S.prep.filelist{f};

        fprintf('\nProcessing %s.\n\n', file);

        % GET FILENAME PARTS
        [pth nme ext] = fileparts(file); 
        fparts = strsplit(nme,'_');

        % LOAD DATA
        EEG = pop_loadset('filename',file,'filepath',S.path.file);

        % collect summary data
        S.prep.outtable.file{S.fn+f} = file;

        % adjust latencies - needed for MNP study
        if isfield(S.prep.epoch,'latency_correction') && S.prep.epoch.latency_correction~=0
            
            if isfield(S.prep.epoch,'latency_correction_subjects') && ~isempty(S.prep.epoch.latency_correction_subjects)
                if ismember(S.(S.func).designtab.subjects{f},S.prep.epoch.latency_correction_subjects)
                    correct_latency=1;
                else
                    correct_latency=0;
                end
            else 
                correct_latency=1;
            end
            
            if correct_latency
                eventi = find(ismember({EEG.event(:).(S.prep.epoch.latency_correction_eventfield{1})},S.prep.epoch.latency_correction_eventfield{2}));
                for e = eventi
                    EEG.event(e).latency = EEG.event(e).latency + S.prep.epoch.latency_correction*EEG.srate;
                    EEG.urevent(e).latency = EEG.event(e).latency; % samples to add to marker latency to correct for delays in stimulus markers in SCIn
                end
            end
        end

        % Channel switching (a patch: unlikely to be needed)
        if isfield(S.prep.chan,'switchhead') && ~isempty(S.prep.chan.switchhead)
            load(fullfile(S.path.main,'chanlocs.mat'));
            newchanind = 32+find(ismember({chanlocs.labels},{EEG.chanlocs.labels}));
            chanlocs(64).labels = 'FCz';  
            EEG.chanlocs = chanlocs(newchanind);
        end
        
        %ADD CHANNEL LOCATIONS 
        if S.prep.chan.addloc && any(cellfun(@isempty,{EEG.chanlocs(:).theta}))
            EEG=pop_chanedit(EEG, 'lookup',S.path.locfile);
        end
        

        % SELECT TIME WINDOW TO ANALYSE
        if ~isempty(S.prep.cont.timewin)
            [~,~,i]=intersect(fparts,S.prep.select.conds);
            if ~isempty(S.prep.cont.timewin(i))
                EEG = pop_select(EEG,'time',S.prep.cont.timewin{i});
            end
        end

        % INTERPOLATE CHANNELS
        if ~isempty(S.prep.chan.interp) && S.prep.chan.interp~=0
            EEG= eeg_interp(EEG, S.prep.chan.interp, 'spherical');
            % make a record in EEG.chanlocs
            if ~isfield(EEG.chanlocs,'interp')
                interp = zeros(EEG.nbchan,1);
            else
                interp = cell2mat({EEG.chanlocs.interp});
            end
            interp(S.prep.chan.interp)=1;
            interp=num2cell(interp);
            [EEG.chanlocs.interp] = interp{:};
        end

        % RE_REFERENCE CHANNELS
        if S.prep.chan.reref == 1
            % re-reference to the common average
            EEG = pop_reref( EEG, []); 
        elseif S.prep.chan.reref==2
            %% PrepPipeline - performs robust re-referencing
            % conventional re-referencing to the common average, prior to
            % cleaning could introduce artefact to otherwise clean channels.
            % Robust averaging as done here gets around this.
            params = struct();
            params.lineFrequencies = [50, 100, 150, 200, 250];
            params.referenceChannels = 1:length(EEG.chanlocs);
            params.evaluationChannels = 1:length(EEG.chanlocs);
            params.rereferencedChannels = 1:length(EEG.chanlocs);
            params.detrendChannels = 1:length(EEG.chanlocs);
            params.lineNoiseChannels = 1:length(EEG.chanlocs);
            params.detrendType = 'high pass';
            params.detrendCutoff = 1;
            params.referenceType = 'robust';
            params.meanEstimateType = 'median';
            params.interpolationOrder = 'post-reference';
            params.keepFiltered = false;
            params.ignoreBoundaryEvents = true;
            params.removeInterpolatedChannels = true;
            [EEG, params, computationTimes] = prepPipeline(EEG, params); 
            fprintf('Computation times (seconds):\n   %s\n', ...
                getStructureString(computationTimes));
%             fprintf('Post-process\n')
%             EEG = prepPostProcess(EEG, params);
%             save(fname, 'EEG', '-mat', '-v7.3'); 
        end

        % NOTCH FILTER
        if ~isempty(S.prep.filter.notch) && iscell(S.prep.filter.notch)
            for i = 1:length(S.prep.filter.notch)
                EEG = pop_eegfiltnew(EEG,S.prep.filter.notch{i}(1),S.prep.filter.notch{i}(2),[],1); 
            end
        end

        % FILTER
        if ~isempty(S.prep.filter.incl)
            if S.prep.filter.incl(1)>0; EEG = pop_eegfiltnew( EEG, S.prep.filter.incl(1), 0, [], 0);end
            if S.prep.filter.incl(2)>0; EEG = pop_eegfiltnew( EEG, 0, S.prep.filter.incl(2), [], 0);end
        end

        % DOWNSAMPLE
        % do this after filtering if max freq to be analysed is more than half
        % of the downsample freq (Nyquist law)
        if S.prep.cont.downsample
            EEG = pop_resample(EEG, S.prep.cont.downsample);
        end

        % EPOCH
        %create epochs if markers exist
        if isfield(S.prep.epoch,'timewin_padding') && ~isempty(S.prep.epoch.timewin_padding) && iscell(S.prep.epoch.markers) % markers set by condition

            % this bit is complicated - need to pad the data for conditions
            % that have fewer data points that other conditions, if they
            % are to be combined later for rejection/ICA. Will use baseline
            % data to pad. Currently assumes the shortest data is due due
            % to including less pre-stim EEG.
            all_tw = vertcat(S.prep.epoch.timewin{:});
            max_tw = [min(all_tw(:,1)), max(all_tw(:,2))];

            % epochs same for all files
            EEG = pop_epoch( EEG, S.prep.epoch.markers, max_tw);

            % for this file
            if ~isempty(S.prep.epoch.timewinby)
                timewin = S.prep.epoch.timewin{ismember(S.prep.select.(S.prep.epoch.timewinby{:}),S.prep.designtab.(S.prep.epoch.timewinby{:}){f})};
            else
                timewin = S.prep.epoch.timewin{:};
            end
            %EEGtemp2 = pop_epoch( EEG, S.prep.epoch.markers, timewin);
        
            if timewin(1)>max_tw(1)
                %figure; plot(mean(EEG.data(1,:,:),3));
                % get pad data
                pad = S.prep.epoch.timewin_padding;
                pad_datapoint1=dsearchn(EEG.times',pad(1)*1000);
                pad_datapoint2=dsearchn(EEG.times',pad(2)*1000)-1;
                pad_data = EEG.data(:,pad_datapoint1:pad_datapoint2,:);
                padlen = size(pad_data,2);
                % back-pad
                startpoint = pad_datapoint1-padlen;
                endpoint = pad_datapoint2-padlen;
                inv_pad_data = flip(pad_data,2); % invert timepoints of data to ensure ampltiudes are aligned
                while endpoint>=1
                    pad_idx = max(1,startpoint):endpoint;
                    EEG.data(:,pad_idx,:) = inv_pad_data(:,1:length(pad_idx),:);
                    startpoint = startpoint-padlen;
                    endpoint = endpoint-padlen;
                    inv_pad_data = flip(inv_pad_data,2); % invert timepoints of data to ensure ampltiudes are aligned
                end
               % figure; plot(mean(EEG.data(1,:,:),3));
            end

        elseif ~isempty(S.prep.epoch.markers)
            % for this file
            if ~isempty(S.prep.epoch.timewinby)
                timewin = S.prep.epoch.timewin{ismember(S.prep.select.(S.prep.epoch.timewinby{:}),S.prep.designtab.(S.prep.epoch.timewinby{:}){f})};
            else
                timewin = S.prep.epoch.timewin{:};
            end
            EEG = pop_epoch( EEG, S.prep.epoch.markers, timewin);
        end

        % if markers don't exist, add markers and epoch
        if isfield(S.prep.epoch,'importmarker') && ~isempty(S.prep.epoch.importmarker)
            load(S.prep.epoch.importmarker);
            [uMark, iM1, iM2] = unique(seq.condnum','rows');
            eventind = find(ismember({EEG.event.code},'Stimulus') & ismember({EEG.event.type},S.prep.epoch.markers));
            markers = cellstr([repmat('S ',length(iM2),1), num2str(iM2)]);
            markers_present = find(~ismember(markers,{'S 16'}));
            markers_missing = find(ismember(markers,{'S 16'}));
            [EEG.event(eventind).type] = deal(markers{markers_present(1:length(eventind))});
            [EEG.event(eventind(markers_missing)).type] = deal(markers{markers_missing});
            EEG = pop_epoch( EEG, unique(markers), S.prep.epoch.timewin);
        elseif S.prep.epoch.addmarker && isempty(EEG.epoch)
            Sr = EEG.srate; % sampling rate of data
%             Ndp = Sr*(S.prep.epoch.timewin(2)-S.prep.epoch.timewin(1));% number of data points per epoch
            Ndp = Sr*S.prep.epoch.timewin(2);% number of data points post-marker
            Tdp = size(EEG.data,2);% total number of data points in file
            Mep = floor(Tdp/Ndp);% max possible number of epochs
            for i = 1:Mep
                EEG.event(1,i).type = S.prep.epoch.markers{1};
                EEG.event(1,i).latency = -Sr*S.prep.epoch.timewin(1)+(i-1)*Ndp;
                EEG.event(1,i).urevent = S.prep.epoch.markers{1};
            end
            EEG = pop_epoch( EEG, S.prep.epoch.markers(1), S.prep.epoch.timewin);
        end

        % study-specific code for CORE study, where DIN2 is use to time
        % the events, but the even information (marker type) is stored in
        % STIM
        if any(strcmp(S.prep.epoch.markers,'DIN2'))
            rm_trials = [];
            for ep = 1:length(EEG.epoch)
                stimevidx = find(strcmp('STIM',EEG.epoch(ep).eventtype));
                if ep<length(EEG.epoch); stimevidx1 = find(strcmp('STIM',EEG.epoch(ep+1).eventtype));end;
                if isempty(stimevidx)
                    if length(stimevidx1)>1
                        event = EEG.epoch(ep+1).eventurevent{1, stimevidx1(1)};
                        event_idx = find([EEG.event.urevent]==event);
                        %EEG = pop_editeventvals(EEG, 'changefield', {event_idx 'latency' eeg_point2lat(EEG.event(event_idx-1).latency,[],EEG.srate, [EEG.xmin EEG.xmax])});
                        EEG.event(event_idx).latency = EEG.event(event_idx-1).latency;
                    else
                        rm_trials = [rm_trials ep];
                    %
                    %    try 
                            
                    %    catch
                    %        error(['not enough stims in epoch ' num2str(ep+1)]);
                    %    end
                    end
                end
            end
            EEG = pop_select(EEG,'notrial',rm_trials);
            % EPOCH AGAIN AFTER STIM CORRECTION
            EEG = pop_epoch( EEG, S.prep.epoch.markers, timewin);
        end

        S.prep.outtable.initial_trials(S.fn+f) = EEG.trials;

        if isfield(S.prep.epoch,'expected_ISIs') && ~isempty(S.prep.epoch.expected_ISIs)
            if isfield(S.prep.epoch,'expected_ISIs_by') && ~isempty(S.prep.epoch.expected_ISIs_by)
                isi = S.prep.epoch.expected_ISIs(ismember(S.prep.select.(S.prep.epoch.expected_ISIs_by{:}),S.prep.designtab.(S.prep.epoch.expected_ISIs_by{:}){f}));
            else
                isi = S.prep.epoch.expected_ISIs;
            end
            [EEG,rmepochs] = rmTrailsISIvar(EEG,S.prep.epoch.expected_ISIs_errormargin,timewin,isi,S.prep.epoch.markers{:});
            S.prep.outtable.remove_Ntrials_ISIerror(S.fn+f) = length(rmepochs);
        end

        % LINEAR DETREND
        if S.prep.epoch.detrend
            for i = 1:EEG.trials, EEG.data(:,:,i) = detrend(EEG.data(:,:,i)')'; end
        end

        % REMOVE BASELINE
        if S.prep.epoch.rmbase
            if iscell(S.prep.epoch.basewin) % markers set by condition/block
                if ~isempty(S.prep.epoch.basewinby)
                    basewin = S.prep.epoch.basewin{ismember(S.prep.select.(S.prep.epoch.basewinby{:}),S.prep.designtab.(S.prep.epoch.basewinby{:}){f})};
                else
                    basewin = S.prep.epoch.basewin{:};
                end
                EEG = pop_rmbase( EEG, [basewin(1)*1000, basewin(2)*1000]);
            else
                EEG = pop_rmbase( EEG, [S.prep.epoch.basewin(1)*1000, S.prep.epoch.basewin(2)*1000]);
            end
        end

        EEG = eeg_checkset( EEG );
        
        % additional manual check of data and opportunity to reject
        % further channels and trials
        if isfield(S.prep.epoch,'manual_check') && S.prep.epoch.manual_check
            [EEG,S] = manualRejectEEG(EEG,S,0);
        end

        % SAVE
        sname = [nme '_' S.prep.save.suffix{:} '.' S.prep.fname.ext{:}];
        EEG = pop_saveset(EEG,'filename',sname,'filepath',S.prep.save.dir); 

        % update
        S.(S.func).file_processed=file;
        S = filehandler(S,'update');
    end

%     if exist('prev_filelist','var')
%         S.(S.func).filelist = [S.(S.func).filelist prev_filelist];
%         S.(S.func).dirlist = [S.(S.func).dirlist prev_dirlist];
%         S.(S.func).subj_pdat_idx = [S.(S.func).subj_pdat_idx prev_subj_pdat_idx];
%         S.(S.func).designtab = [S.(S.func).designtab; prev_designtab];
%     end

    %% COMBINE
    case 'combine'
    

    % combine data files
    if S.prep.combinefiles.on
        clear OUTEEG

        S = filehandler(S,'start');

        for s = 1:length(S.prep.select.subjects)
            if isempty(S.prep.select.sessions)
                S.prep.select.sessions = {''};
            end
            for a = 1:length(S.prep.select.sessions)

                % FIND THE FILES FOR THIS SUBJECT
                if ~isempty(S.prep.select.sessions{a})
                    subfiles = S.prep.filelist(find(not(cellfun('isempty', strfind(S.prep.filelist,S.prep.select.subjects{s}))) & not(cellfun('isempty', strfind(S.prep.filelist,S.prep.select.sessions{a})))));
                else
                    subfiles = S.prep.filelist(find(not(cellfun('isempty', strfind(S.prep.filelist,S.prep.select.subjects{s})))));
                end

                if isempty(subfiles); continue; end

                % CYCLE THROUGH EACH FILE FOR THIS SUBJECT
                for f = 1:length(subfiles)

                    % get filename parts
                    file = subfiles{f};
                    [pth nme ext] = fileparts(file); 
                    fparts = strsplit(nme,'_');

                    % create base name (without extension)
                    basename='';
                    for i = 1:length(fparts)
                        if i==1;c= '';else;c= '_';end
                        basename = [basename c fparts{i}];
                    end

                    % load and merge
                    EEG = pop_loadset('filename',file,'filepath',S.path.file);
                    [EEG.epoch(:).file] = deal(basename);
                    if exist('OUTEEG','var')
                        OUTEEG = pop_mergeset(OUTEEG, EEG);
                    else
                        OUTEEG = EEG;
                        OUTEEG.fileinfo = [];
                    end
                    OUTEEG.fileinfo = [OUTEEG.fileinfo,{EEG.epoch.file}];
                    clear EEG
                end
                EEG=OUTEEG;
                [EEG.epoch.file] = deal(EEG.fileinfo{:});
                clear OUTEEG
                EEG = eeg_checkset( EEG );

                % SAVE
                if ~isempty(S.prep.select.sessions{a})
                    sessionname = [S.prep.select.sessions{a} '_'];
                elseif length(S.prep.select.sessions)>1
                    sessionname = 'allsessions_';
                else
                    sessionname = '';
                end
                if ~isempty(S.prep.study{:})
                    sname = [S.prep.study{:} '_' S.prep.select.subjects{s} '_' sessionname S.prep.save.suffix{:} '.' S.prep.fname.ext{:}];
                else
                    sname = [S.prep.select.subjects{s} '_' sessionname S.prep.save.suffix{:} '.' S.prep.fname.ext{:}];
                end
%                 if ~exist(fullfile(S.path.prep,S.prep.save.suffix{:}),'dir')
%                     mkdir(fullfile(S.path.prep,S.prep.save.suffix{:}));
%                 end
                EEG = pop_saveset(EEG,'filename',sname,'filepath',S.prep.save.dir);

                for f = 1:length(subfiles)
                    % update
                    S.(S.func).file_processed=subfiles{f};
                    S = filehandler(S,'update');
                end

            end
        end
        

    end


    %% NOISY CHANNEL CORRECTION
    case 'chancorrect'

        S = filehandler(S,'start');
           
        for f = S.(S.func).startfile:length(S.prep.filelist)
            file = S.(S.func).filelist{f};
            disp(['loading file index ' num2str(f)])
            try
                EEG = pop_loadset('filename',file,'filepath',S.path.file);
            catch
                continue
            end

            % collect summary data
            S.prep.outtable.file{S.fn+f} = file;
            S.prep.outtable.initial_trials(S.fn+f) = EEG.trials;

            % remove selected ICA components
            if S.(S.func).epoch.ICAremove
                if ~strcmp(S.(S.func).load.suffix{:},S.(S.func).load.ica_suffix{:})

                    if ~exist(fullfile(strrep(S.path.file,S.(S.func).load.suffix{:},S.(S.func).load.ica_suffix{:}), strrep(file,S.(S.func).load.suffix{:},S.(S.func).load.ica_suffix{:})),'file')
                        continue
                    end

                    EEG_ica = pop_loadset('filename',strrep(file,S.(S.func).load.suffix{:},S.(S.func).load.ica_suffix{:}),'filepath',strrep(S.path.file,S.(S.func).load.suffix{:},S.(S.func).load.ica_suffix{:}));
                    EEG.icawinv = EEG_ica.icawinv;
                    EEG.icasphere = EEG_ica.icasphere;
                    EEG.icaweights = EEG_ica.icaweights;
                    EEG.icachansind = EEG_ica.icachansind;
                    EEG.reject = EEG_ica.reject;
                    clear EEG_ica
                end

                S.prep.outtable.ncomp(S.fn+f) = size(EEG.icaweights,1);
                S.prep.outtable.IClabel_nkeep(S.fn+f) = sum(~EEG.reject.gcompreject);
                S.prep.outtable.IClabel_frackeep(S.fn+f) = sum(~EEG.reject.gcompreject)/size(EEG.icaweights,1);

                if any(EEG.reject.gcompreject==0)
                    EEG = pop_subcomp( EEG, [], 0); 
                end
            end

            % detrend
            if S.(S.func).epoch.detrend
%                 data = EEG.data;
%                 trend = nan(size(data));
                for i = 1:EEG.trials
                    EEG.data(:,:,i) = detrend(EEG.data(:,:,i)')';
%                     trend(:,:,i) = data(:,:,i) - EEG.data(:,:,i);
                end
            end
    
            % REMOVE BASELINE
            if S.prep.epoch.rmbase
                if iscell(S.prep.epoch.basewin) % markers set by condition/block
                    if ~isempty(S.prep.epoch.basewinby)
                        basewin = S.prep.epoch.basewin{ismember(S.prep.select.(S.prep.epoch.basewinby{:}),S.prep.designtab.(S.prep.epoch.basewinby{:}){f})};
                    else
                        basewin = S.prep.epoch.basewin{:};
                    end
                    EEG = pop_rmbase( EEG, [basewin(1)*1000, basewin(2)*1000]);
                else
                    EEG = pop_rmbase( EEG, [S.prep.epoch.basewin(1)*1000, S.prep.epoch.basewin(2)*1000]);
                end
            end
    
            % select or deselect channels for cleaning
            S.(S.func).select.chans = S.(S.func).clean.chan;
            S=select_chans(S);

            
            % initial manual check of data and opportunity to reject
            % channels and trials. 
            if isfield(S.prep.clean,'initial_manual_check') && S.prep.clean.initial_manual_check
                pop_eegplot(EEG, 1, 1, 0); % view data first
                ud=get(gcf,'UserData');     %2
                ud.winlength=10;            %3
                set(gcf,'UserData',ud);     %4
                eegplot('draws',0)          %5
                h = gcf;
                waitfor(h) 
            end

%             if isfield(S.prep.clean,'man') && S.prep.clean.man
%                 close all
%                 % initial arrays
%                 rejchan = [];
%                 rejtrial = [];
%                 
%                 S.(S.func).clean.man.hidebad=0;
%                 [EEG,S] = manualRejectEEG(EEG,S,S.(S.func).clean.man.thresh); 
%                 rejchan = [rejchan;S.(S.func).clean.man.rejchan];
%                 rejtrial = [rejtrial;S.(S.func).clean.man.rejtrial];
%                 
%                 S.(S.func).clean.man.hidebad=1;
%                 [EEG,S] = manualRejectEEG(EEG,S,S.(S.func).clean.man.thresh); 
%                 rejchan = [rejchan;S.(S.func).clean.man.rejchan];
%                 rejtrial = [rejtrial;S.(S.func).clean.man.rejtrial];
% 
%                 S.prep.outtable.man_trials(S.fn+f) = numel(rejtrial);
%                 S.prep.outtable.man_chans(S.fn+f) = numel(rejchan);
%                 
%                 close all
%                 % reject bad channels
%                 EEG = reject_channels(EEG,rejchan);
%                 % reject bad trials
%                 EEG = pop_select(EEG, 'notrial', unique(rejtrial));
%             end

            % clean bursts
            if isfield(S.prep.clean,'bursts')
                EEGcont = eeg_epoch2continuous(EEG);
                [BUR,state] = clean_asr(EEGcont,...
                    S.prep.clean.bursts.BurstCriterion,...
                    S.prep.clean.bursts.WindowLength,...
                    [], ...
                    [],...
                    S.prep.clean.bursts.burst_crit_refmaxbadchns,...
                    S.prep.clean.bursts.burst_crit_reftolerances,...
                    [], ...
                    S.prep.clean.bursts.use_GPU,...
                    false,... 
                    S.prep.clean.bursts.max_mem); 
                EEG.data = reshape(BUR.data,size(EEG.data));
                frac_rej = reshape(state.frac_rej,size(EEG.data,2),size(EEG.data,3));
                frac_rej = mean(frac_rej,1);
                fig=figure('Name',file);
                plot(frac_rej); 
                title('BUR fraction components removed')
                xlabel('trials')
                savefig(fig,fullfile(S.path.prep,S.prep.save.suffix{:},strrep(file,'.set','_BURcorrection.fig'))); 
                S.prep.outtable.BUR_fracrej_median(S.fn+f) = nanmedian(frac_rej);
                S.prep.outtable.BUR_fracrej_95th(S.fn+f) = prctile(frac_rej,95);
                clear BUR EEGcont state
                close(fig)
            end


            % custom channel cleaning
            if isfield(S.prep.clean,'noisychan')

                % While loop ensures not to many channels are interpolated
                % at a time, either in full or partially
                partial_interp_allow=0;
                rejchan_all =[];
                while partial_interp_allow==0

                    rejbad = zeros(EEG.nbchan, EEG.trials);

                    % parallel
                    poolobj = gcp('nocreate'); % If no pool, do not create new one.
                    if isempty(poolobj) && S.prep.clean.noisychan.parallel_workers
                        parpool('local',S.prep.clean.noisychan.parallel_workers);
                    elseif ~isempty(poolobj) && S.prep.clean.noisychan.parallel_workers
                        if poolobj.NumWorkers ~= S.prep.clean.noisychan.parallel_workers
                            delete(poolobj)
                            parpool('local',S.prep.clean.noisychan.parallel_workers);
                        end
                    end
    
                    % correlated channel metric
                    S = correlated_channels(S,EEG,file);
                    rejbad = rejbad + S.(S.func).clean.corrchan.bad;
                    S.prep.outtable.chancorr_frac_bad(S.fn+f) = S.(S.func).clean.corrchan.frac_bad;
    
                    S = noisyflat_channel_diagnostics(S,EEG,file);
                    for oi = 1:length(S.(S.func).clean.noisychan.metric)
                        mname = S.(S.func).clean.noisychan.metric(oi).name;
                        S.prep.outtable.([mname '_var_trial'])(S.fn+f) = S.(S.func).clean.noisychan.metric(oi).var_trial;
                        S.prep.outtable.([mname '_var_chan'])(S.fn+f) = S.(S.func).clean.noisychan.metric(oi).var_chan;
                        
                        % combine outliers
                        rejbad = rejbad + S.(S.func).clean.noisychan.metric(oi).OL_upper;
                        rejbad = rejbad + S.(S.func).clean.noisychan.metric(oi).OL_lower;
                    end
                    rejbad = rejbad>=S.prep.clean.noisychan.NbadMetrics;
    
                    % interpolate whole channels if too many trials are affected (defined by S.prep.clean.noisychan.badsegs_wholechannel)
                    interp_allow=0;
                    bs_thresh = S.prep.clean.noisychan.badsegs_wholechannel;
                    while interp_allow==0 % reduce threshold until max number of bad channels reached (S.prep.clean.noisychan.badchans)
                        interp_idx=(sum(rejbad,2)/size(rejbad,2))>bs_thresh;
                        rejchan = find(interp_idx);
                        if sum(interp_idx) <= S.prep.clean.noisychan.badchans
                            interp_allow=1;
                        else
                            bs_thresh = bs_thresh+0.01;
                        end
                    end
                    if interp_allow
                        EEG = reject_channels(EEG,rejchan);
                        rejchan_all = unique([rejchan_all, rejchan']);
                        rejbad(interp_idx,:)=0;
                    end
                    rejbad_all = rejbad;
    
                    % ensure trials are ignored if they have a lot of
                    % badchans still, otherwise this could loop endlessly
                    trials_reject = sum(rejbad,1)>S.prep.clean.noisychan.badchans;
                    trials_reject = repmat(trials_reject,size(rejbad,1),1);
                    rejbad = rejbad.*(~trials_reject);

                    % update according to minimum criteria for partial
                    % correction later
                    % first, ensure channels are only selected if noisy more
                    % than X% of trials
                    chans_reject = (sum(rejbad,2)/size(rejbad,2))>=S.prep.clean.noisychan.badsegs;
                    chans_reject = repmat(chans_reject,1,size(rejbad,2));
                    rejbad = rejbad.*chans_reject;

                    if all(sum(rejbad,1)<=S.prep.clean.noisychan.badchans)
                        partial_interp_allow=1;
                        % otherwise, loop back and correct more whole
                        % channels
                    end
                end

                % summary data
                % ADD CHANNELS INTERPOLATED IN FULL
                S.prep.outtable.Nnoisy_chans_fullinterp(S.fn+f) = length(rejchan_all);
                if length(rejchan_all)>0
                    S.prep.outtable.noisy_chans_fullinterp(S.fn+f) = {{EEG.chanlocs(rejchan_all).labels}};
                end
                % ADD BAD TRIALS IGNORED FOR NOW (THEY WILL BE REMOVED LATER)
                S.prep.outtable.Nnoisy_trials_ignored_partialinterp(S.fn+f) = sum(all(trials_reject,1));
                S.prep.outtable.FracNoisy_trials_ignored_partialinterp(S.fn+f) = sum(all(trials_reject,1))/S.prep.outtable.initial_trials(S.fn+f);
                %S.prep.outtable.Noisy_trials_ignored_partialinterp(S.fn+f) = {find(~all(trials_reject,1))};
                % ADD CHANNELS/TRIALS PARTIALLY INTERPOLATED
                S.prep.outtable.Nnoisy_chans_partialinterp(S.fn+f) = sum(any(rejbad,2));
                S.prep.outtable.Nnoisy_trials_partialinterp(S.fn+f) = sum(any(rejbad,1));
                if any(rejbad,'all')
                    S.prep.outtable.noisy_chans_partialinterp(S.fn+f) = {{EEG.chanlocs(any(rejbad,2)).labels}};
                end

                % plot
                S.fig2=figure('Name',file)
                tiledlayout(1,4)
                nexttile
                rejchanmat = zeros(size(rejbad,1),1);
                rejchanmat(rejchan_all)=1;
                repmat(rejchanmat,1,size(rejbad,2));
                imagesc(rejchanmat,[0 1]); title('noisy chan fully interpolated')
                nexttile
                imagesc(rejbad_all,[0 1]); title('noisy remaining')
                nexttile
                imagesc(rejbad,[0 1]); title('noisy to partially interpolate')
                nexttile
                imagesc(trials_reject,[0 1]); title('noisy trials to be fully removed')
                
                if isfield(S,'fig') savefig(S.fig,fullfile(S.path.prep,S.prep.save.suffix{:},strrep(file,'.set','_noisy.fig'))); end
                if isfield(S,'fig2') savefig(S.fig2,fullfile(S.path.prep,S.prep.save.suffix{:},strrep(file,'.set','_noisy_removed.fig'))); end
                if isfield(S,'fig3') savefig(S.fig3,fullfile(S.path.prep,S.prep.save.suffix{:},strrep(file,'.set','_corrchan.fig'))); end
    
                % REMOVE NOISY CHANNELS using var
                if S.prep.clean.noisychan.clean_data
%                     S.prep.clean.noisychan.parallel_workers = 8;
                    poolobj = gcp('nocreate'); % If no pool, do not create new one.
                    if isempty(poolobj) && S.prep.clean.noisychan.parallel_workers
                        parpool('local',S.prep.clean.noisychan.parallel_workers);
                    elseif ~isempty(poolobj) && S.prep.clean.noisychan.parallel_workers
                        if poolobj.NumWorkers ~= S.prep.clean.noisychan.parallel_workers
                            delete(poolobj)
                            parpool('local',S.prep.clean.noisychan.parallel_workers);
                        end
                    end
                    [EEG,badlist] = reject_channelsbytrial(EEG,rejbad,bs_thresh,S.prep.clean.noisychan.badchans,S.prep.clean.noisychan.plotbads);
                    S.prep.outtable.noisy_chans_Nbadchans(S.fn+f) = badlist.nbadchan;
                    S.prep.outtable.noisy_chans_Nbadtrials(S.fn+f) = badlist.nbadtrial;
                    parfevalOnAll(@clearvars, 0);
                end

                % store bad trials
                EEG.reject.trials = unique([find(all(trials_reject,1)),badlist.badtrial]);

                close all
            end

%             % REMOVE FLAT CHANNELS using 1/var
%             % Only do this after removing very bad noisy channels!
%             % Otherwise interpolation is not reliable
%             if isfield(S.prep.clean,'flatchan') && (any(S.prep.clean.flatchan.trial_weight>0) && any(S.prep.clean.flatchan.chan_weights>0))
%                 S = flat_channel_reject(S,EEG);
%                 rejchan = S.prep.clean.flatchan.rejchan;
%                 rejtrial = S.prep.clean.flatchan.rejtrial;
% 
%                 S.prep.outtable.flat_trials(S.fn+f) = numel(rejtrial);
%                 S.prep.outtable.flat_chans(S.fn+f) = numel(rejchan);
%                 S.prep.outtable.flat_chans_labels(S.fn+f) = {{EEG.chanlocs(rejchan).labels}};
% 
%                 % reject bad channels
%                 EEG = reject_channels(EEG,rejchan);
%                 % reject bad trials
%                 EEG = pop_select(EEG, 'notrial', unique(rejtrial));
%             end

            % semi-auto
            if isfield(S.(S.func).clean,'FTrej') && ~isempty(S.(S.func).clean.FTrej.freq)
                S = FTrejman(EEG,S);%S.(S.func).clean.FTrej.freq{i},S.(S.func).inclchan); 
                rejchan = S.(S.func).clean.FTrej.rejchan;
                rejtrial = S.(S.func).clean.FTrej.rejtrial;

                S.prep.outtable.FTrejman_trials(S.fn+f) = numel(rejtrial);
                S.prep.outtable.FTrejman_chans(S.fn+f) = numel(rejchan);

                % reject bad channels
                EEG = reject_channels(EEG,rejchan);
                % reject bad trials
                EEG = pop_select(EEG, 'notrial', unique(rejtrial));
            end

            % manual check of data and opportunity to reject
            % channels and trials. 
            if isfield(S.prep.clean,'manual_check') && S.prep.clean.manual_check
                pop_eegplot(EEG, 1, 1, 0); % view data first
                ud=get(gcf,'UserData');     %2
                ud.winlength=10;            %3
                set(gcf,'UserData',ud);     %4
                eegplot('draws',0)          %5
                h = gcf;
                waitfor(h) 
            end

            if isfield(S.(S.func),'epoch') && isfield(S.(S.func).epoch,'reref') && S.(S.func).epoch.reref==1
                % re-reference to the common average
                EEG = pop_reref( EEG, []);
            end

            % SAVE
            sname = strrep(file,S.prep.load.suffix{:},S.prep.save.suffix{:});
            disp(['saving file index ' num2str(f)])
            EEG = pop_saveset(EEG,'filename',sname,'filepath',S.prep.save.dir); 
            
            % update
            S.(S.func).file_processed=file;
            S = filehandler(S,'update');
        end

    %% NOISY TRIAL REJECTION
    case 'rejecttrials'

        S = filehandler(S,'start');
           
        for f = S.(S.func).startfile:length(S.prep.filelist)
            file = S.(S.func).filelist{f};
            disp(['loading file index ' num2str(f)])
            try
                EEG = pop_loadset('filename',file,'filepath',S.path.file);
            catch
                continue;
            end

            % collect summary data
            S.prep.outtable.file{S.fn+f} = file;
            S.prep.outtable.initial_trials(S.fn+f) = EEG.trials;

            % remove selected ICA components (MUST HAVE ALREADY SELECTED THESE MANUALLY)
            if S.(S.func).epoch.ICAremove
                if ~strcmp(S.(S.func).load.suffix{:},S.(S.func).load.ica_suffix{:})
                    EEG_ica = pop_loadset('filename',strrep(file,S.(S.func).load.suffix{:},S.(S.func).load.ica_suffix{:}),'filepath',strrep(S.path.file,S.(S.func).load.suffix{:},S.(S.func).load.ica_suffix{:}));
                    EEG.icawinv = EEG_ica.icawinv;
                    EEG.icasphere = EEG_ica.icasphere;
                    EEG.icaweights = EEG_ica.icaweights;
                    EEG.icachansind = EEG_ica.icachansind;
                    EEG.reject = EEG_ica.reject;
                    clear EEG_ica
                end

                S.prep.outtable.ncomp(S.fn+f) = size(EEG.icaweights,1);
                S.prep.outtable.IClabel_nkeep(S.fn+f) = sum(~EEG.reject.gcompreject);
                S.prep.outtable.IClabel_frackeep(S.fn+f) = sum(~EEG.reject.gcompreject)/size(EEG.icaweights,1);

                if any(EEG.reject.gcompreject==0)
                    EEG = pop_subcomp( EEG, [], 0); 
                end
            end

            % detrend
            if S.(S.func).epoch.detrend
                for i = 1:EEG.trials, EEG.data(:,:,i) = detrend(EEG.data(:,:,i)')'; end
            end
    
            % REMOVE BASELINE
            if S.prep.epoch.rmbase
                if iscell(S.prep.epoch.basewin) % markers set by condition/block
                    if ~isempty(S.prep.epoch.basewinby)
                        basewin = S.prep.epoch.basewin{ismember(S.prep.select.(S.prep.epoch.basewinby{:}),S.prep.designtab.(S.prep.epoch.basewinby{:}){f})};
                    else
                        basewin = S.prep.epoch.basewin{:};
                    end
                    EEG = pop_rmbase( EEG, [basewin(1)*1000, basewin(2)*1000]);
                else
                    EEG = pop_rmbase( EEG, [S.prep.epoch.basewin(1)*1000, S.prep.epoch.basewin(2)*1000]);
                end
            end
    
            % select or deselect channels for cleaning
            S.(S.func).select.chans = S.(S.func).clean.chan;
            S=select_chans(S);
            
            % initial manual check of data and opportunity to reject
            % channels and trials. 
            if isfield(S.prep.clean,'initial_manual_check') && S.prep.clean.initial_manual_check
                pop_eegplot(EEG, 1, 1, 0); % view data first
                ud=get(gcf,'UserData');     %2
                ud.winlength=10;            %3
                set(gcf,'UserData',ud);     %4
                eegplot('draws',0)          %5
                h = gcf;
                waitfor(h) 
            end

%             if isfield(S.prep.clean,'man') && S.prep.clean.man
%                 close all
%                 % initial arrays
%                 rejchan = [];
%                 rejtrial = [];
%                 
%                 S.(S.func).clean.man.hidebad=0;
%                 [EEG,S] = manualRejectEEG(EEG,S,S.(S.func).clean.man.thresh); 
%                 rejchan = [rejchan;S.(S.func).clean.man.rejchan];
%                 rejtrial = [rejtrial;S.(S.func).clean.man.rejtrial];
%                 
%                 S.(S.func).clean.man.hidebad=1;
%                 [EEG,S] = manualRejectEEG(EEG,S,S.(S.func).clean.man.thresh); 
%                 rejchan = [rejchan;S.(S.func).clean.man.rejchan];
%                 rejtrial = [rejtrial;S.(S.func).clean.man.rejtrial];
% 
%                 S.prep.outtable.man_trials(S.fn+f) = numel(rejtrial);
%                 S.prep.outtable.man_chans(S.fn+f) = numel(rejchan);
%                 
%                 close all
%                 % reject bad channels
%                 EEG = reject_channels(EEG,rejchan);
%                 % reject bad trials
%                 EEG = pop_select(EEG, 'notrial', unique(rejtrial));
%             end

            % reject trials from previous chancorrect step
            if S.(S.func).clean.reject_trials_from_chancorrect && isfield(EEG.reject,'trials')
                S.(S.func).outtable.Nbadtrials_from_chancorrect(S.fn+f) = length(EEG.reject.trials);
                S.(S.func).outtable.Fracbadtrials_from_chancorrect(S.fn+f) = length(EEG.reject.trials)/S.prep.outtable.initial_trials(S.fn+f);
                EEG = pop_select(EEG, 'notrial', EEG.reject.trials);
            end

            % Auto trial rejection
            % first run without ignoring artefact, then with, and only
            % exclude trials rejected by both in case excluding data causes
            % artificial jump artfacts, e.g. due to EOG.
            if isfield(S.(S.func).clean,'FTrejauto')
                ignore_stim_artefact = S.(S.func).clean.FTrejauto.ignore_stim_artefact;
                if any(ignore_stim_artefact) 
                    S.(S.func).clean.FTrejauto.ignore_stim_artefact = 0;
                    S = FTrejauto(EEG,S);%...
%                     S.(S.func).clean.FTrejauto.cutoffs,...
%                     S.(S.func).clean.FTrejauto.interactive,...
%                     S.(S.func).inclchan); 
                    rejtrial_with_artefact = S.(S.func).clean.FTrejauto.rejtrial;
                    S.prep.outtable.FTrejauto_trials_withstimartefact_jump(S.fn+f) = S.(S.func).clean.FTrejauto.Nrejtrial_by_type.jump;
                    S.prep.outtable.FTrejauto_trials_withstimartefact_muscle(S.fn+f) = S.(S.func).clean.FTrejauto.Nrejtrial_by_type.muscle;
                    S.prep.outtable.FTrejauto_trials_withstimartefact_eog(S.fn+f) = S.(S.func).clean.FTrejauto.Nrejtrial_by_type.eog;
                end
                S.(S.func).clean.FTrejauto.ignore_stim_artefact=ignore_stim_artefact;
                S = FTrejauto(EEG,S);%...
%                     S.(S.func).clean.FTrejauto.cutoffs,...
%                     S.(S.func).clean.FTrejauto.interactive,...
%                     S.(S.func).inclchan); 
                rejtrial = S.(S.func).clean.FTrejauto.rejtrial;

                if any(ignore_stim_artefact) 
                    rejtrial = intersect(rejtrial_with_artefact,rejtrial);
                end

                S.prep.outtable.FTrejauto_trials(S.fn+f) = numel(rejtrial);
                S.prep.outtable.FTrejauto_trials_frac(S.fn+f) = numel(rejtrial)/S.prep.outtable.initial_trials(S.fn+f);
                S.prep.outtable.FTrejauto_trials_jump(S.fn+f) = S.(S.func).clean.FTrejauto.Nrejtrial_by_type.jump;
                S.prep.outtable.FTrejauto_trials_muscle(S.fn+f) = S.(S.func).clean.FTrejauto.Nrejtrial_by_type.muscle;
                S.prep.outtable.FTrejauto_trials_eog(S.fn+f) = S.(S.func).clean.FTrejauto.Nrejtrial_by_type.eog;

                % reject bad trials
                try
                    EEG = pop_select(EEG, 'notrial', unique(rejtrial));
                catch
                    continue % too many trials rejected
                end
            end

            % semi-auto
            if isfield(S.(S.func).clean,'FTrej') && ~isempty(S.(S.func).clean.FTrej.freq)
                S = FTrejman(EEG,S);%S.(S.func).clean.FTrej.freq{i},S.(S.func).inclchan); 
                rejchan = S.(S.func).clean.FTrej.rejchan;
                rejtrial = S.(S.func).clean.FTrej.rejtrial;

                S.prep.outtable.FTrejman_trials(S.fn+f) = numel(rejtrial);
                S.prep.outtable.FTrejman_chans(S.fn+f) = numel(rejchan);

                % reject bad channels
                EEG = reject_channels(EEG,rejchan);
                % reject bad trials
                EEG = pop_select(EEG, 'notrial', unique(rejtrial));
            end

            % manual check of data and opportunity to reject
            % channels and trials. 
            if isfield(S.prep.clean,'manual_check') && S.prep.clean.manual_check
                pop_eegplot(EEG, 1, 1, 0); % view data first
                ud=get(gcf,'UserData');     %2
                ud.winlength=10;            %3
                set(gcf,'UserData',ud);     %4
                eegplot('draws',0)          %5
                h = gcf;
                waitfor(h) 
            end

            if isfield(S.(S.func),'epoch') && isfield(S.(S.func).epoch,'reref') && S.(S.func).epoch.reref==1
                % re-reference to the common average
                EEG = pop_reref( EEG, []);
            end

            % SAVE
            sname = strrep(file,S.prep.load.suffix{:},S.prep.save.suffix{:});
            disp(['saving file index ' num2str(f)])
            EEG = pop_saveset(EEG,'filename',sname,'filepath',S.prep.save.dir); 
            
            % update
            S.(S.func).file_processed=file;
            S = filehandler(S,'update');
        end

%% ICA
    case 'ICA'
        if S.prep.clean.ICA

%             % previously processed files
%             if isfield(S.(S.func),'filelist')
%                 if S.(S.func).overwrite==0
%                     prev_filelist = S.(S.func).filelist;
%                     prev_dirlist = S.(S.func).dirlist;
%                     prev_subj_pdat_idx = S.(S.func).subj_pdat_idx;
%                     prev_designtab = S.(S.func).designtab;
%                 end
%                 S.(S.func) = rmfield(S.(S.func),'dirlist');
%                 S.(S.func) = rmfield(S.(S.func),'subj_pdat_idx');
%                 S.(S.func) = rmfield(S.(S.func),'designtab');
%             end
% 
%             % GET FILE LIST 
%             S.path.file = fullfile(S.path.prep,S.prep.load.suffix{:});
%             S = getfilelist(S,S.prep.load.suffix);
% 
%             % select which to process
%             if S.(S.func).overwrite==0 && exist('prev_filelist','var')
%                 idx_remove = ismember(S.(S.func).filelist,prev_filelist);
%                 S.(S.func).filelist(idx_remove)=[];
%                 S.(S.func).dirlist(idx_remove)=[];
%                 S.(S.func).subj_pdat_idx(idx_remove)=[];
%                 S.(S.func).designtab(idx_remove,:)=[];
%             end
    
            S = filehandler(S,'start');

            for f = S.prep.startfile:length(S.prep.filelist)
                file = S.prep.filelist{f};
                EEG = pop_loadset('filename',file,'filepath',S.path.file);

                % detrend
                if isfield(S.(S.func),'epoch') && isfield(S.(S.func).epoch,'detrend') && S.(S.func).epoch.detrend
                    for i = 1:EEG.trials, EEG.data(:,:,i) = detrend(EEG.data(:,:,i)')'; end
                end
        
                % REMOVE BASELINE
                if S.prep.epoch.rmbase
                    if iscell(S.prep.epoch.basewin) % markers set by condition/block
                        if ~isempty(S.prep.epoch.basewinby)
                            basewin = S.prep.epoch.basewin{ismember(S.prep.select.(S.prep.epoch.basewinby{:}),S.prep.designtab.(S.prep.epoch.basewinby{:}){f})};
                        else
                            basewin = S.prep.epoch.basewin{:};
                        end
                        EEG = pop_rmbase( EEG, [basewin(1)*1000, basewin(2)*1000]);
                    else
                        EEG = pop_rmbase( EEG, [S.prep.epoch.basewin(1)*1000, S.prep.epoch.basewin(2)*1000]);
                    end
                end           
    
                %RUN ICA 
                if isfield(S.prep,'clean') && isfield(S.prep.clean,'ignore_stim_artefact') && ~isempty(S.prep.clean.ignore_stim_artefact)
                    art_time = S.prep.clean.ignore_stim_artefact;
                    art_dp = dsearchn(EEG.times', [S.prep.clean.ignore_stim_artefact*1000]');
                    EEGcut = EEG;
                    EEGcut.times(art_dp(1):art_dp(2))=[];
                    EEGcut.data(:,art_dp(1):art_dp(2),:) = [];
                    EEGcut.pnts = length(EEGcut.times);
                
                    numcomp = numcompeig(EEGcut);
                    EEGcut = pop_runica(EEGcut, 'extended',1,'interupt','off','pca',numcomp);
                
                    EEG.icaact = EEGcut.icaact;
                    EEG.icawinv = EEGcut.icawinv;
                    EEG.icasphere = EEGcut.icasphere;
                    EEG.icaweights = EEGcut.icaweights;
                    EEG.icachansind = EEGcut.icachansind;
                    EEG.icasplinefile = EEGcut.icasplinefile;
                else
                    numcomp = numcompeig(EEG);
                    EEG = pop_runica(EEG, 'extended',1,'interupt','off','pca',numcomp);
                end
    
                % SAVE
                %[pth nme ext] = fileparts(file);
                %sname = [nme '_' S.prep.save.suffix{:} '.' S.prep.fname.ext{:}];
                sname = strrep(file,S.prep.load.suffix{:},S.prep.save.suffix{:});
%                 if ~exist(fullfile(S.path.prep,S.prep.save.suffix{:}),'dir')
%                     mkdir(fullfile(S.path.prep,S.prep.save.suffix{:}));
%                 end
                EEG = pop_saveset(EEG,'filename',sname,'filepath',S.prep.save.dir);

                % update
                S.(S.func).file_processed=file;
                S = filehandler(S,'update');

            end

        end

    case 'SASICA'

%         % previously processed files
%         if isfield(S.(S.func),'filelist')
%             if S.(S.func).overwrite==0
%                 prev_filelist = S.(S.func).filelist;
%                 prev_dirlist = S.(S.func).dirlist;
%                 prev_subj_pdat_idx = S.(S.func).subj_pdat_idx;
%                 prev_designtab = S.(S.func).designtab;
%             end
%             S.(S.func) = rmfield(S.(S.func),'dirlist');
%             S.(S.func) = rmfield(S.(S.func),'subj_pdat_idx');
%             S.(S.func) = rmfield(S.(S.func),'designtab');
%         end
%         
%         % GET FILE LIST 
%         S.path.file = fullfile(S.path.prep,S.prep.load.suffix{:});
%         S = getfilelist(S,S.prep.load.suffix);
%         
%         % select which to process
%         if S.(S.func).overwrite==0 && exist('prev_filelist','var')
%             idx_remove = ismember(S.(S.func).filelist,prev_filelist);
%             S.(S.func).filelist(idx_remove)=[];
%             S.(S.func).dirlist(idx_remove)=[];
%             S.(S.func).subj_pdat_idx(idx_remove)=[];
%             S.(S.func).designtab(idx_remove,:)=[];
%         end
%         
%         % setup summary table
%         if ~isfield(S.(S.func).clean,'ICA') || S.(S.func).overwrite==1
%             S.(S.func).clean.ICA = table;
%         end
%         fn = height(S.(S.func).clean.ICA);

        S = filehandler(S,'start');

        for f = S.prep.startfile:length(S.prep.filelist)
            file = S.(S.func).filelist{f};
            EEG = pop_loadset('filename',file,'filepath',S.path.file);
           
%             %RUN ICA on auto
%             [ALLEEG EEG CURRENTSET ALLCOM] = eeglab('nogui'); 
%             EEG = pop_loadset('filename',file,'filepath',S.path.file);
%             [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, 0); 
%             SASICA(EEG)
%             load SASICA_obj
            def=SASICAsettings_MNP(~S.(S.func).clean.sisica_plot); % definitions at end of this function
            EEG = eeg_SASICA(EEG,def);
%             EEG=SASICA(EEG,'opts_nocompute',1);
% 
            if S.(S.func).clean.sisica_plot
                pause
            end

            % save info to a table
            S.prep.outtable.file{S.fn+f} = file;
            S.prep.outtable.ncomp(S.fn+f) = size(EEG.icaweights,1);
            S.prep.outtable.sasica_nreject(S.fn+f) = sum(EEG.reject.gcompreject);
            S.prep.outtable.sasica_fracreject(S.fn+f) = sum(EEG.reject.gcompreject)/size(EEG.icaweights,1);

            % SAVE
            %[pth nme ext] = fileparts(file);
            %sname = [nme '_' S.prep.save.suffix{:} '.' S.prep.fname.ext{:}];
            sname = strrep(file,S.prep.load.suffix{:},S.prep.save.suffix{:});
%             if ~exist(fullfile(S.path.prep,S.prep.save.suffix{:}),'dir')
%                 mkdir(fullfile(S.path.prep,S.prep.save.suffix{:}));
%             end
            EEG = pop_saveset(EEG,'filename',sname,'filepath',S.prep.save.dir);
            
            % update
            S.(S.func).file_processed=file;
            S = filehandler(S,'update');
        end
    
%         if exist('prev_filelist','var')
%             S.(S.func).filelist = [S.(S.func).filelist prev_filelist];
%             S.(S.func).dirlist = [S.(S.func).dirlist prev_dirlist];
%             S.(S.func).subj_pdat_idx = [S.(S.func).subj_pdat_idx prev_subj_pdat_idx];
%             S.(S.func).designtab = [S.(S.func).designtab; prev_designtab];
%         end

    case 'IClabel'

%         % previously processed files
%         if isfield(S.(S.func),'filelist')
%             if S.(S.func).overwrite==0
%                 prev_filelist = S.(S.func).filelist;
%                 prev_dirlist = S.(S.func).dirlist;
%                 prev_subj_pdat_idx = S.(S.func).subj_pdat_idx;
%                 prev_designtab = S.(S.func).designtab;
%             end
%             S.(S.func) = rmfield(S.(S.func),'dirlist');
%             S.(S.func) = rmfield(S.(S.func),'subj_pdat_idx');
%             S.(S.func) = rmfield(S.(S.func),'designtab');
%         end
% 
%         % GET FILE LIST 
%         S.path.file = fullfile(S.path.prep,S.prep.load.suffix{:});
%         S = getfilelist(S,S.prep.load.suffix);
% 
%         % select which to process
%         if S.(S.func).overwrite==0 && exist('prev_filelist','var')
%             idx_remove = ismember(S.(S.func).filelist,prev_filelist);
%             S.(S.func).filelist(idx_remove)=[];
%             S.(S.func).dirlist(idx_remove)=[];
%             S.(S.func).subj_pdat_idx(idx_remove)=[];
%             S.(S.func).designtab(idx_remove,:)=[];
%         end
% 
%         % setup summary table
%         if ~isfield(S.(S.func).clean,'ICA') || S.(S.func).overwrite==1
%             S.(S.func).clean.ICA = table;
%         end
%         fn = height(S.(S.func).clean.ICA);
% 
%         loadpath = S.path.file;

        S = filehandler(S,'start');

        for f = S.prep.startfile:length(S.prep.filelist)
            file = S.(S.func).filelist{f};
            EEG = pop_loadset('filename',file,'filepath',S.path.file);
           
            EEG = iclabel(EEG);
            getbrain = strcmp(EEG.etc.ic_classification.ICLabel.classes,'Brain');
            brainvals = EEG.etc.ic_classification.ICLabel.classifications(:,getbrain);
            brainmax = brainvals==max(EEG.etc.ic_classification.ICLabel.classifications,[],2);
            if S.prep.clean.brain_threshold>0 
                EEG.reject.gcompreject = [brainvals<S.prep.clean.brain_threshold]';
            else
                EEG.reject.gcompreject = ~brainmax';
            end

            if S.prep.clean.IClabel_plot
                pop_viewprops(EEG, 0) % for component properties
            end

            % save info to a table
            S.prep.outtable.file{S.fn+f} = file;
            S.prep.outtable.ncomp(S.fn+f) = size(EEG.icaweights,1);
            S.prep.outtable.IClabel_nOverThresh(S.fn+f) = sum([brainvals>=S.prep.clean.brain_threshold]);
            S.prep.outtable.IClabel_nBrainMax(S.fn+f) = sum(brainmax);
            S.prep.outtable.IClabel_nreject(S.fn+f) = sum(EEG.reject.gcompreject);
            S.prep.outtable.IClabel_nkeep(S.fn+f) = sum(~EEG.reject.gcompreject);
            S.prep.outtable.IClabel_fracreject(S.fn+f) = sum(EEG.reject.gcompreject)/size(EEG.icaweights,1);
            
            % save all classification values
            S.prep.outtable.IClabel_classes{S.fn+f} = EEG.etc.ic_classification.ICLabel.classes;
            S.prep.outtable.IClabel_classifications{S.fn+f} = EEG.etc.ic_classification.ICLabel.classifications;
            
            
            % SAVE
            %[pth nme ext] = fileparts(file);
            %sname = [nme '_' S.prep.save.suffix{:} '.' S.prep.fname.ext{:}];
            sname = strrep(file,S.prep.load.suffix{:},S.prep.save.suffix{:});
%             if ~exist(fullfile(S.path.prep,S.prep.save.suffix{:}),'dir')
%                 mkdir(fullfile(S.path.prep,S.prep.save.suffix{:}));
%             end
            EEG = pop_saveset(EEG,'filename',sname,'filepath',S.prep.save.dir);

            % update
            S.(S.func).file_processed=file;
            S = filehandler(S,'update');
        end
    
   
    %% separate files
    case 'sep'
%         % previously processed files
%         if isfield(S.(S.func),'filelist')
%             if S.(S.func).overwrite==0
%                 prev_filelist = S.(S.func).filelist;
%                 prev_dirlist = S.(S.func).dirlist;
%                 prev_subj_pdat_idx = S.(S.func).subj_pdat_idx;
%                 prev_designtab = S.(S.func).designtab;
%             end
%             S.(S.func) = rmfield(S.(S.func),'dirlist');
%             S.(S.func) = rmfield(S.(S.func),'subj_pdat_idx');
%             S.(S.func) = rmfield(S.(S.func),'designtab');
%         end
% 
%         % GET FILE LIST
%         S.path.file = fullfile(S.path.prep,S.prep.load.suffix{:});
%         S = getfilelist(S,S.prep.load.suffix);
% 
%         % select which to process
%         if S.(S.func).overwrite==0 && exist('prev_filelist','var')
%             idx_remove = ismember(S.(S.func).filelist,prev_filelist);
%             S.(S.func).filelist(idx_remove)=[];
%             S.(S.func).dirlist(idx_remove)=[];
%             S.(S.func).subj_pdat_idx(idx_remove)=[];
%             S.(S.func).designtab(idx_remove,:)=[];
%         end
% 
%         sname_ext = S.(S.func).save.suffix{:};
%         if ~exist(fullfile(S.path.prep,sname_ext),'dir')
%             mkdir(fullfile(S.path.prep,sname_ext));
%         end
% 
%         loadpath = S.path.file;

        S = filehandler(S,'start');
        sname_ext = S.(S.func).save.suffix{:};

        for f = 1:length(S.(S.func).filelist)
            file = S.(S.func).filelist{f};
            disp(['loading file index ' num2str(f)])

            INEEG = pop_loadset('filename',file,'filepath',S.path.file);
            [sfiles,ia,ib] = unique({INEEG.epoch.file});
            for s = 1:length(sfiles)
                ind = find(ib==s);
                EEGsep1 = pop_select(INEEG,'trial',ind);
                
                % SEPARATE INTO FILES ACCORDING TO MARKER TYPE AND SAVE
                if ~isempty(S.(S.func).separatefiles.markerindex)
                    nfiles = length(S.(S.func).separatefiles.markerindex);
                    INEEGsep1 =EEGsep1; 
                    allmarkers = {INEEGsep1.epoch.eventtype};
                    for n = 1:nfiles
                        selectmarkers = S.(S.func).epoch.markers(S.(S.func).separatefiles.markindex{n});
                        markerindex = find(ismember(allmarkers,selectmarkers));
                        EEG = pop_select(INEEGsep1,'trial',markerindex);
    
                        % save .set LIKELY NEEDS UPDATING TO BE SIMILAR TO
                        % AFTER THE 'ELSE'
                        sname = [sfiles{s} '_' S.(S.func).separatefiles.suffix{n} '.' S.(S.func).fname.ext{:}];
                        pop_saveset(EEG,'filename',sname,'filepath',fullfile(S.path.prep,sname_ext)); 
                    end
                else
                    % save as one file
                    sname = strrep(sfiles{s},'epoched','cleaned');

                    % put into sub-directory according to original
                    % block/condition (needed for MoNoPly study)
                    spath = fullfile(S.prep.save.dir,[S.prep.designtab.blocks{f} '_' S.prep.designtab.conds{f}]);
                    if ~exist(spath,'dir')
                        mkdir(spath);
                    end

                    EEG=EEGsep1;
                    pop_saveset(EEG,'filename',sname,'filepath',spath);
                end 
            end
            
            % update
            S.(S.func).file_processed=file;
            S = filehandler(S,'update');
             
        end

%         if exist('prev_filelist','var')
%             S.(S.func).filelist = [S.(S.func).filelist prev_filelist];
%             S.(S.func).dirlist = [S.(S.func).dirlist prev_dirlist];
%             S.(S.func).subj_pdat_idx = [S.(S.func).subj_pdat_idx prev_subj_pdat_idx];
%             S.(S.func).designtab = [S.(S.func).designtab; prev_designtab];
%         end
end

%% REJECT CHANNELS
function EEG = reject_channels(EEG,rejchan)
% reject bad channels
EEG = eeg_interp(EEG, unique(rejchan));
if length(EEG.chanlocs)>EEG.nbchan
    EEG.chanlocs(unique(rejchan))=[];
end
% make a record in EEG.chanlocs
if ~isfield(EEG.chanlocs,'interp')
    interp = zeros(EEG.nbchan,1);
else
    interp = cell2mat({EEG.chanlocs.interp});
end
interp(unique(rejchan))=1;
interp=num2cell(interp);
[EEG.chanlocs.interp] = interp{:};

%% REJECT CHANNELS BY TRIAL
function [EEG,badlist] = reject_channelsbytrial(EEG,rejbad,badsegs,badchans,plot_bads)
% rejbad is a boolean matrix of chans*trials
bads=tbt_bool2cell(rejbad,EEG);
[EEG, com, badlist] = pop_TBT_cab(EEG,bads,badchans,badsegs,plot_bads);

%% SASICA SETTINGS
function def=SASICAsettings_MNP(plot)

def.autocorr.enable = true;
def.autocorr.dropautocorr = 0.1;
def.autocorr.autocorrint = 20;% will compute autocorrelation with this many milliseconds lag

def.focalcomp.enable = true;
def.focalcomp.focalICAout = 7;

def.trialfoc.enable = true;
def.trialfoc.focaltrialout = 20;

def.resvar.enable = false;
def.resvar.thresh = 15;% %residual variance allowed

def.SNR.enable = true;
def.SNR.snrPOI = [20 Inf];% period of interest (signal) in ms
def.SNR.snrBL = [-500 -80];% period of no interest (noise) in ms
def.SNR.snrcut = 1;% SNR below this threshold will be dropped

def.EOGcorr.enable = false;
def.EOGcorr.corthreshV = 'auto 4';% threshold correlation with vertical EOG
def.EOGcorr.Veogchannames = [];% vertical channel(s)
def.EOGcorr.corthreshH = 'auto 4';% threshold correlation with horizontal EOG
def.EOGcorr.Heogchannames = [];% horizontal channel(s)

def.chancorr.enable = false;
def.chancorr.corthresh = 'auto 4';% threshold correlation
def.chancorr.channames = [];% channel(s)

def.FASTER.enable = true;
def.FASTER.blinkchanname = 'Fp1';

def.ADJUST.enable = true;

def.MARA.enable = false;

def.opts.FontSize = 14;
def.opts.noplot = plot;
def.opts.nocompute = 0;