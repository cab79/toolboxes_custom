function S=eeglab_preprocess(S,part)
%% PREPROCESSING FOR CONTINUOUS EEGLAB .SET FILES
S.func = 'prep';

switch part
    case 'epoch'
        
    % load directory
    if ~isfield(S.prep,'loaddir')
        S.prep.loaddir = fullfile(S.path.(S.func),S.prep.load.suffix{:});
    end

    % GET FILE LIST
    S.path.file = S.prep.loaddir;
    if ~isempty(S.prep.load.suffix{:})
        S = getfilelist(S,S.prep.load.suffix);
    else
        S = getfilelist(S);
    end

    % change to the input directory
    eval(sprintf('%s', ['cd(''' S.path.prep ''')']));

    % report if there are no such files
    if isempty(S.prep.filelist)
        error('No files found!\n');
    end

    % run though all files in a loop
    for f = S.prep.startfile:length(S.prep.filelist)
        filename = S.prep.filelist{f};

        fprintf('\nProcessing %s.\n\n', filename);

        % GET FILENAME PARTS
        [pth nme ext] = fileparts(filename); 
        fparts = strsplit(nme,'_');

        % LOAD DATA
        EEG = pop_loadset('filename',filename,'filepath',S.path.file);

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
            params.lineFrequencies = [];
            params.detrendType = 'high pass';
            params.detrendCutoff = 1;
            params.referenceType = 'robust';
            params.meanEstimateType = 'median';
            params.interpolationOrder = 'post-reference';
            params.keepFiltered = false;
            params.ignoreBoundaryEvents = true;
            [EEG, computationTimes] = prepPipeline(EEG, params); 
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
        if ~isempty(S.prep.epoch.markers)
            try
                EEG = pop_epoch( EEG, S.prep.epoch.markers, S.prep.epoch.timewin);
            end
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


        % LINEAR DETREND
        if S.prep.epoch.detrend
            for i = 1:EEG.trials, EEG.data(:,:,i) = detrend(EEG.data(:,:,i)')'; end
        end

        % REMOVE BASELINE
        if S.prep.epoch.rmbase
            EEG = pop_rmbase( EEG, [S.prep.epoch.basewin(1)*1000, S.prep.epoch.basewin(2)*1000]);
        end

        EEG = eeg_checkset( EEG );
        
        % additional manual check of data and opportunity to reject
        % further channels and trials
        if isfield(S.prep.epoch,'manual_check') && S.prep.epoch.manual_check
            [EEG,S] = manualRejectEEG(EEG,S,0);
        end

        % SAVE
        sname = [nme '_' S.prep.save.suffix{:} '.' S.prep.fname.ext{:}];
        if ~exist(fullfile(S.path.prep,S.prep.save.suffix{:}),'dir')
            mkdir(fullfile(S.path.prep,S.prep.save.suffix{:}));
        end
        EEG = pop_saveset(EEG,'filename',sname,'filepath',fullfile(S.path.prep,S.prep.save.suffix{:})); 
    end

%% COMBINE
    case 'combine'
    

    % combine data files
    if S.prep.combinefiles.on
        clear OUTEEG

        % update last filenameparts
        for fnp = 1:length(S.prep.fname.parts)
            try
                for i = 1:length(S.prep.select.([S.prep.fname.parts{fnp} 's']))

                    S.prep.select.([S.prep.fname.parts{fnp} 's']){i} = strrep(S.prep.select.([S.prep.fname.parts{fnp} 's']){i},'.','_');

                    %S.prep.([S.prep.fname.parts{end} 's']){i} = [S.prep.([S.prep.fname.parts{end} 's']){i} '_' sname_ext];
                end
            end
        end

        % GET FILE LIST
        S.path.file = fullfile(S.path.prep,S.prep.load.suffix{:});
        S = getfilelist(S,S.prep.load.suffix);

        loadpath = fullfile(S.path.prep,S.prep.load.suffix{:});
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
                    EEG = pop_loadset('filename',file,'filepath',loadpath);
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
                if ~exist(fullfile(S.path.prep,S.prep.save.suffix{:}),'dir')
                    mkdir(fullfile(S.path.prep,S.prep.save.suffix{:}));
                end
                EEG = pop_saveset(EEG,'filename',sname,'filepath',fullfile(S.path.prep,S.prep.save.suffix{:}));
            end
        end
    elseif exist(fullfile(S.path.prep,'combined'),'dir')
        S.prep.save.suffix{:} = 'combined';
    end

% FIND THE NEW DATA FILES
%files = dir(fullfile(S.path.prep,['*' S.runsubject '*' S.prep.load.suffix{:}]));
%files_ana = 1:length(files);

%% NOISY TRIAL AND CHANNEL REJECTION USING FIELDTRIP
    case 'rej'
    %if ~isempty(S.prep.clean.FTrej) && iscell(S.prep.clean.FTrej)

        % GET FILE LIST
        S.path.file = fullfile(S.path.prep,S.prep.load.suffix{:});
        if strcmp(S.prep.load.suffix{:},'combined')
            S.prep.fname.parts = {'study','subject','session','suffix','ext'};
            S.prep.select.conds = {};
            S.prep.select.blocks = {};
        end
        S = getfilelist(S,S.prep.load.suffix);

        loadpath = S.path.file;
        for f = S.prep.startfile:length(S.prep.filelist)
            file = S.prep.filelist{f};
            disp(['loading file index ' num2str(f)])
            EEG = pop_loadset('filename',file,'filepath',loadpath);
            
            % select or deselect channels for cleaning
            S.(S.func).select.chans = S.(S.func).clean.chan;
            S=select_chans(S);
            
            % initial arrays
            rejchan = [];
            rejtrial = [];
            
            if isfield(S.prep.clean,'man')
                close all
                
                S.(S.func).clean.man.hidebad=0;
                [EEG,S] = manualRejectEEG(EEG,S,S.(S.func).clean.man.thresh); 
                rejchan = [rejchan;S.(S.func).clean.man.rejchan];
                rejtrial = [rejtrial;S.(S.func).clean.man.rejtrial];
                
                S.(S.func).clean.man.hidebad=1;
                [EEG,S] = manualRejectEEG(EEG,S,S.(S.func).clean.man.thresh); 
                rejchan = [rejchan;S.(S.func).clean.man.rejchan];
                rejtrial = [rejtrial;S.(S.func).clean.man.rejtrial];
                
                close all
            end
                
            % strategy - only remove a very small number of very bad trials / chans
            % before ICA - do further cleaning after ICA
            if isfield(S.(S.func).clean,'FTrej')
                S = FTrejman(EEG,S);%S.(S.func).clean.FTrej.freq{i},S.(S.func).inclchan); 
                rejchan = [rejchan;S.(S.func).clean.FTrej.rejchan];
                rejtrial = [rejtrial;S.(S.func).clean.FTrej.rejtrial];
            end
            
            % REMOVE FLAT CHANNELS using 1/var
            % Only do this after removing very bad noisy channels!
            % Otherwise interpolation is not reliable
            if isfield(S.prep.clean,'flatchan') && (any(S.prep.clean.flatchan.trial_weight>0) && any(S.prep.clean.flatchan.chan_weights>0))
                S = flat_channel_reject(S,EEG);
                rejchan = [rejchan;S.prep.clean.flatchan.rejchan];
                rejtrial = [rejtrial;S.prep.clean.flatchan.rejtrial];
            end
            
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
            
            % Auto trial rejection
            if isfield(S.(S.func).clean,'FTrejauto')
                S = FTrejauto(EEG,S);%...
%                     S.(S.func).clean.FTrejauto.cutoffs,...
%                     S.(S.func).clean.FTrejauto.interactive,...
%                     S.(S.func).inclchan); 
                rejtrial = [rejtrial;S.(S.func).clean.FTrejauto.rejtrial];
            end
            
            % reject bad trials
            EEG = pop_select(EEG, 'notrial', unique(rejtrial));
            
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

            % SAVE
            %[pth nme ext] = fileparts(file); 
            %sname = [nme '_' S.prep.save.suffix{:} '.' S.prep.fname.ext{:}];
            sname = strrep(file,S.prep.load.suffix{:},S.prep.save.suffix{:});
            if ~exist(fullfile(S.path.prep,S.prep.save.suffix{:}),'dir')
                mkdir(fullfile(S.path.prep,S.prep.save.suffix{:}));
            end
            disp(['saving file index ' num2str(f)])
            EEG = pop_saveset(EEG,'filename',sname,'filepath',fullfile(S.path.prep,S.prep.save.suffix{:})); 
        end
    %elseif exist(fullfile(S.path.prep,S.prep.save.suffix{:}),'dir')
    %    S.prep.save.suffix{:} = 'manrej';
    %end
    
%% separate prior to ICA
    case 'sep'
        S=eeglab_separate_files(S)

%% ICA
    case 'ICA'
    if S.prep.clean.ICA
        % GET FILE LIST 
        S.path.file = fullfile(S.path.prep,S.prep.load.suffix{:});
        S = getfilelist(S,S.prep.load.suffix);

        loadpath = S.path.file;
        for f = S.prep.startfile:length(S.prep.filelist)
            file = S.prep.filelist{f};
            EEG = pop_loadset('filename',file,'filepath',loadpath);
           

            %RUN ICA 
            numcomp = numcompeig(EEG);
            EEG = pop_runica(EEG, 'extended',1,'interupt','on','pca',numcomp);

            % SAVE
            %[pth nme ext] = fileparts(file);
            %sname = [nme '_' S.prep.save.suffix{:} '.' S.prep.fname.ext{:}];
            sname = strrep(file,S.prep.load.suffix{:},S.prep.save.suffix{:});
            if ~exist(fullfile(S.path.prep,S.prep.save.suffix{:}),'dir')
                mkdir(fullfile(S.path.prep,S.prep.save.suffix{:}));
            end
            EEG = pop_saveset(EEG,'filename',sname,'filepath',fullfile(S.path.prep,S.prep.save.suffix{:}));
        end
    end
end