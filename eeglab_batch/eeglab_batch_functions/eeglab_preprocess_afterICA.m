function S=eeglab_preprocess_afterICA(S)
% PREPROCESS AFTER ICA
% To be run after ICA components have been manually selected for rejection
% via EEGLAB or SASICA, but not yet removed from the data.

S.func = 'prep2';

% GET FILE LIST
S.path.file = S.path.prep;
sname_ext = 'ICA';
S.path.file = fullfile(S.path.prep,sname_ext);
S = getfilelist(S);

loadpath = fullfile(S.path.prep,sname_ext);
for f = S.(S.func).startfile:length(S.(S.func).filelist)
    file = S.(S.func).filelist{f};
    EEG = pop_loadset('filename',file,'filepath',loadpath);
    
    % remove selected ICA components (MUST HAVE ALREADY SELECTED THESE MANUALLY)
    if S.(S.func).epoch.ICAremove && any(EEG.reject.gcompreject==0)
        EEG = pop_subcomp( EEG, [], 0); 
    end
    
    % detrend
    if S.(S.func).epoch.detrend
        for i = 1:EEG.trials, EEG.data(:,:,i) = detrend(EEG.data(:,:,i)')'; end
    end
    
    % remove baseline
    if S.(S.func).epoch.rmbase
        EEG = pop_rmbase( EEG, [S.(S.func).epoch.basewin(1)*1000 S.(S.func).epoch.basewin(2)*1000]);
    end
    
    % select or deselect channels for cleaning
    S.(S.func).select.chans = S.(S.func).clean.chan;
    S=select_chans(S);
    
    % initial arrays
    rejchan = [];
    rejtrial = [];

    % strategy - only remove a very small number of very bad trials / chans
    % before ICA - do further cleaning after ICA
    if isfield(S.(S.func).clean,'FTrej')
        S = FTrejman(EEG,S);%S.(S.func).clean.FTrej.freq{i},S.(S.func).inclchan); 
        rejchan = [rejchan;S.(S.func).clean.FTrej.rejchan];
        rejtrial = [rejtrial;S.(S.func).clean.FTrej.rejtrial];
    end

    % reject bad channels
    EEG = eeg_interp(EEG, unique(rejchan));
    if length(EEG.chanlocs)>EEG.nbchan
        EEG.chanlocs(unique(rejchan))=[];
    end

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
    
    % additional manual check of data and opportunity to reject
    % further channels and trials
    if isfield(S.(S.func).clean,'manual_check') && S.(S.func).clean.manual_check
        eeglab
        [EEG,S] = manualRejectEEG(EEG,S,0);
    end
    
    % re-reference to the common average
    if S.(S.func).epoch.reref == 1
        % re-reference to the common average
        EEG = pop_reref( EEG, []); 
    elseif S.(S.func).reref==2
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
   
    % save .set
    [pth nme ext] = fileparts(file); 
    sname_ext = 'cleaned';
    sname = [nme '_' sname_ext '.' S.(S.func).fname.ext{:}];
    if ~exist(fullfile(S.path.prep,sname_ext),'dir')
        mkdir(fullfile(S.path.prep,sname_ext));
    end
    EEG = pop_saveset(EEG,'filename',sname,'filepath',fullfile(S.path.prep,sname_ext)); 
    
end

S=eeglab_separate_files(S)
