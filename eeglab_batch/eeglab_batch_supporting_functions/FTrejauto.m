function S = FTrejauto(EEG,S,varargin)
%convert to FT and manually remove chan and trial condidering a particular freq

chan = S.(S.func).inclchan;

if ~isempty(varargin)
    format = varargin{1};
else
    format = 'EEGLAB';
end

if S.(S.func).clean.FTrejauto.interactive
    interactive='yes';
else
    interactive = 'no';
end

if strcmp(format,'EEGLAB')
    FT = convertoft(EEG);
elseif strcmp(format,'SPM')
    FT = EEG.fttimelock;
end

if isfield(S.prep.clean.FTrejauto,'parallel_workers') && S.prep.clean.FTrejauto.parallel_workers
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj) && S.prep.clean.FTrejauto.parallel_workers
        parpool('local',S.prep.clean.FTrejauto.parallel_workers);
    elseif ~isempty(poolobj) && S.prep.clean.FTrejauto.parallel_workers
        if poolobj.NumWorkers ~= S.prep.clean.FTrejauto.parallel_workers
            delete(poolobj)
            parpool('local',S.prep.clean.FTrejauto.parallel_workers);
        end
    end
end

%% define chans
orig_chans=FT.label';
if ~exist('chan','var') || isempty(chan)
    chan = 1:length(orig_chans);
end
if any(chan>length(orig_chans))
    chan = 1:length(orig_chans);
end
channel = orig_chans(chan);

%% remove middle section of data
if any(S.(S.func).clean.FTrejauto.ignore_stim_artefact)
    art_time = S.(S.func).clean.FTrejauto.ignore_stim_artefact;
    cfg.latency = [FT.time{1, 1}(1) art_time(1)]; data1 = ft_selectdata(cfg, FT);
    cfg.latency = [art_time(2) FT.time{1, 1}(end)]; data2 = ft_selectdata(cfg, FT);
    temp = FT;
    for t = 1:length(FT.trial)
        temp.trial{t} = horzcat(data1.trial{t},data2.trial{t});
        temp.time{t} = horzcat(data1.time{t},data2.time{t});
    end
    FT=temp;
end

cutoffs=S.(S.func).clean.FTrejauto.cutoffs;

%% jump
cfg = [];
cfg.continuous = 'no';

% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel = channel; %'';
cfg.artfctdef.zvalue.cutoff = cutoffs.jump;
cfg.artfctdef.zvalue.trlpadding = 0;
cfg.artfctdef.zvalue.artpadding = 0;
cfg.artfctdef.zvalue.fltpadding = 0;

% algorithmic parameters
cfg.artfctdef.zvalue.cumulative = 'yes';
cfg.artfctdef.zvalue.medianfilter = 'yes';
cfg.artfctdef.zvalue.medianfiltord = 9;
cfg.artfctdef.zvalue.absdiff = 'yes';

% make the process interactive
cfg.artfctdef.zvalue.interactive = interactive;

% run
if isfield(S.prep.clean.FTrejauto,'parallel_workers') && S.prep.clean.FTrejauto.parallel_workers
    cfg.parallel = S.prep.clean.FTrejauto.parallel_workers;
    [cfg, artifact_jump] = ft_artifact_zvaluecab(cfg,FT);
else
    [cfg, artifact_jump] = ft_artifact_zvalue(cfg,FT);
end

%% muscle
cfg = [];
cfg.continuous = 'no';

% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel      = channel;
cfg.artfctdef.zvalue.cutoff       = cutoffs.muscle;
cfg.artfctdef.zvalue.trlpadding   = 0;
cfg.artfctdef.zvalue.fltpadding   = 0;
cfg.artfctdef.zvalue.artpadding   = 0.1;

% algorithmic parameters
upper=min(floor(EEG.srate/2.1),140);
cfg.artfctdef.zvalue.bpfilter     = 'yes';
cfg.artfctdef.zvalue.bpfreq       = [upper-30 upper];
cfg.artfctdef.zvalue.bpfiltord    = 9;
cfg.artfctdef.zvalue.bpfilttype   = 'but';
cfg.artfctdef.zvalue.hilbert      = 'yes';
cfg.artfctdef.zvalue.boxcar       = 0.2;

% make the process interactive
cfg.artfctdef.zvalue.interactive = interactive;
if isfield(S.prep.clean.FTrejauto,'parallel_workers') && S.prep.clean.FTrejauto.parallel_workers
    cfg.parallel = S.prep.clean.FTrejauto.parallel_workers;
    [cfg, artifact_muscle] = ft_artifact_zvaluecab(cfg,FT);
else
    [cfg, artifact_muscle] = ft_artifact_zvalue(cfg,FT);
end

% EOG
cfg = [];
cfg.continuous = 'no';

 % channel selection, cutoff and padding
 cfg.artfctdef.zvalue.channel     = channel;
 cfg.artfctdef.zvalue.cutoff      = cutoffs.EOG;
 cfg.artfctdef.zvalue.trlpadding  = 0;
 cfg.artfctdef.zvalue.artpadding  = 0.1;
 cfg.artfctdef.zvalue.fltpadding  = 0;

 % algorithmic parameters
 cfg.artfctdef.zvalue.bpfilter   = 'yes';
 cfg.artfctdef.zvalue.bpfilttype = 'but';
 cfg.artfctdef.zvalue.bpfreq     = [2 15];
 cfg.artfctdef.zvalue.bpfiltord  = 4;
 cfg.artfctdef.zvalue.hilbert    = 'yes';

 % feedback
 cfg.artfctdef.zvalue.interactive = interactive;

 if isfield(S.prep.clean.FTrejauto,'parallel_workers') && S.prep.clean.FTrejauto.parallel_workers
    cfg.parallel = S.prep.clean.FTrejauto.parallel_workers;
    [cfg, artifact_EOG] = ft_artifact_zvaluecab(cfg,FT);
 else
    [cfg, artifact_EOG] = ft_artifact_zvalue(cfg,FT);
 end


%cfg=[];
%cfg.artfctdef.reject = 'complete'; % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
%cfg.artfctdef.eog.artifact = artifact_EOG; %
%cfg.artfctdef.jump.artifact = artifact_jump;
%cfg.artfctdef.muscle.artifact = artifact_muscle;
%FT = ft_rejectartifact(cfg,FT);

%[~, rejchan] = setdiff(orig_chans(chan), FT.label);
trl = cfg.artfctdef.zvalue.trl(:,1:2);
rejtrial = zeros(length(trl),1);
for tr = 1:length(EEG.epoch)
    C(tr,1) = any(intersect(artifact_jump(:),trl(tr,1):trl(tr,2)));
    C(tr,2) = any(intersect(artifact_muscle(:),trl(tr,1):trl(tr,2)));
    C(tr,3) = any(intersect(artifact_EOG(:),trl(tr,1):trl(tr,2)));
    if any(C(tr,:))
        rejtrial(tr)=1;
    end
end

S.(S.func).clean.FTrejauto.rejtrial = find(rejtrial);
S.(S.func).clean.FTrejauto.Nrejtrial_by_type.jump = sum(C(:,1));
S.(S.func).clean.FTrejauto.Nrejtrial_by_type.muscle = sum(C(:,2));
S.(S.func).clean.FTrejauto.Nrejtrial_by_type.eog = sum(C(:,3));

disp(['TOTAL: ' num2str(sum(rejtrial)) ' out of ' num2str(length(rejtrial)) ' trials rejected'])
disp(['JUMP: ' num2str(sum(C(:,1))) ' out of ' num2str(length(rejtrial)) ' trials rejected'])
disp(['MUSCLE: ' num2str(sum(C(:,2))) ' out of ' num2str(length(rejtrial)) ' trials rejected'])
disp(['EOG: ' num2str(sum(C(:,3))) ' out of ' num2str(length(rejtrial)) ' trials rejected'])

% if strcmp(format,'EEGLAB')
%     %EEG = eeg_interp(EEG, rejchan);
%     %if length(EEG.chanlocs)>EEG.nbchan
%     %    EEG.chanlocs(rejchan)=[];
%     %end
%     EEG = pop_select(EEG, 'notrial', find(rejtrial));
% elseif strcmp(format,'SPM')
%     %EEG = badchannels(EEG, rejchan, 1);
%     EEG = badtrials(EEG, find(rejtrial), 1);
% end
