function S = noisyflat_channel_diagnostics(S,EEG,fname)

metric = struct; oi=0;
chan = S.(S.func).inclchan;

if ~S.(S.func).clean.noisychan.visual_mode
    S.fig =figure('Name',fname);
    tiledlayout(length(S.(S.func).clean.noisychan.metrics)*length(S.(S.func).clean.noisychan.filtering),4)
end

% filter x 2
for i = 1:length(S.(S.func).clean.noisychan.filtering)

    % convert to FT
    FT = convertoft(EEG);

    % general settings
    cfg =[];
    cfg.method = 'summary'; % calc for each channel and trial
    cfg.visual_mode = S.(S.func).clean.noisychan.visual_mode;
    cfg.alim = 1e-5;
    cfg.keepchannel='yes';
    cfg.trial='all';
    cfg.keeptrial='yes';

    % select chans
    orig_chans=FT.label';
    if ~exist('chan','var') || isempty(chan)
        chan = 1:length(orig_chans);
    end
    if any(chan>length(orig_chans))
        chan = 1:length(orig_chans);
    end
    cfg.channel = orig_chans(chan);


    % settings for filtering
    switch S.(S.func).clean.noisychan.filtering{i}
        case 'eog'
             %The following settings are useful for identifying EOG artifacts:
              cfg.preproc.bpfilter    = 'yes';
              cfg.preproc.bpfilttype  = 'but';
              cfg.preproc.bpfreq      = [1 15];
              cfg.preproc.bpfiltord   = 4;
              cfg.preproc.rectify     = 'yes';
        case 'muscle'
             %The following settings are useful for identifying muscle artifacts:
              cfg.preproc.bpfilter    = 'yes';
              cfg.preproc.bpfreq      = [100 120];
              cfg.preproc.bpfiltord   =  8;
              cfg.preproc.bpfilttype  = 'but';
              cfg.preproc.rectify     = 'yes';
              cfg.preproc.boxcar      = 0.2;
    end

    % settings for metrics
    %  for each channel in each trial, can be
    %  'var'       variance within each channel (default)
    %  'min'       minimum value in each channel
    %  'max'       maximum value each channel
    %  'maxabs'    maximum absolute value in each channel
    %  'range'     range from min to max in each channel
    %  'kurtosis'  kurtosis, i.e. measure of peakedness of the amplitude distribution
    %  'zvalue'    mean and std computed over all time and trials, per channel
    cfg.metric = S.(S.func).clean.noisychan.metrics;

    [FT,cfg] = ft_rejectvisual_cab(cfg, FT);

    if cfg.visual_mode
        [rejchanact] = setdiff(orig_chans(chan), FT.label);
        [~,rejchanfreq{i}] = intersect(orig_chans,rejchanact);
        [~, ~, rejtrialfreq{i}] = intersect(FT.cfg.artfctdef.summary.artifact(:,1),FT.sampleinfo(:,1)); 
    else
    % reject
        for m = 1:length(cfg.metric)
            oi=oi+1;
            name = [S.(S.func).clean.noisychan.filtering{i} ', ' cfg.metric{m}];
            level = cfg.level{m}; % chan x trials

%             % stats
%             mdn = median(level,'all');
%             OL = isoutlier(level,'quartiles');
%             OL_upper = OL.*[level>mdn];
%             OL_lower = OL.*[level<mdn];
%             drawmax = prctile(level(:),90);

            mdn = median(level,'all');
            if S.prep.clean.noisychan.gesd_ThresholdFactor
                OL = isoutlier(level,'gesd','ThresholdFactor',S.prep.clean.noisychan.gesd_ThresholdFactor);
            else
                OL = ones(size(level));
            end
            if strcmp(cfg.metric,'autocorr')
                OL_upper = OL.*[level<S.prep.clean.noisychan.metric_cutoffs.(S.(S.func).clean.noisychan.filtering{i})(m)];
            else
                OL_upper = OL.*[level>S.prep.clean.noisychan.metric_cutoffs.(S.(S.func).clean.noisychan.filtering{i})(m)];
                OL_lower = OL.*[level<-S.prep.clean.noisychan.metric_cutoffs.(S.(S.func).clean.noisychan.filtering{i})(m)];
            end
            drawmin = prctile(level(:),0);
            drawmax = prctile(level(:),95);

            if 1
                % plot
                nexttile
                boxplot(level(:));
                nexttile
                imagesc(level,[drawmin drawmax]); title(name)
                nexttile
                imagesc(OL_upper,[0 1]); title([name ': upper outlier'])
                nexttile
                imagesc(OL_lower,[0 1]); title([name ': lower outlier'])
            end

            % outputs
            var_name = strrep(name,', ', '_');
            metric(oi).name = var_name;
            metric(oi).var_trials_perchan = std(level,0,2);
            metric(oi).var_chan_pertrial = std(level,0,1);
            metric(oi).var_trial = std(metric(oi).var_trials_perchan);
            metric(oi).var_chan = std(metric(oi).var_chan_pertrial);
            metric(oi).OL_upper = OL_upper;
            metric(oi).OL_lower = OL_lower;

        end
    end
end

%reject
if cfg.visual_mode
    S.(S.func).clean.noisychan.rejchan = unique(vertcat(rejchanfreq{:}));
    S.(S.func).clean.noisychan.rejtrial = unique(vertcat(rejtrialfreq{:}));
else
    S.(S.func).clean.noisychan.metric = metric;
end


%% sub-functions

function [data,cfg] = ft_rejectvisual_cab(cfg, data)

% FT_REJECTVISUAL shows the preprocessed data in all channels and/or trials to
% allow the user to make a visual selection of the data that should be
% rejected. The data can be displayed in a "summary" mode, in which case
% the variance (or another metric) in each channel and each trial is
% computed. Alternatively, all channels can be shown at once allowing
% paging through the trials, or all trials can be shown, allowing paging
% through the channels.
%
% Use as
%   [data] = ft_rejectvisual(cfg, data)
%
% The configuration can contain
%   cfg.method      = string, describes how the data should be shown, this can be
%                     'summary'  show a single number for each channel and trial (default)
%                     'channel'  show the data per channel, all trials at once
%                     'trial'    show the data per trial, all channels at once
%   cfg.channel     = Nx1 cell-array with selection of channels (default = 'all'),
%                     see FT_CHANNELSELECTION for details
%   cfg.keepchannel = string, determines how to deal with channels that are not selected, can be
%                     'no'          completely remove deselected channels from the data (default)
%                     'yes'         keep deselected channels in the output data
%                     'nan'         fill the channels that are deselected with NaNs
%                     'repair'      repair the deselected channels using FT_CHANNELREPAIR
%   cfg.neighbours  = neighbourhood structure, see also FT_PREPARE_NEIGHBOURS (required for repairing channels)
%   cfg.trials      = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.keeptrial   = string, determines how to deal with trials that are
%                     not selected, can be
%                     'no'     completely remove deselected trials from the data (default)
%                     'yes'    keep deselected trials in the output data
%                     'nan'    fill the trials that are deselected with NaNs
%   cfg.metric      = string, describes the metric that should be computed in summary mode
%                     for each channel in each trial, can be
%                     'var'       variance within each channel (default)
%                     'min'       minimum value in each channel
%                     'max'       maximum value each channel
%                     'maxabs'    maximum absolute value in each channel
%                     'range'     range from min to max in each channel
%                     'kurtosis'  kurtosis, i.e. measure of peakedness of the amplitude distribution
%                     'zvalue'    mean and std computed over all time and trials, per channel
%   cfg.latency     = [begin end] in seconds, or 'minperlength', 'maxperlength',
%                     'prestim', 'poststim' (default = 'maxperlength')
%   cfg.alim        = value that determines the amplitude scaling for the
%                     channel and trial display, if empty then the amplitude
%                     scaling is automatic (default = [])
%   cfg.eegscale    = number, scaling to apply to the EEG channels prior to display
%   cfg.eogscale    = number, scaling to apply to the EOG channels prior to display
%   cfg.ecgscale    = number, scaling to apply to the ECG channels prior to display
%   cfg.emgscale    = number, scaling to apply to the EMG channels prior to display
%   cfg.megscale    = number, scaling to apply to the MEG channels prior to display
%   cfg.gradscale   = number, scaling to apply to the MEG gradiometer channels prior to display (in addition to the cfg.megscale factor)
%   cfg.magscale    = number, scaling to apply to the MEG magnetometer channels prior to display (in addition to the cfg.megscale factor)
%
% The scaling to the EEG, EOG, ECG, EMG and MEG channels is optional and can
% be used to bring the absolute numbers of the different channel types in
% the same range (e.g. fT and uV). The channel types are determined from
% the input data using FT_CHANNELSELECTION.
%
% Optionally, the raw data is preprocessed (filtering etc.) prior to
% displaying it or prior to computing the summary metric. The
% preprocessing and the selection of the latency window is NOT applied
% to the output data.
%
% The following settings are usefull for identifying EOG artifacts:
%   cfg.preproc.bpfilter    = 'yes'
%   cfg.preproc.bpfilttype  = 'but'
%   cfg.preproc.bpfreq      = [1 15]
%   cfg.preproc.bpfiltord   = 4
%   cfg.preproc.rectify     = 'yes'
%
% The following settings are usefull for identifying muscle artifacts:
%   cfg.preproc.bpfilter    = 'yes'
%   cfg.preproc.bpfreq      = [110 140]
%   cfg.preproc.bpfiltord   =  8
%   cfg.preproc.bpfilttype  = 'but'
%   cfg.preproc.rectify     = 'yes'
%   cfg.preproc.boxcar      = 0.2
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_REJECTARTIFACT, FT_REJECTCOMPONENT

% Undocumented local options:
% cfg.feedback

% Copyright (C) 2005-2006, Markus Bauer, Robert Oostenveld
% Copyright (C) 2006-2016, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% Undocumented options
% cfg.plotlayout = 'square' (default) or '1col', plotting every channel/trial under each other
% cfg.viewmode   = 'remove', 'toggle' or 'hide', only applies to summary mode (default = 'remove')

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% store original datatype
dtype = ft_datatype(data);

% check if the input data is valid for this function, this will convert it to raw if needed
data = ft_checkdata(data, 'datatype', {'raw+comp', 'raw'}, 'feedback', 'yes', 'hassampleinfo', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamedval',  {'metric',  'absmax',  'maxabs'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'method',  'absmax',  'maxabs'});

% resolve some common typing errors
cfg = ft_checkconfig(cfg, 'renamed',  {'keeptrials',  'keeptrial'});
cfg = ft_checkconfig(cfg, 'renamed',  {'keepchannels',  'keepchannel'});

% set the defaults
cfg.channel     = ft_getopt(cfg, 'channel'    , 'all');
cfg.trials      = ft_getopt(cfg, 'trials'     , 'all', 1);
cfg.latency     = ft_getopt(cfg, 'latency'    , 'maxperlength');
cfg.keepchannel = ft_getopt(cfg, 'keepchannel', 'no');
cfg.keeptrial   = ft_getopt(cfg, 'keeptrial'  , 'no');
cfg.feedback    = ft_getopt(cfg, 'feedback'   , 'textbar');
cfg.method      = ft_getopt(cfg, 'method'     , 'summary');
cfg.metric      = ft_getopt(cfg, 'metric'     , 'var');
cfg.alim        = ft_getopt(cfg, 'alim'       );
cfg.eegscale    = ft_getopt(cfg, 'eegscale'   );
cfg.eogscale    = ft_getopt(cfg, 'eogscale'   );
cfg.ecgscale    = ft_getopt(cfg, 'ecgscale'   );
cfg.emgscale    = ft_getopt(cfg, 'emgscale'   );
cfg.megscale    = ft_getopt(cfg, 'megscale'   );
cfg.gradscale   = ft_getopt(cfg, 'gradscale'  );
cfg.magscale    = ft_getopt(cfg, 'magscale'   );
cfg.plotlayout  = ft_getopt(cfg, 'plotlayout' , 'square');
cfg.viewmode    = ft_getopt(cfg, 'viewmode'   , 'remove');

% ensure that the preproc specific options are located in the cfg.preproc substructure
cfg = ft_checkconfig(cfg, 'createsubcfg',  {'preproc'});

% check required fields at the start, rather than further down in the code
if strcmp(cfg.keepchannel, 'repair')
  cfg = ft_checkconfig(cfg, 'required', 'neighbours');
end

% determine the duration of each trial
for i=1:length(data.time)
  begsamplatency(i) = min(data.time{i});
  endsamplatency(i) = max(data.time{i});
end

% determine the latency window which is possible in all trials
minperlength = [max(begsamplatency) min(endsamplatency)];
maxperlength = [min(begsamplatency) max(endsamplatency)];

% latency window for averaging and variance computation is given in seconds
if (strcmp(cfg.latency, 'minperlength'))
  cfg.latency = [];
  cfg.latency(1) = minperlength(1);
  cfg.latency(2) = minperlength(2);
elseif (strcmp(cfg.latency, 'maxperlength'))
  cfg.latency = [];
  cfg.latency(1) = maxperlength(1);
  cfg.latency(2) = maxperlength(2);
elseif (strcmp(cfg.latency, 'prestim'))
  cfg.latency = [];
  cfg.latency(1) = maxperlength(1);
  cfg.latency(2) = 0;
elseif (strcmp(cfg.latency, 'poststim'))
  cfg.latency = [];
  cfg.latency(1) = 0;
  cfg.latency(2) = maxperlength(2);
end

% apply scaling to the selected channel types to equate the absolute numbers (i.e. fT and uV)
% make a seperate copy to prevent the original data from being scaled
tmpdata = data;
scaled  = false;
if ~isempty(cfg.eegscale)
  scaled = true;
  chansel = match_str(tmpdata.label, ft_channelselection('EEG', tmpdata.label));
  for i=1:length(tmpdata.trial)
    tmpdata.trial{i}(chansel,:) = tmpdata.trial{i}(chansel,:) .* cfg.eegscale;
  end
end
if ~isempty(cfg.eogscale)
  scaled = true;
  chansel = match_str(tmpdata.label, ft_channelselection('EOG', tmpdata.label));
  for i=1:length(tmpdata.trial)
    tmpdata.trial{i}(chansel,:) = tmpdata.trial{i}(chansel,:) .* cfg.eogscale;
  end
end
if ~isempty(cfg.ecgscale)
  scaled = true;
  chansel = match_str(tmpdata.label, ft_channelselection('ECG', tmpdata.label));
  for i=1:length(tmpdata.trial)
    tmpdata.trial{i}(chansel,:) = tmpdata.trial{i}(chansel,:) .* cfg.ecgscale;
  end
end
if ~isempty(cfg.emgscale)
  scaled = true;
  chansel = match_str(tmpdata.label, ft_channelselection('EMG', tmpdata.label));
  for i=1:length(tmpdata.trial)
    tmpdata.trial{i}(chansel,:) = tmpdata.trial{i}(chansel,:) .* cfg.emgscale;
  end
end
if ~isempty(cfg.megscale)
  scaled = true;
  chansel = match_str(tmpdata.label, ft_channelselection('MEG', tmpdata.label));
  for i=1:length(tmpdata.trial)
    tmpdata.trial{i}(chansel,:) = tmpdata.trial{i}(chansel,:) .* cfg.megscale;
  end
end
if ~isempty(cfg.gradscale)
  scaled = true;
  chansel = match_str(tmpdata.label, ft_channelselection('MEGGRAD', tmpdata.label));
  for i=1:length(tmpdata.trial)
    tmpdata.trial{i}(chansel,:) = tmpdata.trial{i}(chansel,:) .* cfg.gradscale;
  end
end
if ~isempty(cfg.magscale)
  scaled = true;
  chansel = match_str(tmpdata.label, ft_channelselection('MEGMAG', tmpdata.label));
  for i=1:length(tmpdata.trial)
    tmpdata.trial{i}(chansel,:) = tmpdata.trial{i}(chansel,:) .* cfg.magscale;
  end
end

switch cfg.method
  case 'channel'
    if scaled
      fprintf('showing the scaled data per channel, all trials at once\n');
    else
      fprintf('showing the data per channel, all trials at once\n');
    end
    [chansel, trlsel, cfg] = rejectvisual_channel_cab(cfg, tmpdata);

  case 'trial'
    if scaled
      fprintf('showing the scaled per trial, all channels at once\n');
    else
      fprintf('showing the data per trial, all channels at once\n');
    end
    [chansel, trlsel, cfg] = rejectvisual_trial_cab(cfg, tmpdata);

  case 'summary'
    if scaled
      fprintf('showing a summary of the scaled data for all channels and trials\n');
    else
      fprintf('computing a summary of the data for all channels and trials\n');
    end
    [chansel, trlsel, cfg] = rejectvisual_summary_cab(cfg, tmpdata);

  otherwise
    ft_error('unsupported method %s', cfg.method);
end % switch method

fprintf('%d trials marked as GOOD, %d trials marked as BAD\n', sum(trlsel), sum(~trlsel));
fprintf('%d channels marked as GOOD, %d channels marked as BAD\n', sum(chansel), sum(~chansel));


if ~all(chansel)
  switch cfg.keepchannel
    case 'yes'
      % keep all channels, also when they are not selected
      fprintf('no channels were removed from the data\n');

    case 'no'
      % show the user which channels are removed
      removed = find(~chansel);
      fprintf('the following channels were removed: ');
      for i=1:(length(removed)-1)
        fprintf('%s, ', data.label{removed(i)});
      end
      fprintf('%s\n', data.label{removed(end)});

    case 'nan'
      % show the user which channels are removed
      removed = find(~chansel);
      fprintf('the following channels were filled with NANs: ');
      for i=1:(length(removed)-1)
        fprintf('%s, ', data.label{removed(i)});
      end
      fprintf('%s\n', data.label{removed(end)});
      % mark the selection as nan
      for i=1:length(data.trial)
        data.trial{i}(~chansel,:) = nan;
      end

    case 'repair'
      % show which channels are to be repaired
      removed = find(~chansel);
      fprintf('the following channels were repaired using FT_CHANNELREPAIR: ');
      for i=1:(length(removed)-1)
        fprintf('%s, ', data.label{removed(i)});
      end
      fprintf('%s\n', data.label{removed(end)});

      % create cfg struct for call to ft_channelrepair
      orgcfg = cfg;
      tmpcfg = [];
      tmpcfg.trials = 'all';
      tmpcfg.badchannel = data.label(~chansel);
      tmpcfg.neighbours = cfg.neighbours;
      if isfield(cfg, 'grad')
          tmpcfg.grad = cfg.grad;
      end
      if isfield(cfg, 'elec')
          tmpcfg.elec = cfg.elec;
      end
      % repair the channels that were selected as bad
      data = ft_channelrepair(tmpcfg, data);
      % restore the provenance information
      [cfg, data] = rollback_provenance(cfg, data);
      % restore the original trials parameter, it should not be 'all'
      cfg = copyfields(orgcfg, cfg, {'trials'});

    otherwise
      ft_error('invalid specification of cfg.keepchannel')
  end % case
end % if ~all(chansel)

if ~all(trlsel)
  switch cfg.keeptrial
    case 'yes'
      % keep all trials, also when they are not selected
      fprintf('no trials were removed from the data\n');

    case 'no'
      % show the user which channels are removed
      removed = find(~trlsel);
      fprintf('the following trials were removed: ');
      for i=1:(length(removed)-1)
        fprintf('%d, ', removed(i));
      end
      fprintf('%d\n', removed(end));

    case 'nan'
      % show the user which trials are removed
      removed = find(~trlsel);
      fprintf('the following trials were filled with NANs: ');
      for i=1:(length(removed)-1)
        fprintf('%d, ', removed(i));
      end
      fprintf('%d\n', removed(end));
      % mark the selection as nan
      for i=removed
        data.trial{i}(:,:) = nan;
      end

    otherwise
      ft_error('invalid specification of cfg.keeptrial')
  end % case
end % if ~all(trlsel)

if isfield(data, 'sampleinfo')
  % construct the matrix with sample numbers prior to making the selection
  cfg.artfctdef.(cfg.method).artifact = data.sampleinfo(~trlsel,:);
end

% perform the selection of channels and trials
orgcfg = cfg;
tmpcfg = [];
if strcmp(cfg.keepchannel, 'no')
  tmpcfg.channel = find(chansel);
end
if strcmp(cfg.keeptrial, 'no')
  tmpcfg.trials = find(trlsel); % note that it is keeptrial without S and trials with S
end
data = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);
% restore the original channels and trials parameters
cfg = copyfields(orgcfg, cfg, {'channel', 'trials'});

% convert back to input type if necessary
switch dtype
  case 'timelock'
    data = ft_checkdata(data, 'datatype', 'timelock');
  otherwise
    % keep the output as it is
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   data
ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data

function [chansel, trlsel, cfg] = rejectvisual_summary_cab(cfg, data)

% SUBFUNCTION for ft_rejectvisual

% determine the initial selection of trials
ntrl = length(data.trial);
if isequal(cfg.trials, 'all') % support specification like 'all'
  cfg.trials = 1:ntrl;
end
trlsel = false(1, ntrl);
trlsel(cfg.trials) = true;

% determine the initial selection of channels
nchan = length(data.label);
cfg.channel = ft_channelselection(cfg.channel, data.label); % support specification like 'all'
chansel = false(1, nchan);
chansel(match_str(data.label, cfg.channel)) = true;

% compute the sampling frequency from the first two timepoints
fsample = 1/mean(diff(data.time{1}));

% select the specified latency window from the data
% here it is done BEFORE filtering and metric computation
for i=1:ntrl
  begsample = nearest(data.time{i}, cfg.latency(1));
  endsample = nearest(data.time{i}, cfg.latency(2));
  data.time{i} = data.time{i}(begsample:endsample);
  data.trial{i} = data.trial{i}(:, begsample:endsample);
end

% compute the offset from the time axes
offset = zeros(ntrl, 1);
for i=1:ntrl
  offset(i) = time2offset(data.time{i}, fsample);
end

% set up info
info                = [];
info.data           = data;
info.cfg            = cfg;
info.previousmetric = 'none';
info.level          = nan(nchan, ntrl);
info.ntrl           = ntrl;
info.nchan          = nchan;
info.trlsel         = trlsel;
info.chansel        = chansel;
info.fsample        = fsample;
info.offset         = offset;
info.quit           = 0;

if ~cfg.visual_mode
    
    % Compute all metrics
    for i = 1:length(cfg.metric)
        info.metric = cfg.metric{i};
        info = compute_metric_auto(info);
        cfg.level{i} = info.level;
    end

else
    info.metric = cfg.metric{1}; %cab
    h = figure();
    guidata(h, info);
    
    % set up display
    interactive = true;
    
    % make the figure large enough to hold stuff
    set(h, 'Position', [50 350 800 500]);
    
    % define three axes
    info.axes(1) = axes('position', [0.100 0.650 0.375 0.300]);  % summary
    info.axes(2) = axes('position', [0.575 0.650 0.375 0.300]);  % channels
    info.axes(3) = axes('position', [0.100 0.250 0.375 0.300]);  % trials
    % callback function (for toggling trials/channels) is set later, so that
    % the user cannot try to toggle trials while nothing is present on the
    % plots
    
    % instructions
    instructions = sprintf('Drag the mouse over the channels or trials you wish to reject');
    uicontrol(h, 'Units', 'normalized', 'position', [0.520 0.520 0.400 0.050], 'Style', 'text', 'HorizontalAlignment', 'left', 'backgroundcolor', get(h, 'color'), 'string', instructions, 'FontWeight', 'bold', 'ForegroundColor', 'r');
    
    % set up radio buttons for choosing metric
    bgcolor = get(h, 'color');
    g = uibuttongroup('Position', [0.520 0.220 0.375 0.250 ], 'bordertype', 'none', 'backgroundcolor', bgcolor);
    r(1) = uicontrol('Units', 'normalized', 'parent', g, 'position', [ 0.0  8/9 0.40 0.15 ], 'Style', 'Radio', 'backgroundcolor', bgcolor, 'string', 'var',       'HandleVisibility', 'off');
    r(2) = uicontrol('Units', 'normalized', 'parent', g, 'position', [ 0.0  7/9 0.40 0.15 ], 'Style', 'Radio', 'backgroundcolor', bgcolor, 'String', 'min',       'HandleVisibility', 'off');
    r(3) = uicontrol('Units', 'normalized', 'parent', g, 'position', [ 0.0  6/9 0.40 0.15 ], 'Style', 'Radio', 'backgroundcolor', bgcolor, 'String', 'max',       'HandleVisibility', 'off');
    r(4) = uicontrol('Units', 'normalized', 'parent', g, 'position', [ 0.0  5/9 0.40 0.15 ], 'Style', 'Radio', 'backgroundcolor', bgcolor, 'String', 'maxabs',    'HandleVisibility', 'off');
    r(5) = uicontrol('Units', 'normalized', 'parent', g, 'position', [ 0.0  4/9 0.40 0.15 ], 'Style', 'Radio', 'backgroundcolor', bgcolor, 'String', 'range',     'HandleVisibility', 'off');
    r(6) = uicontrol('Units', 'normalized', 'parent', g, 'position', [ 0.0  3/9 0.40 0.15 ], 'Style', 'Radio', 'backgroundcolor', bgcolor, 'String', 'kurtosis',  'HandleVisibility', 'off');
    r(7) = uicontrol('Units', 'normalized', 'parent', g, 'position', [ 0.0  2/9 0.40 0.15 ], 'Style', 'Radio', 'backgroundcolor', bgcolor, 'String', '1/var',     'HandleVisibility', 'off');
    r(8) = uicontrol('Units', 'normalized', 'parent', g, 'position', [ 0.0  1/9 0.40 0.15 ], 'Style', 'Radio', 'backgroundcolor', bgcolor, 'String', 'zvalue',    'HandleVisibility', 'off');
    r(9) = uicontrol('Units', 'normalized', 'parent', g, 'position', [ 0.0  0/9 0.40 0.15 ], 'Style', 'Radio', 'backgroundcolor', bgcolor, 'String', 'maxzvalue', 'HandleVisibility', 'off');


    % pre-select appropriate metric, if defined
    set(g, 'SelectionChangeFcn', @change_metric);
    for i=1:length(r)
      if strcmp(get(r(i), 'string'), cfg.metric)
        set(g, 'SelectedObject', r(i));
      end
    end


    % editboxes for manually specifying which channels/trials to toggle on or off
    uicontrol(h, 'Units', 'normalized', 'position', [0.630 0.470 0.14 0.05], 'Style', 'text', 'HorizontalAlignment', 'left', 'backgroundcolor', get(h, 'color'), 'string', 'Toggle trial:');
    uicontrol(h, 'Units', 'normalized', 'position', [0.630 0.430 0.14 0.05], 'Style', 'edit', 'HorizontalAlignment', 'left', 'backgroundcolor', [1 1 1], 'callback', @toggle_trials);
    uicontrol(h, 'Units', 'normalized', 'position', [0.630 0.340 0.14 0.05], 'Style', 'text', 'HorizontalAlignment', 'left', 'backgroundcolor', get(h, 'color'), 'string', 'Toggle channel:');
    uicontrol(h, 'Units', 'normalized', 'position', [0.630 0.300 0.14 0.05], 'Style', 'edit', 'HorizontalAlignment', 'left', 'backgroundcolor', [1 1 1], 'callback', @toggle_channels);
    
    % editbox for trial plotting
    uicontrol(h, 'Units', 'normalized', 'position', [0.630 0.210 0.14 0.05], 'Style', 'text', 'HorizontalAlignment', 'left', 'backgroundcolor', get(h, 'color'), 'string', 'Plot trial:');
    info.plottrltxt = uicontrol(h, 'Units', 'normalized', 'position', [0.630 0.170 0.14 0.05], 'Style', 'edit', 'HorizontalAlignment', 'left', 'backgroundcolor', [1 1 1], 'callback', @display_trial);
    
    info.badtrllbl  = uicontrol(h, 'Units', 'normalized', 'position', [0.795 0.470 0.230 0.05], 'Style', 'text', 'HorizontalAlignment', 'left', 'backgroundcolor', get(h, 'color'), 'string', sprintf('Rejected trials: %i/%i', sum(info.trlsel==0), info.ntrl));
    info.badtrltxt  = uicontrol(h, 'Units', 'normalized', 'position', [0.795 0.430 0.230 0.05], 'Style', 'text', 'HorizontalAlignment', 'left', 'backgroundcolor', get(h, 'color'));
    info.badchanlbl = uicontrol(h, 'Units', 'normalized', 'position', [0.795 0.340 0.230 0.05], 'Style', 'text', 'HorizontalAlignment', 'left', 'backgroundcolor', get(h, 'color'), 'string', sprintf('Rejected channels: %i/%i', sum(info.chansel==0), info.nchan));
    info.badchantxt = uicontrol(h, 'Units', 'normalized', 'position', [0.795 0.300 0.230 0.05], 'Style', 'text', 'HorizontalAlignment', 'left', 'backgroundcolor', get(h, 'color'));
    
    
    % "show rejected" button
    % ui_tog = uicontrol(h, 'Units', 'normalized', 'position', [0.55 0.200 0.25 0.05], 'Style', 'checkbox', 'backgroundcolor', get(h, 'color'), 'string', 'Show rejected?', 'callback', @toggle_rejected);
    % if strcmp(cfg.viewmode, 'toggle')
    %     set(ui_tog, 'value', 1);
    % end
    
    % logbox
    info.output_box = uicontrol(h, 'Units', 'normalized', 'position', [0.00 0.00 1.00 0.15], 'Style', 'edit', 'HorizontalAlignment', 'left', 'Max', 3, 'Min', 1, 'Enable', 'inactive', 'FontName', get(0, 'FixedWidthFontName'), 'FontSize', 9, 'ForegroundColor', [0 0 0], 'BackgroundColor', [1 1 1]);
    
    % quit button
    uicontrol(h, 'Units', 'normalized', 'position', [0.80 0.175 0.10 0.05], 'string', 'quit', 'callback', @quit);
    
    guidata(h, info);
    
    % disable trial plotting if cfg.layout not present
    if ~isfield(info.cfg, 'layout')
      set(info.plottrltxt, 'Enable', 'off');
      update_log(info.output_box, sprintf('NOTE: "cfg.layout" parameter required for trial plotting!'));
    end

    % Compute initial metric...
    compute_metric(h);

    while interactive && ishandle(h)
      redraw(h);
      info = guidata(h);
      if info.quit == 0
        uiwait(h);
      else
        chansel = info.chansel;
        trlsel  = info.trlsel;
        cfg     = info.cfg;
        delete(h);
        break;
      end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function info = compute_metric_auto(info)

tmp = info.level(info.chansel, info.trlsel);

if strcmp(info.metric, 'zvalue') || strcmp(info.metric, 'maxzvalue')
  % cellmean and cellstd (see ft_denoise_pca) would work instead of for-loops, but they are too memory-intensive
  runsum = zeros(info.nchan, 1);
  runss  = zeros(info.nchan, 1);
  runnum = 0;

  tic
  tempdat=cell(1,info.ntrl);
  poolobj = gcp('nocreate');
  if ~isempty(poolobj)
      parfor i=1:info.ntrl
        tempdat{i} = preproc(info.data.trial{i}, info.data.label, offset2time(info.offset(i), info.fsample, size(info.data.trial{i}, 2)), info.cfg.preproc); % not entirely sure whether info.data.time{i} is correct, so making it on the fly
      end
  else
      for i=1:info.ntrl
        tempdat{i} = preproc(info.data.trial{i}, info.data.label, offset2time(info.offset(i), info.fsample, size(info.data.trial{i}, 2)), info.cfg.preproc); % not entirely sure whether info.data.time{i} is correct, so making it on the fly
      end
  end
  toc

  for i=1:info.ntrl
    runsum = runsum + nansum(tempdat{i}, 2);
    runss  = runss  + nansum(tempdat{i}.^2, 2);
    runnum = runnum + sum(isfinite(tempdat{i}), 2);
  end

  mval = runsum./runnum;
  sd   = sqrt(runss./runnum - (runsum./runnum).^2);
else
  % otherwise next parfor fails
  mval = zeros(info.nchan,1);
  sd   = ones(info.nchan,1);
end
level = cell(1,info.ntrl);
tic
poolobj = gcp('nocreate');
  if ~isempty(poolobj)
    parfor i=1:info.ntrl
      dat = preproc(info.data.trial{i}, info.data.label, offset2time(info.offset(i), info.fsample, size(info.data.trial{i}, 2)), info.cfg.preproc); % not entirely sure whether info.data.time{i} is correct, so making it on the fly
      switch info.metric
        case 'var'
          level{i} = nanstd(dat, [], 2).^2;
        case 'min'
          level{i} = nanmin(dat, [], 2);
        case 'max'
          level{i} = nanmax(dat, [], 2);
        case 'maxabs'
          level{i} = nanmax(abs(dat), [], 2);
        case 'range'
          level{i} = nanmax(dat, [], 2) - nanmin(dat, [], 2);
        case 'kurtosis'
          level{i} = kurtosis(dat, [], 2);
        case '1/var'
          level{i} = 1./(nanstd(dat, [], 2).^2);
        case 'zvalue'
          level{i} = nanmean( (dat-repmat(mval, 1, size(dat, 2)) )./repmat(sd, 1, size(dat, 2)) , 2);
        case 'maxzvalue'
          level{i} = nanmax( ( dat-repmat(mval, 1, size(dat, 2)) )./repmat(sd, 1, size(dat, 2)) , [], 2);
        case 'autocorr'
            Ncorrint=round(20/(1000/info.fsample)); % number of samples for lag
            for k=1:size(dat,1)
                yy=xcorr(dat(k,:),Ncorrint,'coeff');
                level{i}(k,1) = yy(1);
            end
        otherwise
          ft_error('unsupported method');
      end
    end
  else
    for i=1:info.ntrl
      dat = preproc(info.data.trial{i}, info.data.label, offset2time(info.offset(i), info.fsample, size(info.data.trial{i}, 2)), info.cfg.preproc); % not entirely sure whether info.data.time{i} is correct, so making it on the fly
      switch info.metric
        case 'var'
          level{i} = nanstd(dat, [], 2).^2;
        case 'min'
          level{i} = nanmin(dat, [], 2);
        case 'max'
          level{i} = nanmax(dat, [], 2);
        case 'maxabs'
          level{i} = nanmax(abs(dat), [], 2);
        case 'range'
          level{i} = nanmax(dat, [], 2) - nanmin(dat, [], 2);
        case 'kurtosis'
          level{i} = kurtosis(dat, [], 2);
        case '1/var'
          level{i} = 1./(nanstd(dat, [], 2).^2);
        case 'zvalue'
          level{i} = nanmean( (dat-repmat(mval, 1, size(dat, 2)) )./repmat(sd, 1, size(dat, 2)) , 2);
        case 'maxzvalue'
          level{i} = nanmax( ( dat-repmat(mval, 1, size(dat, 2)) )./repmat(sd, 1, size(dat, 2)) , [], 2);
        case 'autocorr'
            Ncorrint=round(20/(1000/info.fsample)); % number of samples for lag
            for k=1:size(dat,1)
                yy=xcorr(dat(k,:),Ncorrint,'coeff');
                level{i}(k,1) = yy(1);
            end
        otherwise
          ft_error('unsupported method');
      end
    end
  end
toc
info.level = horzcat(level{:});


function compute_metric(h)
info = guidata(h);
tmp = info.level(info.chansel, info.trlsel);
if isequal(info.metric, info.previousmetric) && all(~isnan(tmp(:)))
  % there is no reason to recompute the metric
  return
end

update_log(info.output_box, 'Computing metric...');
ft_progress('init', info.cfg.feedback, 'computing metric');
level = zeros(info.nchan, info.ntrl);
if strcmp(info.metric, 'zvalue') || strcmp(info.metric, 'maxzvalue')
  % cellmean and cellstd (see ft_denoise_pca) would work instead of for-loops, but they are too memory-intensive
  runsum = zeros(info.nchan, 1);
  runss  = zeros(info.nchan, 1);
  runnum = 0;
  for i=1:info.ntrl
    dat = preproc(info.data.trial{i}, info.data.label, offset2time(info.offset(i), info.fsample, size(info.data.trial{i}, 2)), info.cfg.preproc); % not entirely sure whether info.data.time{i} is correct, so making it on the fly
    runsum = runsum + nansum(dat, 2);
    runss  = runss  + nansum(dat.^2, 2);
    runnum = runnum + sum(isfinite(dat), 2);
  end
  mval = runsum./runnum;
  sd   = sqrt(runss./runnum - (runsum./runnum).^2);
end
tic
for i=1:info.ntrl
  ft_progress(i/info.ntrl, 'computing metric %d of %d\n', i, info.ntrl);
  dat = preproc(info.data.trial{i}, info.data.label, offset2time(info.offset(i), info.fsample, size(info.data.trial{i}, 2)), info.cfg.preproc); % not entirely sure whether info.data.time{i} is correct, so making it on the fly
  switch info.metric
    case 'var'
      level(:, i) = nanstd(dat, [], 2).^2;
    case 'min'
      level(:, i) = nanmin(dat, [], 2);
    case 'max'
      level(:, i) = nanmax(dat, [], 2);
    case 'maxabs'
      level(:, i) = nanmax(abs(dat), [], 2);
    case 'range'
      level(:, i) = nanmax(dat, [], 2) - nanmin(dat, [], 2);
    case 'kurtosis'
      level(:, i) = kurtosis(dat, [], 2);
    case '1/var'
      level(:, i) = 1./(nanstd(dat, [], 2).^2);
    case 'zvalue'
      level(:, i) = nanmean( (dat-repmat(mval, 1, size(dat, 2)) )./repmat(sd, 1, size(dat, 2)) , 2);
    case 'maxzvalue'
      level(:, i) = nanmax( ( dat-repmat(mval, 1, size(dat, 2)) )./repmat(sd, 1, size(dat, 2)) , [], 2);
    case 'autocorr'
        Ncorrint=round(20/(1000/info.fsample)); % number of samples for lag
        for k=1:size(dat,1)
            yy=xcorr(dat(k,:),Ncorrint,'coeff');
            level(k,i) = yy(1);
        end
    otherwise
      ft_error('unsupported method');
  end
end
toc
ft_progress('close');
update_log(info.output_box, 'Done.');
info.level = level;
info.previousmetric = info.metric;
guidata(h, info);

function redraw(h)
info  = guidata(h);
% work with a copy of the data
level = info.level;

[maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = set_maxper(level, info.chansel, info.trlsel, strcmp(info.metric, 'min'));

% make the three figures
if gcf~=h, figure(h); end
% datacursormode on;

set(h, 'CurrentAxes', info.axes(1))
cla(info.axes(1));
switch info.cfg.viewmode
  case {'remove', 'toggle'}
    tmp = level;
    tmp(~info.chansel, :) = nan;
    tmp(:, ~info.trlsel)  = nan;
    imagesc(tmp);
  case 'hide'
    imagesc(level(info.chansel==1, info.trlsel==1));
    if ~all(info.trlsel)
      set(info.axes(1), 'Xtick', []);
    end
    if ~all(info.chansel)
      set(info.axes(1), 'Ytick', []);
    end
end % switch
axis ij;
% colorbar;
title(info.cfg.method);
ylabel('channel number');
xlabel('trial number');

set(h, 'CurrentAxes', info.axes(2))
cla(info.axes(2));
switch info.cfg.viewmode
  case 'remove'
    plot(maxperchan(info.chansel==1),     find(info.chansel==1), '.');
    xmax = max(maxperchan);
    xmin = min(maxperchan);
    ymax = info.nchan;
  case 'toggle'
    plot(maxperchan_all(info.chansel==1), find(info.chansel==1), '.');
    hold on;
    plot(maxperchan_all(info.chansel==0), find(info.chansel==0), 'o');
    hold off;
    xmax = max(maxperchan_all);
    xmin = min(maxperchan_all);
    ymax = info.nchan;
  case 'hide'
    xmax = max(maxperchan);
    xmin = min(maxperchan);
    ymax = sum(info.chansel==1);
    plot(maxperchan(info.chansel==1), 1:ymax, '.');
    if ~all(info.chansel)
      set(info.axes(2), 'Ytick', []);
    end
end % switch
if any(info.chansel) && any(info.trlsel)
  % don't try to rescale the axes if they are empty
  % have to use 0 as lower limit because in the single channel case ylim([1 1]) will be invalid
  axis([0.8*xmin 1.2*xmax 0.5 ymax+0.5]);
end
axis ij;
set(info.axes(2), 'ButtonDownFcn', @toggle_visual);  % needs to be here; call to axis resets this property
ylabel('channel number');

set(h, 'CurrentAxes', info.axes(3))
cla(info.axes(3));
switch info.cfg.viewmode
  case 'remove'
    plot(find(info.trlsel==1), maxpertrl(info.trlsel==1), '.');
    xmax = info.ntrl;
    ymax = max(maxpertrl);
    ymin = min(maxpertrl);
  case 'toggle'
    plot(find(info.trlsel==1), maxpertrl_all(info.trlsel==1), '.');
    hold on;
    plot(find(info.trlsel==0), maxpertrl_all(info.trlsel==0), 'o');
    hold off;
    xmax = info.ntrl;
    ymax = max(maxpertrl_all);
    ymin = min(maxpertrl_all);
  case 'hide'
    xmax = sum(info.trlsel==1);
    ymax = max(maxpertrl);
    ymin = min(maxpertrl);
    plot(1:xmax, maxpertrl(info.trlsel==1), '.');
    if ~all(info.trlsel)
      set(info.axes(3), 'Xtick', []);
    end
end % switch
if any(info.chansel) && any(info.trlsel)
  % don't try to rescale the axes if they are empty
  % the 0.8-1.2 is needed to deal with the single trial case
  % note that both ymin and ymax can be negative
  axis([0.5 xmax+0.5 (1-sign(ymin)*0.2)*ymin (1+sign(ymax)*0.2)*ymax]);
end
set(info.axes(3), 'ButtonDownFcn', @toggle_visual);  % needs to be here; call to axis resets this property
xlabel('trial number');

% put rejected trials/channels in their respective edit boxes
set(info.badchanlbl, 'string', sprintf('Rejected channels: %i/%i', sum(info.chansel==0), info.nchan));
set(info.badtrllbl, 'string', sprintf('Rejected trials: %i/%i', sum(info.trlsel==0), info.ntrl));
if ~isempty(find(info.trlsel==0, 1))
  set(info.badtrltxt, 'String', num2str(find(info.trlsel==0)), 'FontAngle', 'normal');
else
  set(info.badtrltxt, 'String', 'No trials rejected', 'FontAngle', 'italic');
end
if ~isempty(find(info.chansel==0, 1))
  if isfield(info.data, 'label')
    chanlabels = info.data.label(info.chansel==0);
    badchantxt = '';
    for i=find(info.chansel==0)
      if ~isempty(badchantxt)
        badchantxt = [badchantxt ', ' info.data.label{i} '(' num2str(i) ')'];
      else
        badchantxt = [info.data.label{i} '(' num2str(i) ')'];
      end
    end
    set(info.badchantxt, 'String', badchantxt, 'FontAngle', 'normal');
  else
    set(info.badtrltxt, 'String', num2str(find(info.chansel==0)), 'FontAngle', 'normal');
  end
else
  set(info.badchantxt, 'String', 'No channels rejected', 'FontAngle', 'italic');
end

function toggle_trials(h, eventdata)
info = guidata(h);
% extract trials from string
rawtrls = get(h, 'string');
if ~isempty(rawtrls)
  spltrls = regexp(rawtrls, '\s+', 'split');
  trls = [];
  for n = 1:length(spltrls)
    trls(n) = str2num(cell2mat(spltrls(n)));
  end
else
  update_log(info.output_box, sprintf('Please enter one or more trials'));
  uiresume;
  return;
end
try
  toggle = trls;
  info.trlsel(toggle) = ~info.trlsel(toggle);
catch
  update_log(info.output_box, sprintf('ERROR: Trial value too large!'));
  uiresume;
  return;
end
% recalculate the metric
compute_metric(h)
guidata(h, info);
uiresume;

% process input from the "toggle channels" textbox
function toggle_channels(h, eventdata)
info = guidata(h);
rawchans = get(h, 'string');
if ~isempty(rawchans)
  splchans = regexp(rawchans, '\s+', 'split');
  chans = zeros(1, length(splchans));
  % determine whether identifying channels via number or label
  [junk, junk, junk, procchans] = regexp(rawchans, '([A-Za-z]+|[0-9]{4, })');
  clear junk;
  if isempty(procchans)
    % if using channel numbers
    for n = 1:length(splchans)
      chans(n) = str2num(splchans{n});
    end
  else
    % if using channel labels
    for n = 1:length(splchans)
      try
        chans(n) = find(ismember(info.data.label, splchans(n)));
      catch
        update_log(info.output_box, sprintf('ERROR: Please ensure the channel name is correct (case-sensitive)!'));
        uiresume;
        return;
      end
    end
  end
else
  update_log(info.output_box, sprintf('Please enter one or more channels'));
  uiresume;
  return;
end
try
  toggle = chans;
  info.chansel(toggle) = ~info.chansel(toggle);
catch
  update_log(info.output_box, sprintf('ERROR: Channel value too large!'));
  uiresume;
  return
end
% recalculate the metric
compute_metric(h)
guidata(h, info);
uiresume;

function toggle_visual(h, eventdata)
% copied from select2d, without waitforbuttonpress command
point1 = get(gca, 'CurrentPoint');    % button down detected
finalRect = rbbox;                    % return figure units
point2 = get(gca, 'CurrentPoint');    % button up detected
point1 = point1(1, 1:2);              % extract x and y
point2 = point2(1, 1:2);
x = sort([point1(1) point2(1)]);
y = sort([point1(2) point2(2)]);

g     = get(gca, 'Parent');
info  = guidata(g);

[maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = set_maxper(info.level, info.chansel, info.trlsel, strcmp(info.metric, 'min'));

switch gca
  case info.axes(1)
    % visual selection in the summary plot is not supported
    
  case info.axes(2)
    % the visual selection was made in the channels plot
    switch info.cfg.viewmode
      case 'toggle'
        chanlabels = 1:info.nchan;
        toggle = ...
          chanlabels >= y(1) & ...
          chanlabels <= y(2) & ...
          maxperchan_all(chanlabels)' >= x(1) & ...
          maxperchan_all(chanlabels)' <= x(2);
        info.chansel(toggle) = ~info.chansel(toggle);
        
      case 'remove'
        chanlabels     = 1:info.nchan;
        toggle = ...
          chanlabels >= y(1) & ...
          chanlabels <= y(2) & ...
          maxperchan(chanlabels)' >= x(1) & ...
          maxperchan(chanlabels)' <= x(2);
        info.chansel(toggle) = false;
        
      case 'hide'
        chanlabels = 1:sum(info.chansel==1);
        [junk, origchanlabels] = find(info.chansel==1);
        toggle = ...
          chanlabels >= y(1) & ...
          chanlabels <= y(2) & ...
          maxperchan(origchanlabels)' >= x(1) & ...
          maxperchan(origchanlabels)' <= x(2);
        info.chansel(origchanlabels(toggle)) = false;
    end
    
  case info.axes(3)
    % the visual selection was made in the trials plot
    switch info.cfg.viewmode
      case 'toggle'
        trllabels = 1:info.ntrl;
        toggle = ...
          trllabels >= x(1) & ...
          trllabels <= x(2) & ...
          maxpertrl_all(trllabels) >= y(1) & ...
          maxpertrl_all(trllabels) <= y(2);
        info.trlsel(toggle) = ~info.trlsel(toggle);
        
      case 'remove'
        trllabels = 1:info.ntrl;
        toggle = ...
          trllabels >= x(1) & ...
          trllabels <= x(2) & ...
          maxpertrl(trllabels) >= y(1) & ...
          maxpertrl(trllabels) <= y(2);
        info.trlsel(toggle) = false;
        
      case 'hide'
        trllabels = 1:sum(info.trlsel==1);
        [junk, origtrllabels] = find(info.trlsel==1);
        toggle = ...
          trllabels >= x(1) & ...
          trllabels <= x(2) & ...
          maxpertrl(origtrllabels) >= y(1) & ...
          maxpertrl(origtrllabels) <= y(2);
        info.trlsel(origtrllabels(toggle)) = false;
    end
    
    
end % switch gca

% recalculate the metric
compute_metric(h);
guidata(h, info);
uiresume;

% function display_trial(h, eventdata)
% info = guidata(h);
% rawtrls = get(h, 'string');
% if ~isempty(rawtrls)
%   spltrls = regexp(rawtrls, ' ', 'split');
%   trls = [];
%   for n = 1:length(spltrls)
%     trls(n) = str2num(cell2mat(spltrls(n)));
%   end
% else
%   update_log(info.output_box, sprintf('Please enter one or more trials'));
%   uiresume;
%   return;
% end
% if all(trls==0)
%   % use visual selection
%   update_log(info.output_box, sprintf('make visual selection of trials to be plotted seperately...'));
%   [x, y] = select2d;
%   maxpertrl  = max(info.origlevel, [], 1);
%   toggle = find(1:ntrl>=x(1) & ...
%     1:ntrl<=x(2) & ...
%     maxpertrl(:)'>=y(1) & ...
%     maxpertrl(:)'<=y(2));
% else
%   toggle = trls;
% end
% for i=1:length(trls)
%   figure
%   % the data being displayed here is NOT filtered
%   %plot(data.time{toggle(i)}, data.trial{toggle(i)}(chansel, :));
%   tmp = info.data.trial{toggle(i)}(info.chansel, :);
%   tmp = tmp - repmat(mean(tmp, 2), [1 size(tmp, 2)]);
%   plot(info.data.time{toggle(i)}, tmp);
%   title(sprintf('trial %d', toggle(i)));
% end

function quit(h, eventdata)
info = guidata(h);
info.quit = 1;
guidata(h, info);
uiresume;

function change_metric(h, eventdata)
info = guidata(h);
info.metric = get(eventdata.NewValue, 'string');
guidata(h, info);
compute_metric(h);
uiresume;

function toggle_rejected(h, eventdata)
info = guidata(h);
toggle = get(h, 'value');
if toggle == 0
  info.cfg.viewmode = 'remove';
else
  info.cfg.viewmode = 'toggle';
end
guidata(h, info);
uiresume;

function update_log(h, new_text)
new_text        = [datestr(now, 13) '# ' new_text];
curr_text       = get(h, 'string');
size_curr_text  = size(curr_text, 2);
size_new_text   = size(new_text, 2);
if size_curr_text > size_new_text
  new_text        = [new_text blanks(size_curr_text-size_new_text)];
else
  curr_text   = [curr_text repmat(blanks(size_new_text-size_curr_text), size(curr_text, 1), 1)];
end
set(h, 'String', [new_text; curr_text]);
drawnow;

function [maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = set_maxper(level, chansel, trlsel, minflag)
if minflag
  % take the negative maximum, i.e. the minimum
  level = -1 * level;
end
% determine the maximum value
maxperchan_all = max(level, [], 2);
maxpertrl_all  = max(level, [], 1);
% determine the maximum value over the remaining selection
level(~chansel, :) = nan;
level(:, ~trlsel)  = nan;
maxperchan     = max(level, [], 2);
maxpertrl      = max(level, [], 1);
if minflag
  maxperchan     = -1 * maxperchan;
  maxpertrl      = -1 * maxpertrl;
  maxperchan_all = -1 * maxperchan_all;
  maxpertrl_all  = -1 * maxpertrl_all;
  level          = -1 * level;
end

function display_trial(h, eventdata)
info = guidata(h);
update_log(info.output_box, 'Making multiplot of individual trials ...');
rawtrls = get(h, 'string');
if isempty(rawtrls)
  return;
else
  spltrls = regexp(rawtrls, '\s+', 'split');
  trls = [];
  for n = 1:length(spltrls)
    trls(n) = str2num(cell2mat(spltrls(n)));
  end
end
cfg_mp = [];
% disable hashing of input data (speeds up things)
cfg_mp.trackcallinfo = 'no';
cfg_mp.layout  = info.cfg.layout;
cfg_mp.channel = info.data.label(info.chansel);
currfig = gcf;
for n = 1:length(trls)
  % ft_multiplotER should be able to make the selection, but fails due to http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2978
  % that bug is hard to fix, hence it is solved here with a work-around
  cfg_sd = [];
  cfg_sd.trials = trls(n);
  cfg_sd.avgoverrpt = 'yes';
  cfg_sd.keeprpt = 'no';
  tmpdata = ft_selectdata(cfg_sd, info.data);
  
  figure()
  cfg_mp.interactive = 'yes';
  ft_multiplotER(cfg_mp, tmpdata);
  title(sprintf('Trial %i', trls(n)));
end
figure(currfig);
update_log(info.output_box, 'Done.');
return;

function offset = time2offset(time, fsample)

offset = round(time(1)*fsample);

function time = offset2time(offset, fsample, nsamples)

% ensure that these are not integers
offset   = double(offset);
nsamples = double(nsamples);

time = (offset + (0:(nsamples-1)))/fsample;

function [dat, label, time, cfg] = preproc(dat, label, time, cfg, begpadding, endpadding)

% PREPROC applies various preprocessing steps on a piece of EEG/MEG data
% that already has been read from a data file.
%
% This function can serve as a subfunction for all FieldTrip modules that
% want to preprocess the data, such as PREPROCESSING, ARTIFACT_XXX,
% TIMELOCKANALYSIS, etc. It ensures consistent handling of both MEG and EEG
% data and consistency in the use of all preprocessing configuration
% options.
%
% Use as
%   [dat, label, time, cfg] = preproc(dat, label, time, cfg, begpadding, endpadding)
%
% The required input arguments are
%   dat         Nchan x Ntime data matrix
%   label       Nchan x 1 cell-array with channel labels
%   time        Ntime x 1 vector with the latency in seconds
%   cfg         configuration structure, see below
% and the optional input arguments are
%   begpadding  number of samples that was used for padding (see below)
%   endpadding  number of samples that was used for padding (see below)
%
% The output is
%   dat         Nchan x Ntime data matrix
%   label       Nchan x 1 cell-array with channel labels
%   time        Ntime x 1 vector with the latency in seconds
%   cfg         configuration structure, optionally with extra defaults set
%
% Note that the number of input channels and the number of output channels
% can be different, for example when the user specifies that he/she wants
% to add the implicit EEG reference channel to the data matrix.
%
% The filtering of the data can introduce artifacts at the edges, hence it
% is better to pad the data with some extra signal at the begin and end.
% After filtering, this padding is removed and the other preprocessing
% steps are applied to the remainder of the data. The input fields
% begpadding and endpadding should be specified in samples. You can also
% leave them empty, which implies that the data is not padded.
%
% The configuration can contain
%   cfg.lpfilter      = 'no' or 'yes'  lowpass filter
%   cfg.hpfilter      = 'no' or 'yes'  highpass filter
%   cfg.bpfilter      = 'no' or 'yes'  bandpass filter
%   cfg.bsfilter      = 'no' or 'yes'  bandstop filter
%   cfg.dftfilter     = 'no' or 'yes'  line noise removal using discrete fourier transform
%   cfg.medianfilter  = 'no' or 'yes'  jump preserving median filter
%   cfg.lpfreq        = lowpass  frequency in Hz
%   cfg.hpfreq        = highpass frequency in Hz
%   cfg.bpfreq        = bandpass frequency range, specified as [low high] in Hz
%   cfg.bsfreq        = bandstop frequency range, specified as [low high] in Hz
%   cfg.dftfreq       = line noise frequencies for DFT filter, default [50 100 150] Hz
%   cfg.lpfiltord     = lowpass  filter order (default set in low-level function)
%   cfg.hpfiltord     = highpass filter order (default set in low-level function)
%   cfg.bpfiltord     = bandpass filter order (default set in low-level function)
%   cfg.bsfiltord     = bandstop filter order (default set in low-level function)
%   cfg.medianfiltord = length of median filter
%   cfg.lpfilttype    = digital filter type, 'but' (default) or 'firws' or 'fir' or 'firls'
%   cfg.hpfilttype    = digital filter type, 'but' (default) or 'firws' or 'fir' or 'firls'
%   cfg.bpfilttype    = digital filter type, 'but' (default) or 'firws' or 'fir' or 'firls'
%   cfg.bsfilttype    = digital filter type, 'but' (default) or 'firws' or 'fir' or 'firls'
%   cfg.lpfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse' or 'onepass-zerophase' (default for firws) or 'onepass-minphase' (firws, non-linear!)
%   cfg.hpfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse' or 'onepass-zerophase' (default for firws) or 'onepass-minphase' (firws, non-linear!)
%   cfg.bpfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse' or 'onepass-zerophase' (default for firws) or 'onepass-minphase' (firws, non-linear!)
%   cfg.bsfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse' or 'onepass-zerophase' (default for firws) or 'onepass-minphase' (firws, non-linear!)
%   cfg.lpinstabilityfix = deal with filter instability, 'no', 'reduce', 'split' (default  = 'no')
%   cfg.hpinstabilityfix = deal with filter instability, 'no', 'reduce', 'split' (default  = 'no')
%   cfg.bpinstabilityfix = deal with filter instability, 'no', 'reduce', 'split' (default  = 'no')
%   cfg.bsinstabilityfix = deal with filter instability, 'no', 'reduce', 'split' (default  = 'no')
%   cfg.lpfiltdf      = lowpass transition width (firws, overrides order, default set in low-level function)
%   cfg.hpfiltdf      = highpass transition width (firws, overrides order, default set in low-level function)
%   cfg.bpfiltdf      = bandpass transition width (firws, overrides order, default set in low-level function)
%   cfg.bsfiltdf      = bandstop transition width (firws, overrides order, default set in low-level function)
%   cfg.lpfiltwintype = lowpass window type, 'hann' or 'hamming' (default) or 'blackman' or 'kaiser' (firws)
%   cfg.hpfiltwintype = highpass window type, 'hann' or 'hamming' (default) or 'blackman' or 'kaiser' (firws)
%   cfg.bpfiltwintype = bandpass window type, 'hann' or 'hamming' (default) or 'blackman' or 'kaiser' (firws)
%   cfg.bsfiltwintype = bandstop window type, 'hann' or 'hamming' (default) or 'blackman' or 'kaiser' (firws)
%   cfg.lpfiltdev     = lowpass max passband deviation (firws with 'kaiser' window, default 0.001 set in low-level function)
%   cfg.hpfiltdev     = highpass max passband deviation (firws with 'kaiser' window, default 0.001 set in low-level function)
%   cfg.bpfiltdev     = bandpass max passband deviation (firws with 'kaiser' window, default 0.001 set in low-level function)
%   cfg.bsfiltdev     = bandstop max passband deviation (firws with 'kaiser' window, default 0.001 set in low-level function)
%   cfg.dftreplace    = 'zero' or 'neighbour', method used to reduce line noise, 'zero' implies DFT filter, 'neighbour' implies spectrum interpolation (default = 'zero')
%   cfg.dftbandwidth  = bandwidth of line noise frequencies, applies to spectrum interpolation, in Hz (default = [1 2 3])
%   cfg.dftneighbourwidth = bandwidth of frequencies neighbouring line noise frequencies, applies to spectrum interpolation, in Hz (default = [2 2 2])
%   cfg.plotfiltresp  = 'no' or 'yes', plot filter responses (firws, default = 'no')
%   cfg.usefftfilt    = 'no' or 'yes', use fftfilt instead of filter (firws, default = 'no')
%   cfg.demean        = 'no' or 'yes'
%   cfg.baselinewindow = [begin end] in seconds, the default is the complete trial
%   cfg.detrend       = 'no' or 'yes', this is done on the complete trial
%   cfg.polyremoval   = 'no' or 'yes', this is done on the complete trial
%   cfg.polyorder     = polynome order (default = 2)
%   cfg.derivative    = 'no' (default) or 'yes', computes the first order derivative of the data
%   cfg.hilbert       = 'no', 'abs', 'complex', 'real', 'imag', 'absreal', 'absimag' or 'angle' (default = 'no')
%   cfg.rectify       = 'no' or 'yes'
%   cfg.precision     = 'single' or 'double' (default = 'double')
%   cfg.absdiff       = 'no' or 'yes', computes absolute derivative (i.e.first derivative then rectify)
%
% Preprocessing options that you should only use for EEG data are
%   cfg.reref         = 'no' or 'yes' (default = 'no')
%   cfg.refchannel    = cell-array with new EEG reference channel(s)
%   cfg.refmethod     = 'avg', 'median', or 'bipolar' (default = 'avg')
%   cfg.implicitref   = 'label' or empty, add the implicit EEG reference as zeros (default = [])
%   cfg.montage       = 'no' or a montage structure (default = 'no')
%
% See also FT_READ_DATA, FT_READ_HEADER

% TODO implement decimation and/or resampling

% Copyright (C) 2004-2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% compute fsample
fsample = 1./nanmean(diff(time));

if nargin<5 || isempty(begpadding)
  begpadding = 0;
end
if nargin<6 || isempty(endpadding)
  endpadding = 0;
end

if iscell(cfg)
  % recurse over the subsequent preprocessing stages
  if begpadding>0 || endpadding>0
    ft_error('multiple preprocessing stages are not supported in combination with filter padding');
  end
  for i=1:length(cfg)
    tmpcfg = cfg{i};
    if nargout==1
      [dat                     ] = preproc(dat, label, time, tmpcfg, begpadding, endpadding);
    elseif nargout==2
      [dat, label              ] = preproc(dat, label, time, tmpcfg, begpadding, endpadding);
    elseif nargout==3
      [dat, label, time        ] = preproc(dat, label, time, tmpcfg, begpadding, endpadding);
    elseif nargout==4
      [dat, label, time, tmpcfg] = preproc(dat, label, time, tmpcfg, begpadding, endpadding);
      cfg{i} = tmpcfg;
    end
  end
  % ready with recursing over the subsequent preprocessing stages
  return
end

% set the defaults for the rereferencing options
cfg.reref =                ft_getopt(cfg, 'reref', 'no');
cfg.refchannel =           ft_getopt(cfg, 'refchannel', {});
cfg.refmethod =            ft_getopt(cfg, 'refmethod', 'avg');
cfg.implicitref =          ft_getopt(cfg, 'implicitref', []);
% set the defaults for the signal processing options
cfg.polyremoval =          ft_getopt(cfg, 'polyremoval', 'no');
cfg.polyorder =            ft_getopt(cfg, 'polyorder', 2);
cfg.detrend =              ft_getopt(cfg, 'detrend', 'no');
cfg.demean =               ft_getopt(cfg, 'demean', 'no');
cfg.baselinewindow =       ft_getopt(cfg, 'baselinewindow', 'all');
cfg.dftfilter =            ft_getopt(cfg, 'dftfilter', 'no');
cfg.lpfilter =             ft_getopt(cfg, 'lpfilter', 'no');
cfg.hpfilter =             ft_getopt(cfg, 'hpfilter', 'no');
cfg.bpfilter =             ft_getopt(cfg, 'bpfilter', 'no');
cfg.bsfilter =             ft_getopt(cfg, 'bsfilter', 'no');
cfg.lpfiltord =            ft_getopt(cfg, 'lpfiltord', []);
cfg.hpfiltord =            ft_getopt(cfg, 'hpfiltord', []);
cfg.bpfiltord =            ft_getopt(cfg, 'bpfiltord', []);
cfg.bsfiltord =            ft_getopt(cfg, 'bsfiltord', []);
cfg.lpfilttype =           ft_getopt(cfg, 'lpfilttype', 'but');
cfg.hpfilttype =           ft_getopt(cfg, 'hpfilttype', 'but');
cfg.bpfilttype =           ft_getopt(cfg, 'bpfilttype', 'but');
cfg.bsfilttype =           ft_getopt(cfg, 'bsfilttype', 'but');
if strcmp(cfg.lpfilttype, 'firws'), cfg.lpfiltdir = ft_getopt(cfg, 'lpfiltdir', 'onepass-zerophase'); else, cfg.lpfiltdir = ft_getopt(cfg, 'lpfiltdir', 'twopass'); end
if strcmp(cfg.hpfilttype, 'firws'), cfg.hpfiltdir = ft_getopt(cfg, 'hpfiltdir', 'onepass-zerophase'); else, cfg.hpfiltdir = ft_getopt(cfg, 'hpfiltdir', 'twopass'); end
if strcmp(cfg.bpfilttype, 'firws'), cfg.bpfiltdir = ft_getopt(cfg, 'bpfiltdir', 'onepass-zerophase'); else, cfg.bpfiltdir = ft_getopt(cfg, 'bpfiltdir', 'twopass'); end
if strcmp(cfg.bsfilttype, 'firws'), cfg.bsfiltdir = ft_getopt(cfg, 'bsfiltdir', 'onepass-zerophase'); else, cfg.bsfiltdir = ft_getopt(cfg, 'bsfiltdir', 'twopass'); end
cfg.lpinstabilityfix =     ft_getopt(cfg, 'lpinstabilityfix', 'no');
cfg.hpinstabilityfix =     ft_getopt(cfg, 'hpinstabilityfix', 'no');
cfg.bpinstabilityfix =     ft_getopt(cfg, 'bpinstabilityfix', 'no');
cfg.bsinstabilityfix =     ft_getopt(cfg, 'bsinstabilityfix', 'no');
cfg.lpfiltdf =             ft_getopt(cfg, 'lpfiltdf',  []);
cfg.hpfiltdf =             ft_getopt(cfg, 'hpfiltdf', []);
cfg.bpfiltdf =             ft_getopt(cfg, 'bpfiltdf', []);
cfg.bsfiltdf =             ft_getopt(cfg, 'bsfiltdf', []);
cfg.lpfiltwintype =        ft_getopt(cfg, 'lpfiltwintype', 'hamming');
cfg.hpfiltwintype =        ft_getopt(cfg, 'hpfiltwintype', 'hamming');
cfg.bpfiltwintype =        ft_getopt(cfg, 'bpfiltwintype', 'hamming');
cfg.bsfiltwintype =        ft_getopt(cfg, 'bsfiltwintype', 'hamming');
cfg.lpfiltdev =            ft_getopt(cfg, 'lpfiltdev', []);
cfg.hpfiltdev =            ft_getopt(cfg, 'hpfiltdev', []);
cfg.bpfiltdev =            ft_getopt(cfg, 'bpfiltdev', []);
cfg.bsfiltdev =            ft_getopt(cfg, 'bsfiltdev', []);
cfg.plotfiltresp =         ft_getopt(cfg, 'plotfiltresp', 'no');
cfg.usefftfilt =           ft_getopt(cfg, 'usefftfilt', 'no');
cfg.medianfilter =         ft_getopt(cfg, 'medianfilter ', 'no');
cfg.medianfiltord =        ft_getopt(cfg, 'medianfiltord', 9);
cfg.dftfreq =              ft_getopt(cfg, 'dftfreq', [50 100 150]);
cfg.hilbert =              ft_getopt(cfg, 'hilbert', 'no');
cfg.derivative =           ft_getopt(cfg, 'derivative', 'no');
cfg.rectify =              ft_getopt(cfg, 'rectify', 'no');
cfg.boxcar =               ft_getopt(cfg, 'boxcar', 'no');
cfg.absdiff =              ft_getopt(cfg, 'absdiff', 'no');
cfg.precision =            ft_getopt(cfg, 'precision', []);
cfg.conv =                 ft_getopt(cfg, 'conv', 'no');
cfg.montage =              ft_getopt(cfg, 'montage', 'no');
cfg.dftinvert =            ft_getopt(cfg, 'dftinvert', 'no');
cfg.standardize =          ft_getopt(cfg, 'standardize', 'no');
cfg.denoise =              ft_getopt(cfg, 'denoise', '');
cfg.subspace =             ft_getopt(cfg, 'subspace', []);
cfg.custom =               ft_getopt(cfg, 'custom', '');
cfg.resample =             ft_getopt(cfg, 'resample', '');

% test whether the MATLAB signal processing toolbox is available
if strcmp(cfg.medianfilter, 'yes') && ~ft_hastoolbox('signal')
  ft_error('median filtering requires the MATLAB signal processing toolbox');
end

% do a sanity check on the filter configuration
if strcmp(cfg.bpfilter, 'yes') && ...
    (strcmp(cfg.hpfilter, 'yes') || strcmp(cfg.lpfilter,'yes'))
  ft_error('you should not apply both a bandpass AND a lowpass/highpass filter');
end

% do a sanity check on the hilbert transform configuration
if strcmp(cfg.hilbert, 'yes') && ~strcmp(cfg.bpfilter, 'yes')
  ft_warning('Hilbert transform should be applied in conjunction with bandpass filter')
end

% do a sanity check on hilbert and rectification
if strcmp(cfg.hilbert, 'yes') && strcmp(cfg.rectify, 'yes')
  ft_error('Hilbert transform and rectification should not be applied both')
end

% do a sanity check on the rereferencing/montage
if ~strcmp(cfg.reref, 'no') && ~strcmp(cfg.montage, 'no')
  ft_error('cfg.reref and cfg.montage are mutually exclusive')
end

% lnfilter is no longer used
if isfield(cfg, 'lnfilter') && strcmp(cfg.lnfilter, 'yes')
  ft_error('line noise filtering using the option cfg.lnfilter is not supported any more, use cfg.bsfilter instead')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do the rereferencing in case of EEG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(cfg.implicitref) && ~any(match_str(cfg.implicitref,label))
  label = {label{:} cfg.implicitref};
  dat(end+1,:) = 0;
end

if strcmp(cfg.reref, 'yes')
  if strcmp(cfg.refmethod, 'bipolar')
    % this is implemented as a montage that the user does not get to see
    tmpcfg = keepfields(cfg, {'refmethod', 'implicitref', 'refchannel', 'channel'});
    tmpcfg.showcallinfo = 'no';
    montage = ft_prepare_montage(tmpcfg);
    % convert the data temporarily to a raw structure
    tmpdata.trial = {dat};
    tmpdata.time  = {time};
    tmpdata.label = label;
    % apply the montage to the data
    tmpdata = ft_apply_montage(tmpdata, montage, 'feedback', 'none');
    dat   = tmpdata.trial{1}; % the number of channels can have changed
    label = tmpdata.label;    % the output channels can be different than the input channels
    clear tmpdata
  else
    % mean or median based derivation of specified or all channels
    cfg.refchannel = ft_channelselection(cfg.refchannel, label);
    refindx = match_str(label, cfg.refchannel);
    if isempty(refindx)
      ft_error('reference channel was not found')
    end
    dat = ft_preproc_rereference(dat, refindx, cfg.refmethod);
  end
end

if ~strcmp(cfg.montage, 'no') && ~isempty(cfg.montage)
  % convert the data temporarily to a raw structure
  tmpdata.trial = {dat};
  tmpdata.time  = {time};
  tmpdata.label = label;
  % apply the montage to the data
  tmpdata = ft_apply_montage(tmpdata, cfg.montage, 'feedback', 'none');
  dat   = tmpdata.trial{1}; % the number of channels can have changed
  label = tmpdata.label;    % the output channels can be different than the input channels
  clear tmpdata
end

if any(any(isnan(dat)))
  % filtering is not possible for at least a selection of the data
  ft_warning('data contains NaNs, no filtering or preprocessing applied');
  
else
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % do the filtering on the padded data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ~isempty(cfg.denoise)
    hflag    = isfield(cfg.denoise, 'hilbert') && strcmp(cfg.denoise.hilbert, 'yes');
    datlabel = match_str(label, cfg.denoise.channel);
    reflabel = match_str(label, cfg.denoise.refchannel);
    tmpdat   = ft_preproc_denoise(dat(datlabel,:), dat(reflabel,:), hflag);
    dat(datlabel,:) = tmpdat;
  end
  
  % The filtering should in principle be done prior to the demeaning to
  % ensure that the resulting mean over the baseline window will be
  % guaranteed to be zero (even if there are filter artifacts).
  % However, the filtering benefits from the data being pulled towards zero,
  % causing less edge artifacts. That is why we start by removing the slow
  % drift, then filter, and then repeat the demean/detrend/polyremove.
  if strcmp(cfg.polyremoval, 'yes')
    nsamples  = size(dat,2);
    begsample = 1        + begpadding;
    endsample = nsamples - endpadding;
    dat = ft_preproc_polyremoval(dat, cfg.polyorder, begsample, endsample); % this will also demean and detrend
  elseif strcmp(cfg.detrend, 'yes')
    nsamples  = size(dat,2);
    begsample = 1        + begpadding;
    endsample = nsamples - endpadding;
    dat = ft_preproc_polyremoval(dat, 1, begsample, endsample); % this will also demean
  elseif strcmp(cfg.demean, 'yes')
    nsamples  = size(dat,2);
    begsample = 1        + begpadding;
    endsample = nsamples - endpadding;
    dat = ft_preproc_polyremoval(dat, 0, begsample, endsample);
  end
  
  if strcmp(cfg.medianfilter, 'yes'), dat = ft_preproc_medianfilter(dat, cfg.medianfiltord); end
  if strcmp(cfg.lpfilter, 'yes'),     dat = ft_preproc_lowpassfilter(dat, fsample, cfg.lpfreq, cfg.lpfiltord, cfg.lpfilttype, cfg.lpfiltdir, cfg.lpinstabilityfix, cfg.lpfiltdf, cfg.lpfiltwintype, cfg.lpfiltdev, cfg.plotfiltresp, cfg.usefftfilt); end
  if strcmp(cfg.hpfilter, 'yes'),     dat = ft_preproc_highpassfilter(dat, fsample, cfg.hpfreq, cfg.hpfiltord, cfg.hpfilttype, cfg.hpfiltdir, cfg.hpinstabilityfix, cfg.hpfiltdf, cfg.hpfiltwintype, cfg.hpfiltdev, cfg.plotfiltresp, cfg.usefftfilt); end
  if strcmp(cfg.bpfilter, 'yes'),     dat = ft_preproc_bandpassfilter(dat, fsample, cfg.bpfreq, cfg.bpfiltord, cfg.bpfilttype, cfg.bpfiltdir, cfg.bpinstabilityfix, cfg.bpfiltdf, cfg.bpfiltwintype, cfg.bpfiltdev, cfg.plotfiltresp, cfg.usefftfilt); end
  if strcmp(cfg.bsfilter, 'yes')
    for i=1:size(cfg.bsfreq,1)
      % apply a bandstop filter for each of the specified bands, i.e. cfg.bsfreq should be Nx2
      dat = ft_preproc_bandstopfilter(dat, fsample, cfg.bsfreq(i,:), cfg.bsfiltord, cfg.bsfilttype, cfg.bsfiltdir, cfg.bsinstabilityfix, cfg.bsfiltdf, cfg.bsfiltwintype, cfg.bsfiltdev, cfg.plotfiltresp, cfg.usefftfilt);
    end
  end
  if strcmp(cfg.polyremoval, 'yes')
    % the begin and endsample of the polyremoval period correspond to the complete data minus padding
    nsamples  = size(dat,2);
    begsample = 1        + begpadding;
    endsample = nsamples - endpadding;
    dat = ft_preproc_polyremoval(dat, cfg.polyorder, begsample, endsample);
  end
  if strcmp(cfg.detrend, 'yes')
    % the begin and endsample of the detrend period correspond to the complete data minus padding
    nsamples  = size(dat,2);
    begsample = 1        + begpadding;
    endsample = nsamples - endpadding;
    dat = ft_preproc_detrend(dat, begsample, endsample);
  end
  if strcmp(cfg.demean, 'yes')
    if ischar(cfg.baselinewindow) && strcmp(cfg.baselinewindow, 'all')
      % the begin and endsample of the baseline period correspond to the complete data minus padding
      nsamples  = size(dat,2);
      begsample = 1        + begpadding;
      endsample = nsamples - endpadding;
      dat       = ft_preproc_baselinecorrect(dat, begsample, endsample);
    else
      % determine the begin and endsample of the baseline period and baseline correct for it
      begsample = nearest(time, cfg.baselinewindow(1));
      endsample = nearest(time, cfg.baselinewindow(2));
      dat       = ft_preproc_baselinecorrect(dat, begsample, endsample);
    end
  end
  if strcmp(cfg.dftfilter, 'yes')
    datorig = dat;
    optarg = {};
    if isfield(cfg, 'dftreplace')
      optarg = cat(2, optarg, {'dftreplace', cfg.dftreplace});
      if strcmp(cfg.dftreplace, 'neighbour') && (begpadding>0 || endpadding>0)
        ft_error('Padding by data mirroring is not supported for spectrum interpolation.');
      end
    end
    if isfield(cfg, 'dftbandwidth')
      optarg = cat(2, optarg, {'dftbandwidth', cfg.dftbandwidth});
    end
    if isfield(cfg, 'dftneighbourwidth')
      optarg = cat(2, optarg, {'dftneighbourwidth', cfg.dftneighbourwidth});
    end
    dat     = ft_preproc_dftfilter(dat, fsample, cfg.dftfreq, optarg{:});
    if strcmp(cfg.dftinvert, 'yes')
      dat = datorig - dat;
    end
  end
  if ~strcmp(cfg.hilbert, 'no')
    dat = ft_preproc_hilbert(dat, cfg.hilbert);
  end
  if strcmp(cfg.rectify, 'yes')
    dat = ft_preproc_rectify(dat);
  end
  if isnumeric(cfg.boxcar)
    numsmp = round(cfg.boxcar*fsample);
    if ~rem(numsmp,2)
      % the kernel should have an odd number of samples
      numsmp = numsmp+1;
    end
    % kernel = ones(1,numsmp) ./ numsmp;
    % dat    = convn(dat, kernel, 'same');
    dat = ft_preproc_smooth(dat, numsmp); % better edge behaviour
  end
  if isnumeric(cfg.conv)
    kernel = (cfg.conv(:)'./sum(cfg.conv));
    if ~rem(length(kernel),2)
      kernel = [kernel 0];
    end
    dat = convn(dat, kernel, 'same');
  end
  if strcmp(cfg.derivative, 'yes')
    dat = ft_preproc_derivative(dat, 1);
  end
  if strcmp(cfg.absdiff, 'yes')
    % this implements abs(diff(data), which is required for jump detection
    dat = abs([diff(dat, 1, 2) zeros(size(dat,1),1)]);
  end
  if strcmp(cfg.standardize, 'yes')
    dat = ft_preproc_standardize(dat, 1, size(dat,2));
  end
  if ~isempty(cfg.subspace)
    dat = ft_preproc_subspace(dat, cfg.subspace);
  end
  if ~isempty(cfg.custom)
    if ~isfield(cfg.custom, 'nargout')
      cfg.custom.nargout = 1;
    end
    if cfg.custom.nargout==1
      dat = feval(cfg.custom.funhandle, dat, cfg.custom.varargin);
    elseif cfg.custom.nargout==2
      [dat, time] = feval(cfg.custom.funhandle, dat, cfg.custom.varargin);
    end
  end
  if strcmp(cfg.resample, 'yes')
    if ~isfield(cfg, 'resamplefs')
      cfg.resamplefs = fsample./2;
    end
    if ~isfield(cfg, 'resamplemethod')
      cfg.resamplemethod = 'resample';
    end
    [dat               ] = ft_preproc_resample(dat,  fsample, cfg.resamplefs, cfg.resamplemethod);
    [time, dum, fsample] = ft_preproc_resample(time, fsample, cfg.resamplefs, cfg.resamplemethod);
  end
  if ~isempty(cfg.precision)
    % convert the data to another numeric precision, i.e. double, single or int32
    dat = cast(dat, cfg.precision);
  end
end % if any(isnan)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove the filter padding and do the preprocessing on the remaining trial data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if begpadding~=0 || endpadding~=0
  dat = ft_preproc_padding(dat, 'remove', begpadding, endpadding);
  if strcmp(cfg.demean, 'yes') || nargout>2
    time = ft_preproc_padding(time, 'remove', begpadding, endpadding);
  end
end

function [cfg, varargout] = rollback_provenance(cfg, varargin)

% ROLLBACK_PROVENANCE rolls the provenance one step back and should
% be used whenever a FT function calls another FT function without
% the user being (or having to be) aware of this.
%
% Some examples for use
%
%   tmpcfg            = [];
%   tmpcfg.downsample = cfg.downsample;  % simply copy this option
%   tmpcfg.smooth     = 'no';            % override the default for this option
%   mri = ft_volumedownsample(tmpcfg, mri);
%   [cfg, mri] = rollback_provenance(cfg, mri);
%
%   tmpcfg           = [];
%   tmpcfg.parameter = cfg.parameter;
%   [varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
%   [cfg, varargin{:}] = rollback_provenance(cfg, varargin{:});
%
% See also FT_PREAMBLE, FT_POSTAMBLE

% Copyright (C) 2013-2014, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

for i=1:(nargin-1)
  
  if ~isfield(varargin{i}, 'cfg')
    % nothing to do
    continue
  end
  
  if isempty(cfg)
    % allow for [] as the input cfg
    fn0 = {};
  else
    fn0 = fieldnames(cfg);
  end
  
  if ~isfield(varargin{i}, 'cfg')
    % input does not contain cfg, so no rollback to be performed
    continue;
  end
  
  fn1 = fieldnames(varargin{i}.cfg);
  
  % only work on the fields that are explicitly present in the cfg
  fn = intersect(fn0, fn1);
  
  % ignore the provenance fields themselves
  fn = setdiff(fn, { ...
    'trackconfig'
    'checkconfig'
    'checksize'
    'trackusage'
    'trackdatainfo'
    'trackcallinfo'
    'showcallinfo'
    'callinfo'
    'version'
    'warning'
    'debug'
    'previous'
    });
  
  for j=1:length(fn)
    cfg.(fn{j}) = varargin{i}.cfg.(fn{j});
  end % for all fields that overlap
  
  if isfield(varargin{i}.cfg, 'previous')
    % take it one step back
    varargin{i}.cfg = varargin{i}.cfg.previous;
  else
    % there is no previous information
    varargin{i} = rmfield(varargin{i}, 'cfg');
  end
  
end % for all input data structures

% return the updated data structures
varargout = varargin;
