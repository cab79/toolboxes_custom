function D = eegstats_dataprep(S)
% inputs: S is a settings structure (see example script)
% outputs: D is the prepared EEG data and predictors as a struct

%% preliminaries
D = struct;
if strcmp(S.prep.parallel.option,'condor')  
    delete('input*.mat');
    delete('output*.mat');
    delete('*.out');
    delete('*.log');
    delete('*.err');
end
% change to the input directory
eval(sprintf('%s', ['cd(''' S.prep.path.inputs.eeg ''')']));

%% find data using getfilelist.m function
GFL=S;
GFL.func = 'prep';
GFL.path.file = S.prep.path.inputs.eeg;
GFL.path.datfile = S.prep.path.inputs.subjects;
if isfield(S.prep.select,'prefixes') && ~isempty(S.prep.select.prefixes); GFL.(GFL.func).load.prefix = S.prep.select.prefixes; end
if isfield(S.prep.select,'suffixes') && ~isempty(S.prep.select.suffixes); GFL.(GFL.func).load.suffix = S.prep.select.suffixes; end
GFL = getfilelist(GFL);
if isfield(GFL.(GFL.func),'designtab')
    designmat=GFL.(GFL.func).designtab;
else
    designmat=cell2table(GFL.(GFL.func).designmat(2:end,:));
    designmat.Properties.VariableNames = GFL.(GFL.func).designmat(1,:);
end
filelist = designmat.file;
subjects = GFL.(GFL.func).select.subjects;
if isempty(filelist)
    error('No files found to import!\n');
end

% run though all EEG files in a loop
for d = 1:length(subjects)
    
    % create output D
    D(d).prep.subname = subjects(d);
        
    % FIND THE FILES FOR THIS SUBJECT
    subfiles = filelist(find(not(cellfun('isempty', strfind(filelist,D(d).prep.subname{:})))));
    
    if isempty(subfiles)
        continue
    end
    
    disp(['subject ' num2str(d) '/' num2str(length(subjects))])  

    % CYCLE THROUGH EACH FILE FOR THIS SUBJECT
    % add data to D structure which has a common format for any data type
    for f = 1:length(subfiles)
        filename = subfiles{f};
        loadsuffix = S.prep.select.suffixes{1};
        fprintf('\nImporting %s.\n\n', filename);

        switch S.prep.fname.ext{:}
            case 'mat'
                temp = load(filename);
            case 'set'
                temp = pop_loadset(filename);
        end
        fename = fieldnames(temp(f));

         % CORE study: apply filtering to simplify the data for PCA/CCA
        % later
        if ~isempty(S.prep.calc.eeg.filter)
            % Assuming you have already loaded your epoched EEG data into the variable 'EEG'

            % Set the filter parameters
            lowpass_cutoff = S.prep.calc.eeg.filter(2);  % Low-pass cutoff frequency (in Hz)
            highpass_cutoff = S.prep.calc.eeg.filter(1);  % High-pass cutoff frequency (in Hz)
            filter_order = 4;     % Filter order
            
            % Create a copy of the EEG data for buffering
            EEG=temp;
            EEG_filtered = temp;
            
            % Define the length of the buffer zone (in data points)
            buffer_length = S.prep.calc.eeg.filter_addbuffer;
            
            % Apply the buffer zone to each epoch
            for epoch = 1:length(EEG.epoch)
                % Create the buffer zone by flipping the first and last 100 samples
                buffer_zone = [fliplr(EEG.data(:, 1:buffer_length, epoch)), EEG.data(:, :, epoch), fliplr(EEG.data(:, end-buffer_length+1:end, epoch))];
                
                % Apply the filter to the buffer zone
                filtered_buffer_zone = eegfilt(buffer_zone, EEG.srate, highpass_cutoff, lowpass_cutoff);
                
                % Remove the buffer zone after filtering
                filtered_epoch = filtered_buffer_zone(:, (buffer_length+1):(end-buffer_length));
                
                % Update the filtered epoch in the EEG_filtered structure
                EEG_filtered.data(:, :, epoch) = filtered_epoch;
            end
            temp = EEG_filtered;
           
            
%             % Display the filtered data
%             pop_eegplot(EEG);
%             pop_eegplot(EEG_filtered);


        end

        if ~isempty(S.prep.calc.eeg.TF.band)

            % Specify channels to include in IAF calculation
            channels = find(ismember({temp.chanlocs.labels}, S.prep.calc.eeg.TF.iaf.channels));
            
            % Read IAF from file, or calculate IAF using specified channels
            if ~isempty(S.prep.calc.eeg.TF.iaf.file)
                loaded_dtab = readtable(S.prep.calc.eeg.TF.iaf.file);
                if S.prep.calc.eeg.TF.iaf.use_regularised
                    iaf = loaded_dtab.iaf_reg(find(ismember(loaded_dtab.ID,subjects{d})));
                else
                    iaf = loaded_dtab.iaf(find(ismember(loaded_dtab.ID,subjects{d})));
                end
            else
                iaf = calculateIAF(temp.data, temp.srate, channels, S.prep.calc.eeg.TF.iaf.plot);
            end
            
            % Define frequency bands relative to IAF
            freq_bands = {
                'Delta1', [0.5 2];
                'Delta2', [2 4];
                'Theta1', [4 6];
                'Theta2', [6 8];
                'Alpha1', [iaf - 2, iaf];
                'Alpha2', [iaf, iaf + 2];
                'Beta1', [iaf + 2, iaf + 5];
                'Beta2', [iaf + 5, iaf + 10];
                'Beta3', [iaf + 10, iaf + 20];
                'Gamma1', [iaf + 20, iaf + 30];
                'Gamma2', [iaf + 30, iaf + 50];
                'Gamma3', [iaf + 50, iaf + 70];
            };
            
            % Time window settings (adjust as needed)
            time_windows = {
                [2 0.5]; % Delta1
                [2 0.5]; % Delta2
                [0.5 0.05]; % Theta1
                [0.5 0.05]; % Theta2
                [0.5 0.05]; % Alpha1
                [0.5 0.05]; % Alpha2
                [0.5 0.05]; % Beta1
                [0.5 0.05]; % Beta2
                [0.5 0.05]; % Beta3
                [0.2 0.02]; % Gamma1
                [0.2 0.02]; % Gamma2
                [0.2 0.02]; % Gamma3
            };

            % Define taper settings for low and high frequency bands
            taper_settings = {
                [1 1]; % Delta1
                [1 1]; % Delta2
                [2 3]; % Theta1
                [2 3]; % Theta2
                [2 3]; % Alpha1
                [2 3]; % Alpha2
                [2 3]; % Beta1
                [2 3]; % Beta2
                [2 3]; % Beta3
                [4 7]; % Gamma1
                [4 7]; % Gamma2
                [4 7]; % Gamma3
            };
            
            % Select relevant frequency bands and time windows
            freq_bands_select = freq_bands(ismember(freq_bands(:,1), S.prep.calc.eeg.TF.band), :);
            time_windows_select = time_windows(ismember(freq_bands(:,1), S.prep.calc.eeg.TF.band), :);
            taper_settings_select = taper_settings(ismember(freq_bands(:,1), S.prep.calc.eeg.TF.band), :);

            
            % Example code to use these frequency bands for time-frequency analysis
            chronux_params = struct('Fs', temp.srate, 'tapers', [],'pad',S.prep.calc.eeg.TF.fft_pad);
            
            % Initialize matrix to hold time-frequency data for each trial
            num_channels = size(temp.data, 1); % Number of channels
            num_timepoints = size(temp.data, 2); % Number of time points
            num_trials = size(temp.data, 3); % Number of trials
            
            % Define padding length (e.g., 1 second of padding)
            padding_length = S.prep.calc.eeg.TF.padding * temp.srate; 
            
            % Initialize maximum number of time points from all spectrograms
            max_timepoints = 0;
            
            % Temporary storage for time-frequency data
            tf_data_temp = cell(num_channels, num_trials, length(freq_bands_select));
            time_vectors = cell(num_channels, num_trials, length(freq_bands_select));


            if S.prep.calc.eeg.TF.remove_ERP
                temp.data = temp.data - mean(temp.data, 3);
            end
            
            % Parallel pool setup
            if isempty(gcp('nocreate'))
                parpool; % Start parallel pool if not already started
            end
            
            % Loop through each channel and compute time-frequency representation
            parfor ch = 1:num_channels
                disp(['Channel ' num2str(ch)])
                local_tf_data_temp = cell(num_trials, length(freq_bands_select));
                local_time_vectors = cell(num_trials, length(freq_bands_select));
                local_max_timepoints = 0;
            
                for trial = 1:num_trials
                    trial_data = squeeze(temp.data(ch, :, trial)); % Data for single trial

                    
                    % Add zero padding
                    padded_data = zeros(1, num_timepoints + 2 * padding_length);
                    padded_data(padding_length + 1:padding_length + num_timepoints) = trial_data;
            
                    % Create a map to store results of mtspecgramc for each unique movingwin
                    mtspecgramc_results = containers.Map('KeyType', 'char', 'ValueType', 'any');
                    
                    % Identify unique movingwin parameters
                    unique_movingwin = unique(cellfun(@(x) mat2str([x(1), x(2)]), time_windows_select, 'UniformOutput', false));
                    
                    for uwin_idx = 1:length(unique_movingwin)
                        movingwin_str = unique_movingwin{uwin_idx};
                        movingwin = str2num(movingwin_str); %#ok<ST2NM> % Convert string back to numeric
          
                        % Find the corresponding taper settings
                        taper_setting_idx = find(cellfun(@(x) isequal([x(1), x(2)], movingwin), time_windows_select));
                        local_chronux_params = chronux_params;
                        local_chronux_params.tapers = taper_settings_select{taper_setting_idx(1)};
                        
                        % Compute spectrogram for this unique movingwin
                        [Sp, t, fr] = mtspecgramc(padded_data', movingwin, local_chronux_params);
                        
                        % Remove time points that were added due to padding
                        valid_indices = t >= (padding_length / temp.srate) & t <= ((padding_length + num_timepoints) / temp.srate);
                        t = t(valid_indices);
                        Sp = Sp(valid_indices, :);
                        
                        % Store the results in the map
                        mtspecgramc_results(movingwin_str) = {Sp, t, fr};
                    end
                    
                    for band_idx = 1:size(freq_bands_select, 1)
                        freq_range = freq_bands_select{band_idx, 2};
                        window_length = time_windows_select{band_idx}(1);
                        step_size = time_windows_select{band_idx}(2);
                        movingwin_str = mat2str([window_length, step_size]);
            
                        % Retrieve the precomputed spectrogram for this movingwin
                        result = mtspecgramc_results(movingwin_str);
                        Sp = result{1};
                        t = result{2};
                        fr = result{3};
                        
                        % Store temporary time-frequency data and time vectors
                        local_tf_data_temp{trial, band_idx} = Sp;
                        local_time_vectors{trial, band_idx} = t;
                        
                        % Update maximum number of time points
                        if length(t) > local_max_timepoints
                            local_max_timepoints = length(t);
                        end
                    end
                end
                
                % Update the shared maximum number of time points
                max_timepoints = max(max_timepoints, local_max_timepoints);
                
                % Store the local results in the global cell arrays
                tf_data_temp(ch, :, :) = local_tf_data_temp;
                time_vectors(ch, :, :) = local_time_vectors;
            end

            
            % Option to resample the number of time points to the maximum across bands
            resample_to_max = S.prep.calc.eeg.TF.resample_to_max;
            
            % Initialize final tf_data as a cell array
            tf_data = cell(size(freq_bands_select, 1), 1);
            
            for band_idx = 1:size(freq_bands_select, 1)
                if resample_to_max
                    % Initialize temporary storage for resampled data
                    tf_data_resampled = zeros(num_channels, max_timepoints, num_trials);
                    for trial = 1:num_trials
                        for ch = 1:num_channels
                            Sp = tf_data_temp{ch, trial, band_idx};
                            t = time_vectors{ch, trial, band_idx};
            
                            % Average over the frequency band of interest
                            freq_range = freq_bands_select{band_idx, 2};
                            fr = linspace(0, chronux_params.Fs / 2, size(Sp, 2));
                            [~, min_idx] = min(abs(fr - freq_range(1)));
                            [~, max_idx] = min(abs(fr - freq_range(2)));
                            freq_indices = min_idx:max_idx;
                            S_band = mean(Sp(:, freq_indices), 2); % Average over selected frequencies
            
                            % Resample S_band to match the maximum number of time points
                            S_band_resampled = interp1(t, S_band, linspace(t(1), t(end), max_timepoints));
                            time_vectors{ch, trial, band_idx} = interp1(t, t, linspace(t(1), t(end), max_timepoints));
            
                            % Store the result in tf_data_resampled
                            tf_data_resampled(ch, :, trial) = S_band_resampled; % Ensure correct dimensions
                        end
                    end
                    tf_data{band_idx} = tf_data_resampled;
                else
                    % Store the original data without resampling
                    tf_data_original = [];
                    for trial = 1:num_trials
                        for ch = 1:num_channels
                            Sp = tf_data_temp{ch, trial, band_idx};
                            % t = time_vectors{ch, trial, band_idx};
            
                            % Average over the frequency band of interest
                            freq_range = freq_bands_select{band_idx, 2};
                            fr = linspace(0, chronux_params.Fs / 2, size(Sp, 2));
                            [~, min_idx] = min(abs(fr - freq_range(1)));
                            [~, max_idx] = min(abs(fr - freq_range(2)));
                            freq_indices = min_idx:max_idx;
                            S_band = mean(Sp(:, freq_indices), 2); % Average over selected frequencies
            
                            % Store the result in tf_data_original
                            tf_data_original(ch, :, trial) = S_band; % Ensure correct dimensions
                        end
                    end
                    tf_data{band_idx} = tf_data_original;
                end
            end
            

            original_times = temp.times; % Time points in ms
            temp.times = {};
            temp.pnts = {};
            for band_idx = 1:size(freq_bands_select, 1)
                temp.times{band_idx} = (time_vectors{1, 1, band_idx}-(padding_length / temp.srate))*1000+original_times(1);
                temp.pnts{band_idx} = length(temp.times{band_idx});
            end
            
            if S.prep.calc.eeg.TF.plot_eachband
                for band_idx = 1:size(freq_bands_select, 1)
                    % Average the data over trials
                    avg_tf_data = mean(tf_data{band_idx}, 3);
                    
                    % Plot the time-frequency representation averaged over trials for the current frequency band
                    figure;
                    subplot(2, 1, 1);
                    imagesc(temp.times{band_idx}, 1:num_channels, 10*log10(avg_tf_data)); % Example
                    axis xy;
                    xlabel('Time (ms)');
                    ylabel('Channel');
                    title(['Multitaper Time-Frequency Averaged Over Trials - Band: ', freq_bands_select{band_idx, 1}]);
                    colorbar;
        
                    % Calculate the average power across selected IAF channels
                    avg_power = mean(10*log10(tf_data{band_idx}(channels, :, :)), 1); % Example
        
                    % Plot the average power across IAF channels
                    subplot(2, 1, 2);
                    imagesc(temp.times{band_idx}, 1:num_trials, squeeze(avg_power)');
                    axis xy;
                    xlabel('Time (ms)');
                    ylabel('Trials');
                    title(['Average Power Across IAF Channels - Band: ', freq_bands_select{band_idx, 1}]);
                    colorbar;
                    
                end
            end

            temp.data = tf_data;

            if S.prep.calc.eeg.TF.plot_allbands

                % Option to display values relative to baseline
                baseline_correct = true;
                
                % New figure to plot data from all frequency bands on one plot using imagesc
                figure;
                
                % Initialize a matrix to hold the averaged data for all frequency bands
                all_band_data = [];
                
                % Collect all time points from different bands into a single array
                all_times = [];
                for i = 1:length(temp.times)
                    all_times = [all_times; temp.times{i}(:)];
                end
                
                % Find the unique and sorted time points
                unique_times = unique(all_times);
                
                % Find the maximum number of time points across all bands for uniformity
                max_time_points = length(unique_times);
                
                % Initialize time vector for the x-axis (uniform across all bands)
                uniform_times = linspace(min(unique_times), max(unique_times), max_time_points);
                
                for band_idx = 1:size(freq_bands_select, 1)
                    % Average the data over trials and selected IAF channels
                    avg_tf_data_band = mean(mean(tf_data{band_idx}(channels, :, :), 1), 3);
                    
                    % Check if there are at least two time points for interpolation
                    if length(temp.times{band_idx}) >= 2
                        % Interpolate to ensure uniform time points across all bands
                        interpolated_data = interp1(temp.times{band_idx}, 10*log10(squeeze(avg_tf_data_band)), uniform_times, 'linear', 'extrap');
                    else
                        % If there are fewer than two points, replicate the available data
                        interpolated_data = repmat(10*log10(squeeze(avg_tf_data_band)), 1, max_time_points);
                    end
                    
                    % Ensure interpolated_data is not empty
                    if isempty(interpolated_data)
                        interpolated_data = NaN(1, max_time_points);
                    end
                    
                    % Baseline correction if the option is enabled
                    if baseline_correct
                        baseline_indices = find(uniform_times >= -200 & uniform_times <= 0);
                        if ~isempty(baseline_indices)
                            baseline_mean = mean(interpolated_data(baseline_indices));
                            interpolated_data = interpolated_data - baseline_mean;
                        else
                            warning('No baseline period found for band %s', freq_bands_select{band_idx, 1});
                        end
                    end
                
                    % Append the interpolated data to the all_band_data matrix
                    all_band_data = [all_band_data; interpolated_data];
                end
                
                % Ensure the all_band_data matrix is the correct size for imagesc
                if size(all_band_data, 1) == size(freq_bands_select, 1)
                    % Plot the data using imagesc
                    imagesc(uniform_times, 1:size(freq_bands_select, 1), flipud(all_band_data));
                    set(gca, 'YTick', 1:size(freq_bands_select, 1), 'YTickLabel', flipud({freq_bands_select{:, 1}}'));
                    xlabel('Time (ms)');
                    ylabel('Frequency Band');
                    title('Average Power Across IAF Channels for All Frequency Bands');
                    c = colorbar;
                    ylabel(c, 'Power (dB)'); % Label the colorbar
                else
                    error('The size of all_band_data does not match the number of frequency bands selected.');
                end








            end



        end
        
        for fn = 1:length(fename)
            D(d).prep.(fename{fn})(f).dat = temp.(fename{fn});
        end

        % get other filetypes with same name
        for fn = 1:length(S.prep.select.suffixes)-1
            loadsuffix = S.prep.select.suffixes{fn+1};
            filename = strrep(filename,S.prep.select.suffixes{1},loadsuffix);

            temp = load(filename);
            fename = fieldnames(temp(f));
            for fn = 1:length(fename)
                if isstruct(temp.(fename{fn}))
                    D(d).prep.(fename{fn})(f) = temp.(fename{fn});
                else
                    D(d).prep.(fename{fn})(f).dat = temp.(fename{fn});
                end
            end
        end
    end

    % get data information
    switch S.prep.fname.ext{:}
        case 'mat'
            if isfield(S.prep.select,'freq')
                timecourse = squeeze(D(d).prep.fdata.dat{1, 1}.powspctrm(:,:,S.prep.select.freq,:));
                n_chans = length(D(d).prep.fdata.dat{1, 1}.label);
                n_samples = length(D(d).prep.fdata.dat{1, 1}.time{1});
                n_trials = length(D(d).prep.fdata.dat{1, 1}.time);
                timecourse = reshape(permute(timecourse,[2 3 1]),n_chans,[]);
                eventType = D(d).prep.fdata.dat{1, 1}.conds; tnums = D(d).prep.fdata.dat{1, 1}.tnums; fnums = D(d).prep.fdata.dat{1, 1}.fnums; bnums = D(d).prep.fdata.dat{1, 1}.bnums;
            elseif isfield(S.prep.fname,'struct')
                eeg = D(d).prep.(fename{1})(1).dat.(S.prep.fname.struct{2}); 
                n_chans = size(eeg,1);
                n_samples = size(eeg,2);
                n_trials = size(eeg,3);
                eventType = ones(1,n_trials);
                tnums=1:length(eventType); % assumes all trials present
                D(d).prep.dim = [n_chans,n_samples,n_trials];
                timecourse = reshape(eeg,n_chans,[]); % make 2D: chans x (time*trials)
                % save other info
                other = fieldnames(D(d).prep.(fename{1})(1).dat);
                other = setdiff(other,S.prep.fname.struct{2});
                for fo = 1:length(other)
                    D(d).prep.other.(other{:}) = D(d).prep.(fename{1})(1).dat.(other{:});
                end
            else % this is for processing group ICA data files
                try
                    topography = D(d).prep.topography.dat;
                end
                timecourse = D(d).prep.timecourse.dat;
                eventType = D(d).prep.eventType.dat;
                n_chans = D(d).prep.n_chans.dat;
                n_samples = D(d).prep.n_samples.dat;
                n_trials = D(d).prep.n_trials.dat;
                tnums=1:length(eventType); % assumes all trials present
            end
        case 'set'
            % For data recorded from an EGI system with STIM markers for
            % timing events
            if any(find(strcmp('STIM',temp.epoch(1).eventtype)))
                [eventType,tnums, fnums, bnums] = get_markers(temp,'EGI');
            else % from other systems, e.g. Brainamps
                [eventType,tnums, fnums, bnums] = get_markers(temp,'BV');
            end
            
            % This is a study-specific patch for CORE - can remove
            % eventually.
            % Make correction for CORE006, part2 data
                % block 2 tnums incorrectly restart from 0, when should start
                % from 452
            if strcmp(D(d).prep.subname,'CORE006') && length(tnums)<1000 % i.e. it is part2 data
                tnums(bnums==2)=tnums(bnums==2)+451;
            end

            % 
            n_chans = D(d).prep.nbchan.dat;
            n_trials = D(d).prep.trials.dat;

            if iscell(D(d).prep.data.dat)
                for band_idx = 1:size(freq_bands_select, 1)
                    n_samples = D(d).prep.pnts.dat{band_idx};
                    D(d).prep.dim{band_idx} = [n_chans,n_samples,n_trials];
                    timecourse{band_idx} = reshape(D(d).prep.data.dat{band_idx},n_chans,[]); % make 2D: chans x (time*trials)
                end
            else
                n_samples = D(d).prep.pnts.dat;
                D(d).prep.dim = [n_chans,n_samples,n_trials];
                timecourse = reshape(D(d).prep.data.dat,n_chans,[]); % make 2D: chans x (time*trials)
            end
    end
    D(d).prep.tnums = tnums;
    
    if exist('topography','var')
        comps = 1:size(timecourse,1);
        data_tc= topography(:,comps) * timecourse(comps,:);  
    else
        data_tc= timecourse;  
    end
    if isfield(D(d).prep,'times')
        total_samples = D(d).prep.times.dat;
    else
        total_samples=S.prep.original_samples;
    end
    if ~isempty(S.prep.select.samples)
        select_samples=S.prep.select.samples;
    else
        select_samples=total_samples;
    end
    
    
    %% Predictors: factors from EEG markers - ASSUMES THEY ARE NUMERIC IN EEG DATA FILE but this can be updated if needed
    % create predictor variables
    dtab=table; % initialise data table for this subject
    
    if ~isempty(S.prep.pred.factor.markers) 
        factor_markers=S.prep.pred.factor.markers;
        
        if ~iscell(factor_markers{1})
            factor_markers = {factor_markers}; % for backwards compatability with some old scripts
        end
        for fac = 1:length(factor_markers)
        
            % for each level of the factor
            for fac_lev = 1:length(factor_markers{fac})
                condidx{fac_lev}=[];
                % for each marker, get indices over EEG trials
                for marker = 1:length(factor_markers{fac}{fac_lev})
                    condidx{fac_lev} = [condidx{fac_lev} find(ismember(eventType,factor_markers{fac}{fac_lev}(marker)))];
                end
                % pool into conditions defined by factor_markers, indexing
                % from zero
                predtemp(condidx{fac_lev},1)=fac_lev-1;
            end
            if any(S.prep.pred.factor.ordinal==fac)
                pred = categorical(predtemp,unique(predtemp),'ordinal',1);
            else
                pred = categorical(predtemp,unique(predtemp));
            end
            if isfield(S.prep.pred.factor,'label') && ~isempty(S.prep.pred.factor.label)
                pred_label=S.prep.pred.factor.label{fac};
            else
                pred_label=['fac' num2str(fac)];
            end
            dtab.(pred_label) = pred;
            clear predtemp predlabel
        end
    end
    
    %% select factor levels
    orig_tnums = tnums;
    if ~isempty(S.prep.pred.factor.select)
        for f = 1:size(S.prep.pred.factor.select,1)
            pred_label = S.prep.pred.factor.select{f,1};
            disp(['selecting conditions from factor ' pred_label])
            fac_levels = categorical(S.prep.pred.factor.select{f,2});
            dat = categorical(dtab.(pred_label));
            idx = ismember(dat,fac_levels);
            dtab = dtab(idx,:);
            tnums = tnums(idx);
            fnums = fnums(idx);
            bnums = bnums(idx);
            eventType = eventType(idx);
            if length(unique(dtab.(pred_label)))==1
                dtab.(pred_label) =[];
            elseif length(unique(dtab.(pred_label)))==2
                dtab.(pred_label) = categorical(dtab.(pred_label)); % remove ordinality
                dtab.(pred_label) = removecats(dtab.(pred_label));
            end
        end
    end
    D(d).prep.tnums = tnums;
    
    
    %% EEG data operations / transformations / calculations
    % this is done after defining predictors in case subtractions are
    % required

    if iscell(data_tc)
        idx = 1:size(freq_bands_select, 1);
        D(d).prep.dim = {};
        D(d).prep.data = {};
    else
        idx=1;
    end

    for i=idx

        if iscell(data_tc)
            data = data_tc{i};
        else
            data = data_tc;
        end
        if iscell(total_samples)
            tsamples = total_samples{i};
        else
            tsamples = total_samples;
        end
        if iscell(select_samples)
            ssamples = select_samples{i};
        else
            ssamples = select_samples;
        end
    
        % smooth
        if S.prep.calc.eeg.smooth_samples
            disp('smoothing...')
            for i = 1:size(data,1)
                data(i,:) = smooth(data(i,:),S.prep.calc.eeg.smooth_samples,'moving');
            end
        end
    
        % downsample
        if S.prep.calc.eeg.dsample
            orig_samples = tsamples;
            tsamples = downsample(tsamples',S.prep.calc.eeg.dsample)';
            ssamples = downsample(ssamples',S.prep.calc.eeg.dsample)';
            idx = repmat(ismember(orig_samples,tsamples),1,n_trials);
            data = data(:,idx);
        end
        
        % make 3D
        data=reshape(data,size(data,1),[],n_trials);
        
        % select data
        data = data(:,:,ismember(orig_tnums,tnums));
        n_trials = size(data,3);
    
        % flip channels right to left
        if ~isempty(S.prep.calc.eeg.flipchan)
            load(S.prep.path.inputs.chanlocs);
            flipidx = ismember(eventType,S.prep.calc.eeg.flipchan);
            data(:,:,flipidx) = flipchan(data(:,:,flipidx),chanlocs,S.prep.path.inputs.GSNlocs);
        end
        
        % subtractions
    %     if ~isempty(S.prep.pred.factor.subtract)
    %         ...
    %     end
        
        % Data transformation
        if ~isempty(S.prep.calc.eeg.transform)
            x=data;
            if strcmp(S.prep.calc.eeg.transform,'arcsinh')
                data=log(x+sqrt(x.^2+1));
            elseif strcmp(S.prep.calc.eeg.transform,'log')
                data=log(x);
            elseif strcmp(S.prep.calc.eeg.transform,'pow')
                data=10*log10(x);
            end
        end
        
        % baseline correct
        if ~isempty(S.prep.calc.eeg.baseline_correct)
            bp = dsearchn(tsamples',S.prep.calc.eeg.baseline_correct');
            data = bsxfun(@minus, data, mean(data(:,bp(1):bp(2),:),2));
        end
    
        % select samples
        if exist('select_samples','var') && ~isempty(ssamples)
            select = dsearchn(tsamples',[ssamples(1),ssamples(end)]');
            data = data(:,select(1):select(2),:);
            total_samples{i} = total_samples{i}(select(1):select(2));
        end
        
        % trim
    %     S.prep.calc.eeg.ndec=8; % trim data to a number of decimal places
            
 
        if S.prep.calc.eeg.TF.zscore

            % Get the size of the input matrix
            [dim1, dim2, dim3] = size(data);
            
            % Reshape the 3D matrix into a 2D matrix
            data = reshape(data, [], 1);
            
            % Compute the mean and standard deviation of the original data
            mean_data = mean(data);
            std_data = std(data);
            
            % Z-score the data
            data = (data - mean_data) / std_data;
            
            % Reshape the z-scored 2D matrix back to 3D matrix
            data = reshape(data, dim1, dim2, dim3);

        elseif S.prep.calc.eeg.TF.norm

            data = data/norm(data(:));
        end

        %% reshape data into samples and combine over subjects

        if iscell(data_tc)
            D(d).prep.dim{i} = size(data);
            D(d).prep.data{i}=reshape(data,[],n_trials)';
            if S.prep.calc.eeg.TF.zscore
                D(d).prep.other.data_mean{i}=mean_data;
                D(d).prep.other.data_std{i}=std_data;
            end
        else
            D(d).prep.dim = size(data);
            D(d).prep.data=reshape(data,[],n_trials)';
        end


        clear data % to free memory

    end
    
     

    if ~S.prep.calc.eeg.TF.save_as_cell
        data = [];
        for band_idx = 1:size(freq_bands_select, 1)
            data = cat(2,data,D(d).prep.data{band_idx});
        end
        D(d).prep.data=data;
        dim(1)= D(d).prep.dim{1}(1);
        dim(2)= sum(cellfun(@(x) x(2), D(d).prep.dim));
        dim(3)= D(d).prep.dim{1}(3);
        freq_pnts = cellfun(@(x) x(2), D(d).prep.dim);
        D(d).prep.dim = dim;
    end

    %% separate EEG/predictor data into training and testing fractions, only for testing decoders
    trainfrac=S.prep.traintest.frac;
    switch S.prep.traintest.type
        case 'random'
            % Split into training (encoding) and testing (decoding) fractions
            % Produces S.trainidx and S.testidx. These are indices
            % of conData-ordered trials that will be used for training and
            % testing.
            trainidx = {[]};
            testidx = {[]};
            if S.prep.traintest.balance_conds
                ucond = unique(eventType);
                for u = 1:length(ucond)
                    cond_idx = find(eventType==ucond(u)); % for each condition, it's index within conData
                    % DO THIS FOR A NUMBER OF RUNS - MULTIPLE RANDOMISED INDICES
                    for run = 1:S.prep.traintest.num_runs
                        rng(run); % seed random number generator for consistency
                        trainidx{run} = [trainidx{run} randsample(cond_idx,round(trainfrac*length(cond_idx)))]; 
                    end
                end
            else
                % DO THIS FOR A NUMBER OF RUNS - MULTIPLE RANDOMISED INDICES
                for run = 1:S.prep.traintest.num_runs
                    rng(run); % seed random number generator for consistency
                    trainidx{run} = randsample(length(eventType),round(trainfrac*length(eventType)));
                end
            end
            % create test index
            train_matrix = zeros(length(tnums),length(trainidx));
            test_matrix=zeros(length(tnums),length(testidx));
            for run = 1:length(trainidx)
                trainidx{run} = sort(trainidx{run});
                if trainfrac<1
                    testidx{run} = 1:length(eventType);
                    testidx{run}(trainidx{con}) = [];
                else
                    % otherwise, duplicate - same trials for each
                    testidx{run} = trainidx{run};
                end
                % convert to design matrix
                train_matrix(trainidx{run},run) = 1;
                test_matrix(testidx{run},run) = 1;
            end

    end
    
    % store info in table format
    D(d).prep.dtab=dtab;
    D(d).prep.dtab.ID=repmat(D(d).prep.subname,length(tnums),1);
    D(d).prep.dtab.group=repmat(designmat.groups(d),length(tnums),1);
    D(d).prep.dtab.eventTypes=categorical(eventType');
    if exist('fnums','var'); D(d).prep.dtab.fnums=categorical(fnums'); end
    if exist('bnums','var'); D(d).prep.dtab.bnums=categorical(bnums'); end
    if exist('tnums','var'); D(d).prep.dtab.tnums=tnums'; end
    D(d).prep.dtab.train=categorical(train_matrix);
    D(d).prep.dtab.test=categorical(test_matrix);
    if exist('iaf','var')
        D(d).prep.other.IAF=iaf;
    end
    if exist('freq_bands_select','var')
        D(d).prep.other.freq=freq_bands_select;
    end
    if exist('freq_pnts','var')
        D(d).prep.other.freq_pnts=total_samples;
    end
    
    %D(d).prep.tnums = tnums;
    
end

% clean up D
E=struct;
for d = 1:length(D)
    E(d).prep.subname = D(d).prep.subname;
    E(d).prep.data = D(d).prep.data;
    E(d).prep.dtab = D(d).prep.dtab;
    E(d).prep.dim = D(d).prep.dim;
    if isfield(D(1).prep,'other')
        E(d).prep.other = D(d).prep.other;
    end
    % if exist('iaf','var')
    %     E(d).prep.IAF=iaf;
    % end
    % if exist('freq_bands_select','var')
    %     E(d).prep.freq=freq_bands_select;
    % end
    % if exist('freq_pnts','var')
    %     E(d).prep.freq_pnts=freq_pnts;
    % end
    %E(d).prep.tnums = D(d).prep.tnums;
end
D=E;
clear E

if S.prep.output.save
    disp('saving data to disk...')
    save(fullfile(S.prep.path.outputs,[S.prep.sname '.mat']),'D','S','-v7.3')
    disp('...done')
end

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
end


function iaf = calculateIAF(data, srate, channels, plot_psd)
    % Default to no plotting if not provided
    if nargin < 4
        plot_psd = false;
    end

    % Subtract the average event-related response from each trial to obtain induced responses
    mean_response = mean(data, 3);
    induced_responses = data - mean_response;

    % Calculate power spectral density (PSD) for each specified channel
    num_channels = length(channels);
    psd_all = [];

    % Multitaper parameters
    params.Fs = srate;
    params.tapers = [3 5]; % Time-bandwidth product and number of tapers

    % Loop through each specified channel and calculate PSD using multitaper method
    for ch = 1:num_channels
        channel_data = reshape(induced_responses(channels(ch), :, :), 1, []);
        [S, f] = mtspectrumc(channel_data, params);
        psd_all = [psd_all; S']; % Collect PSD for each channel
    end
    
    % Compute median PSD across specified channels
    median_psd = median(psd_all, 1);
    
    % Apply Gaussian smoothing to the PSD
    smooth_psd = smoothdata(median_psd, 'gaussian', max(5,round(length(median_psd)*0.01))); % Adjust the window size if needed

    % Limit frequency range to 0-20 Hz
    freq_limit = 20;
    freq_indices = f <= freq_limit;
    f = f(freq_indices);
    median_psd = median_psd(freq_indices);
    smooth_psd = smooth_psd(freq_indices);
    
    % Define alpha range
    alpha_range = [6 14];
    alpha_indices = find(f >= alpha_range(1) & f <= alpha_range(2));
    
    % Find peak frequency within alpha range for both unsmoothed and smoothed PSD
    [~, max_idx_unsmoothed] = max(median_psd(alpha_indices));
    iaf_unsmoothed = f(alpha_indices(max_idx_unsmoothed));
    
    [~, max_idx_smoothed] = max(smooth_psd(alpha_indices));
    iaf_smoothed = f(alpha_indices(max_idx_smoothed));
    
    % Output the smoothed IAF
    iaf = iaf_smoothed;
    
    % Optional plotting of the PSD
    if plot_psd
        figure;
        plot(f, 10*log10(median_psd), 'b');
        hold on;
        plot(f, 10*log10(smooth_psd), 'r');
        plot(iaf_unsmoothed, 10*log10(median_psd(alpha_indices(max_idx_unsmoothed))), 'bo');
        plot(iaf_smoothed, 10*log10(smooth_psd(alpha_indices(max_idx_smoothed))), 'ro');
        xlabel('Frequency (Hz)');
        ylabel('Power Spectral Density (dB)');
        title(['Power Spectral Density with Identified IAF: ' num2str(iaf_smoothed)]);
        legend('Unsmoothed PSD', 'Smoothed PSD', 'Unsmoothed IAF', 'Smoothed IAF');
        xlim([0 freq_limit]);
        hold off;
    end
end






