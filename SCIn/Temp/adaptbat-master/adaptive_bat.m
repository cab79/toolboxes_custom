function adaptive_bat(subject, type, ptb, params)

% Adaptive threshold version of the Beat Alignment Test.

% Subjects are presented musical stimuli superimposed by a metronome. The
% metronome is either aligned with the beat or perturbed in phase or tempo
% with varying magnitude from trial to trial. 20% of trials are catch
% trials randomly distributed throughout the test. 50% of test trials are
% negative deviations and 50% are positive (e.g., if tempo, half are reduced
% tempo and half are increased tempo).
%
% Perturbation magnitude is determined trial-to-trial by the ZEST threshold
% algorithm. Trial presentation ends when the threshold stopping criteria
% are met and the threshold converges on a final estimate.
%
% The most efficient way to run this function is to call it from the
% wrapper script BAT.m
%
% REQUIRED INPUT:
%   * subject: subject ID string (e.g., 'JR1')
%   * type: deviant type string. 'B' for ibi (tempo) test; 'P' for phase test
%   * ptb: structure array containing PsychToolbox parameters
%       (ptb params created by the function bat_initialize_ptb.m)
%   * params: structure array with the following fields:
%       - stim_fpath: directory path where stim files live (string)
%       - stim_names: names of all stims (cell of strings)
%       - data_fpath: directory path where data files will be written
%       - outdata_fname: name for data file
%       - zest: structure array containing parameters that control behavior
%           of the ZEST algorithm
%       - stim_names: cell array of excert names
%       - repetitions: number of times the excerpt should repeat upon
%           playback
%
%
% Written by Brian K. Hurley, Sept 2016.
%

%% TASK/TRIAL HANDLING
clear mean_pdf sd_pdf

max_offbeat_trials = params.zest.max_trials; % doesn't include onbeat catch trials
offbeat_trial_counter = 0;

% Initialize data table for output.
% Initializing with more than enough rows. At the end of this function we
% trim of empty rows.
out_data_vars = {'response', 'aug_cond', 'deviation', 'score', ...
    'thresh', 'sd_pdf', 'converged', 'excerpt'};
resp_tbl = cell2table(cell(40,length(out_data_vars)), ...
    'VariableNames', out_data_vars);

% Initialize ZEST
% if (second) params input given, zest ignores (first) response input
% position and intializes zest...hence NaN
mean_pdf = zest(NaN, params.zest);

for trial_number = 1:max_offbeat_trials
    
    % select random song index
    k = randperm(12,1);
    
    % first trial starts with on-beat stim, subsequent trials, probability
    % of off-beat stim is 0.80
    switch trial_number
        case 1 % for first trial
            song = params.stim_names{k}(1:3);
            aug_cond = 'onbeat';
        otherwise % all subsequent trials
            % Get stimulus and make sure it isn't a repeat.
            % If current selection matches previous trial, re-select
            % until we have a new song. Insures we never have back-to-back repeats.
            while strcmp(resp_tbl.excerpt{trial_number-1}, song)
                k = randperm(12,1);
                song = params.stim_names{k}(1:3);
            end
            
            % use RAND to implement 0.80 probability of getting an off-beat trial.
            % This will randomize placement of catch trials
            x = rand;
            if x < 0.80
                aug_cond = 'offbeat';
                % add this trial to offbeat trial counter
                offbeat_trial_counter = offbeat_trial_counter + 1;
            else
                aug_cond = 'onbeat';
            end
            
    end % switch trial_number
    
    % for on-beat trials, generate stim w/ on-beat metronome
    if strcmp(aug_cond,'onbeat')
        deviation = 0;
        stim_wvf = mkbat3(deviation, type, song, params);
        stim_wvf = stim_wvf';
        
        % not tracking thresh metrics on catch trials
        resp_tbl.thresh{trial_number} = NaN;
        resp_tbl.sd_pdf{trial_number} = NaN;
        resp_tbl.converged{trial_number} = NaN;
        
    else % else, generate stim w/ deviation specified by threshold algorithm
        
        % randomize whether deviation is positive or
        % negative
        dev_sign = randperm(2,1);
        switch dev_sign
            case 1
                deviation = mean_pdf;
            case 2
                deviation = -mean_pdf;
        end
        
        % generate stimulus
        stim_wvf = mkbat3(deviation, type, song, params);
        stim_wvf = stim_wvf';
        
    end
    
    %% PLAY AUDIO AND GET RESPONSE
    
    % Wait 0.5 s before starting
    WaitSecs(0.5);
    
    % Fill the audio playback buffer with the audio data
    PsychPortAudio('FillBuffer', ptb.pahandle, [stim_wvf; stim_wvf]);
    
    % Start audio playback
    PsychPortAudio('Start', ptb.pahandle, params.repetitions, ptb.startCue, ...
        ptb.waitForDeviceStart);
    
    % Collect keyboard response: on beat/YES = p, off beat/NO = q
    respToBeMade = true;
    while respToBeMade == true
        [~,~, keyCode] = KbCheck;
        if keyCode(ptb.YES)
            response = 'yes';
            if strcmp(aug_cond, 'onbeat')
                score = 'correct';
            else
                score = 'incorrect';
            end
            respToBeMade = false;
        elseif keyCode(ptb.NO)
            response = 'no';
            if strcmp(aug_cond, 'offbeat')
                score = 'correct';
            else
                score = 'incorrect';
            end
            respToBeMade = false;
        elseif keyCode(ptb.escapeKey)
            sca;
        end
    end
    
    %% Populate data table and update ZEST algorithm for next trial (or final threshold)
    
    % populate this trial in response table
    resp_tbl.response{trial_number} = response;
    resp_tbl.aug_cond{trial_number} = aug_cond;
    resp_tbl.deviation{trial_number} = deviation;
    resp_tbl.score{trial_number} = score;
    resp_tbl.excerpt{trial_number} = song;
    
    % if it's an off-beat trial, estimate threshold
    % for the ZEST function, 1 is correct response, 0 is incorrect response
    if strcmp(aug_cond, 'offbeat')
        if strcmp(response, 'no')
            zest_resp = 1;
        elseif strcmp(response, 'yes')
            zest_resp = 0;
        end
        
        % estmate thresh
        [mean_pdf, sd_pdf] = zest(zest_resp);
        
        resp_tbl.thresh{trial_number} = mean_pdf;
        resp_tbl.sd_pdf{trial_number} = sd_pdf;
        
        % check for convergence.
        % if converged, stop the trial_number loop (i.e. stop delivering stims)
        % we won't let it converge until at least 7 trials have passed
        if max_offbeat_trials >= 7 && sd_pdf <= params.zest.sd_stop || ...
            offbeat_trial_counter == max_offbeat_trials
            resp_tbl.converged{trial_number} = 1;
            break
        else
            resp_tbl.converged{trial_number} = 0;
        end
    end
    
end % for trial_number

% since we can't predict ahead of time which trial the threshold will
% converge on, we need to trim empty rows from the data table.
empty_row_idx = cellfun(@isempty,resp_tbl.response);
resp_tbl(empty_row_idx,:) = [];

% add column for subject_id
resp_tbl.subject_id = cell(length(resp_tbl.response),1);
resp_tbl.subject_id(:) = deal({subject});

%% OUTPUT RESULTS TO FILE

% write data to file
out_fname = fullfile(params.data_fpath, sprintf(['%s_' params.outdata_fname], subject));
writetable(resp_tbl, out_fname);

end
