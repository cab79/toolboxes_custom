function outData = adaptive_bat_load_data(~,params)

% imports subjects' individual data files and combines concatenates into 
% single data table.
% called from adaptive_bat_analyses.m
% 
% 11 Oct 2016 BKH

% find csv file names in data directory
batData_path = fullfile(params.data_fpath,'*.csv');
batData = dir(batData_path);
batData_fnames = {batData.name};
nsubs = length(batData_fnames);

% cycle through subject data & concatenate to single data table
sub_data = cell(nsubs,1);
sub_thresh = cell(nsubs,1);
for isub = 1:nsubs
    
    % load subject data as table
    curr_fname = fullfile(params.data_fpath,batData_fnames{isub});
    this_sub_data = readtable(curr_fname);
    
    % if there is no trial_num variable in table, add trial_num column
    if ~any(strcmp(this_sub_data.Properties.VariableNames, 'trial_num'))
        n_trials = length(this_sub_data.response);
        this_sub_data.trial_num = (1:n_trials)';
    end
    
    % if there is no subject variable in table, add subject column using file name
    if ~any(strcmp(this_sub_data.Properties.VariableNames,'subject_id'))
        
        % find whether this is an IBI or phase data file
        ibi_fstub_idx =  strfind(batData_fnames{isub},params.ibi.outdata_fname);
        phase_fstub_idx = strfind(batData_fnames{isub},params.phase.outdata_fname);
        % get subject name from file name
        if ~isempty(ibi_fstub_idx)
            subject_id = batData_fnames{isub}(1:ibi_fstub_idx-2);
        else
            subject_id = batData_fnames{isub}(1:phase_fstub_idx-2);
        end  
        
        this_sub_data.subject_id = cell(size(this_sub_data,1),1);
        this_sub_data.subject_id(:) = deal({subject_id});
        
    end 
    
    % if there is no 'test' variable in table (i.e. was this a tempo test or 
    % a phase test?), add 'test' column as determined by data file name
    if ~any(strcmp(this_sub_data.Properties.VariableNames,'test'))
        
        % find whether this is an IBI or phase data file and populate
        % test column accordingly
        ibi_fstub_idx =  strfind(batData_fnames{isub},params.ibi.outdata_fname);
        phase_fstub_idx = strfind(batData_fnames{isub},params.phase.outdata_fname);
        if ~isempty(ibi_fstub_idx)
            test = 'tempo';
        elseif ~isempty(phase_fstub_idx)
            test = 'phase';
        end        
        this_sub_data.test = cell(size(this_sub_data,1),1);
        this_sub_data.test(:) = deal({test});
        
    end 
    
    % cell array of each subject's threshold data table
    sub_thresh{isub} = this_sub_data(this_sub_data.converged == 1,:);    
    
    % cell array of each subject's trial data table
    sub_data{isub} = this_sub_data;
    
end % for isub

% concatenate subject data tables to single table for (1) trials and (2) thresholds
outData.vars{1} = 'adaptBAT_trial_data';
outData.data{1} = vertcat(sub_data{:});

outData.vars{2} = 'adaptBAT_thresh_data';
outData.data{2} = vertcat(sub_thresh{:});

end