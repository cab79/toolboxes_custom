% Code to identify non-completed outputs, re-submit them, and integrate with
% old outputs.
% Needed in the situation when Condor jobs are stuck in "idle", e.g. due to
% computers going to sleep.

%% step 1
% remove current job from Condor

%% step 2: Run this Matlab code (eegstats_condor_resubmit) from Linux script

% create variables listing input files, output files that don't exist but
% have input files, and output files that exist but are zero sized. Save
% this info in a file as a record.

% move all completed input and output files to another folder. delete any remaining output files. re-name
% input files consecutively and keep reference information relative to
% original index - save in file.

% end of Matlab function

%% Step 3: Linux script
% submit job

%% Step 4: Run new Matlab function (eegstats_condor_monitor_resubmit)
% Function runs load_condor_data and monitors outputs.
% when complete, loads save file with index of output file names. Re-names
% outputs, then moves the original outputs completed back to main folder.

%% Step 5: Continue compiling outputs.

function eegstats_condor_resubmit

S.encode.method = 'LMM';

pth=pwd;
pthData=fullfile(pwd,'Data');

if ~exist(fullfile(pthData,'input_files'),'dir')
    mkdir(fullfile(pthData,'input_files'))
else
    delete(fullfile(pthData,'input_files','*.mat'));
end
if ~exist(fullfile(pthData,'output_files'),'dir')
    mkdir(fullfile(pthData,'output_files'))
else
    delete(fullfile(pthData,'output_files','*.mat'));
end

% create variables listing input files, output files that don't exist but
% have input files, and output files that exist but are zero sized. 
input_files = dir(fullfile(pth,'input*.mat'));
input_filenames = {input_files(:).name};
output_files = dir(fullfile(pth,'output*.mat'));
output_filenames = {output_files(:).name};
expected_output_filenames = regexprep(input_filenames, 'input', 'output');
output_filenames_missing = setdiff(expected_output_filenames,output_filenames);
output_filenames_zero = output_filenames([output_files(:).bytes]==0);
output_filenames_resubmit_original = sort(horzcat(output_filenames_missing,output_filenames_zero));
output_filenames_completed = setdiff(expected_output_filenames,output_filenames_resubmit_original);
input_filenames_resubmit_original = regexprep(output_filenames_resubmit_original, 'output', 'input');
input_filenames_resubmit_new={};
for i = 1:length(input_filenames_resubmit_original)
    input_filenames_resubmit_new{i} = ['input' num2str(i-1) '.mat'];
end
output_filenames_resubmit_new = regexprep(input_filenames_resubmit_new, 'input', 'output');

% Save this info in a file as a record.
save('resubmit_files.mat','input_filenames_resubmit_original','input_filenames_resubmit_new','output_filenames_resubmit_original','output_filenames_resubmit_new')

% move all input and completed output files to another folder. delete any remaining output files. 
movefile('input*.mat',fullfile(pthData,'input_files'));
for i = 1:length(output_filenames_completed)
    movefile(output_filenames_completed{i},fullfile(pthData,'output_files'));
end
delete('output*.mat');

% copy input files to main folder while re-naming consecutively.
for i = 1:length(input_filenames_resubmit_original)
    copyfile(fullfile(pthData,'input_files',input_filenames_resubmit_original{i}),fullfile(pth,input_filenames_resubmit_new{i}));
end

% create_job_submission_file
methodfunc=['eegstats_encode_' S.encode.method];
pth=pwd;
nf=length(input_filenames_resubmit_original);
disp('creating job submission file')
A = {
    ['executable=' methodfunc '.exe']
    'indexed_input_files=input.mat'
    'indexed_output_files=output.mat'
    ['indexed_stdout=' methodfunc '.out']
    ['indexed_stderr=' methodfunc '.err']
    ['indexed_log=' methodfunc '.log']
    'runtime=240'
    ['total_jobs=' num2str(nf)]
};

fid = fopen(fullfile(pth, [methodfunc '_run.sub']),'w');
for i = 1:size(A,1)
    fprintf(fid,'%s\n',A{i});
end
fclose(fid);

% end of Matlab function
quit


