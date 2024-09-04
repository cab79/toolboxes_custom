function eegstats_condor_move_resubmit

pth=pwd;
pthData=fullfile(pwd,'Data');
load('resubmit_files.mat');

% copy input and output files to data folder while re-naming.
for i = 1:length(output_filenames_resubmit_new)
    if exist(fullfile(pth,output_filenames_resubmit_new{i}),'file')
        copyfile(fullfile(pth,output_filenames_resubmit_new{i}),fullfile(pthData,'output_files',output_filenames_resubmit_original{i}));
    end
    % if exist(fullfile(pth,input_filenames_resubmit_new{i}),'file')
    %     delete(fullfile(pth,input_filenames_resubmit_new{i}));
    % end
end

% remove them back
movefile(fullfile(pthData,'output_files','output*'),pth);
movefile(fullfile(pthData,'input_files','input*'),pth);

quit