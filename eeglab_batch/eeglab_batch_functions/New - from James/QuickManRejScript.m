%% Quick script
% Run this to have quick access to the manual rejection step of data
% analysis

% load existing S structure
set_path = 'Y:\James Henshaw\WP3 Analysis\Preprocessed'; % WP3 set path
load(fullfile(set_path,'S'))

% add paths
addpath(S.fieldtrip_path); ft_defaults;
addpath(S.eeglab_path); eeglab; 
delete(findall(0,'Type','figure')) % close all open figur
addpath(genpath(S.batchfun_path));
addpath(genpath(S.support_path));

fprintf('DEBUG MODE ACTIVE...\n\n...STARTING NOW\n')   % DELETE
sname_ext = 'combined';                                          % DELETE

% separate into single and double figures
partnum = [1:10];
ind_sing = partnum<10; ind_dbl = partnum>=10&partnum<100; % handles single and double figures IDs
singleCell = strseq('P00',partnum(ind_sing));
dblCell = strseq('P0',partnum(ind_dbl));

% create cell of subject IDs
strCell = [singleCell;dblCell;{}];
S.subjects = strCell;
S.conds = {'EO1','EO2','EC1','EC2','NFB1','NFB2'};              % DELETE
S.blocks = {'B1','B2','B3'};                                    % DELETE

% Choose whether to overwrite or ignore files if one exists
S.overwriteManrej = 0; % overwrite manual rejection (0 = no, 1 = yes)
S.overwriteICA = 0; % overwrite ICA file (0 = no, 1 = yes)

% NOISY TRIAL AND CHANNEL REJECTION USING FIELDTRIP
if iscell(S.FTrej)
    
    % GET FILE LIST
    S.filepath = fullfile(S.setpath,sname_ext);
    if strcmp(sname_ext,'combined')
        S.blocks={};S.conds={};
    end
    S = getfilelist(S,sname_ext);
    
    loadpath = fullfile(S.setpath,sname_ext);
    
    
    global EEG
    for f = 1:length(S.filelist)
        
        file = S.filelist{f};
        
        % check if file exists already
        [~,nme,~] = fileparts(file);
        checker = fullfile(S.setpath,'manrej',sprintf('%s*',nme));
        
        if isempty(dir(checker)) || S.overwriteManrej % skip existing files, depending on settings
            [ALLEEG EEG CURRENTSET ALLCOM] = eeglab; % open new eeglab session (allows pop_eegplot to function properly!)
            
            EEG = pop_loadset('filename',file,'filepath',loadpath);
            [EEG,changes] = eeg_checkset(EEG);
            
            % reject channels
            [~, indElec] = pop_rejchan(EEG, 'elec',[1:EEG.nbchan] ,'threshold',5,'norm','on','measure','kurt');
            
            EEG = manualRejectEEG(EEG);
            
            % strategy - only remove a very small number of very bad trials / chans
            % before ICA - do further cleaning after ICA
            for i = 1:length(S.FTrej)
                EEG = FTrejman(EEG,S.FTrej{i});
            end
            
            % SAVE
            [pth nme ext] = fileparts(file);
            sname_ext = 'manrej';
            sname = [nme '_' sname_ext '.' S.loadext];
            if ~exist(fullfile(S.setpath,sname_ext),'dir')
                mkdir(fullfile(S.setpath,sname_ext));
            end
            EEG = pop_saveset(EEG,'filename',sname,'filepath',fullfile(S.setpath,sname_ext));
            
        else
            fprintf('%s already exists, skipping to next file...\n',checker)
        end
        
    end
end

%% Quick ICA script

S.ICA = 1; % 1 = run ICA, 0 = skip ICA
sname_ext = 'manrej'; % use data from this folder

tic
if S.ICA
    % GET FILE LIST
    S.filepath = fullfile(S.setpath,sname_ext);
    S = getfilelist(S,sname_ext);
    
    loadpath = fullfile(S.setpath,sname_ext);
    for f = 1:length(S.filelist)
        file = S.filelist{f};
        % check if file exists already
        [~,nme,~] = fileparts(file);
        checker = fullfile(S.setpath,'ICA',sprintf('%s*',nme));
        
        if isempty(dir(checker)) || S.overwriteICA % skip existing files, depending on settings
            
            
            EEG = pop_loadset('filename',file,'filepath',loadpath);
            
            EEG = pop_runica(EEG,'extended',1,'interupt','on');
            
            % SAVE
            [pth nme ext] = fileparts(file);
            sname_ext = 'ICA';
            sname = [nme '_' sname_ext '.' S.loadext];
            if ~exist(fullfile(S.setpath,sname_ext),'dir')
                mkdir(fullfile(S.setpath,sname_ext));
            end
            EEG = pop_saveset(EEG,'filename',sname,'filepath',fullfile(S.setpath,sname_ext));
        else
            fprintf('%s already exists, skipping to next file...\n',checker)
        end
    end
end
toc
