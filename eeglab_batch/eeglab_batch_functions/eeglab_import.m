function S=eeglab_import(S)

S.func = 'import';

% previously imported files
if isfield(S.(S.func),'filelist')
    if S.(S.func).overwrite==0
        prev_filelist = S.(S.func).filelist;
        prev_dirlist = S.(S.func).dirlist;
        prev_subj_pdat_idx = S.(S.func).subj_pdat_idx;
        prev_designtab = S.(S.func).designtab;
    end
    S.(S.func) = rmfield(S.(S.func),'dirlist');
    S.(S.func) = rmfield(S.(S.func),'subj_pdat_idx');
    S.(S.func) = rmfield(S.(S.func),'designtab');
end

% GET FILE LIST
S.path.file = S.path.raw;
S = getfilelist(S);

% select which to import
if S.(S.func).overwrite==0 && exist('prev_filelist','var')
    idx_remove = ismember(S.(S.func).filelist,prev_filelist);
    S.(S.func).filelist(idx_remove)=[];
    S.(S.func).dirlist(idx_remove)=[];
    S.(S.func).subj_pdat_idx(idx_remove)=[];
    S.(S.func).designtab(idx_remove,:)=[];
end

% change to the input directory
eval(sprintf('%s', ['cd(''' S.path.raw ''')']));

% report if there are no such files
if isempty(S.(S.func).filelist)
    error('No files found to import!\n');
end

% indices of S.filelist for each subject
col_ind = find(ismember(S.(S.func).designtab(1,:),{'subjects'}));
designtab = cell2table(S.(S.func).designtab(2:end,col_ind)); % convert to table because unique with rows does not work on cell arrays!
[~,first_ind,file_ind]=unique(designtab,'rows','stable');
uni_ind = unique(file_ind);

for sub = 1:length(uni_ind)
    
    % FIND THE FILES
    subfiles = S.(S.func).filelist(file_ind==uni_ind(sub));
    subdirs = S.(S.func).dirlist(file_ind==uni_ind(sub));

    % run though subject files in a loop and convert
    for f = 1:length(subfiles)
        filename = subfiles{f};
        dirname = subdirs{f};
        %if length(file) > 1
        %    error('Expected 1 recording file. Found %d.\n',length(file));
        %else
            %filename = file.name;
            fprintf('\nProcessing %s.\n\n', filename);
        %end

        switch S.(S.func).fname.ext{:}
            case 'vhdr' % Brainvision
                if isfield(S.import,'module') && any(strcmp(S.import.module,'fileio'))
                    EEG = pop_fileio(fullfile(dirname,filename));
                else
                    EEG = pop_loadbv_CAB(dirname,filename); % CAB version changes filesnames to be the same as those on disk
                end
            case 'cnt' % Neuroscan
                EEG = pop_loadcnt(fullfile(dirname,filename));
            case {'mff',''}
                % requires toolbox: mffimport2.0 (Srivas Chennu)
                EEG = pop_readegimff(fullfile(dirname,filename));
            case 'bdf'
                EEG = pop_biosig(fullfile(dirname,filename),[],'BDF'); % NOT TESTED. May need to add in channel range instead of []
            case 'csv'
                % NEEDS CODE
        end
        EEG = eeg_checkset(EEG);


        %add channel locations
        if S.(S.func).chan.addloc
            EEG=pop_chanedit(EEG, 'lookup',S.path.locfile);
        end


        % select chans
        S.(S.func).chanlocs=EEG.chanlocs;
        S=select_chans(S);
        fprintf('Removing excluded channels.\n');
        EEG = pop_select(EEG,'channel',S.(S.func).inclchan);
        S.(S.func).chanlocs = S.(S.func).chanlocs(S.(S.func).inclchan);

        [pth nme ext] = fileparts(filename); 
        % update to ideal file name (currently only works for block names)
        if isfield(S.import.select,'blocks_changenames') && ~isempty(S.import.select.blocks_changenames)
            idx_change = find(~cellfun(@isempty,regexp(nme, regexptranslate('wildcard', S.import.select.blocks))));
            [s1, s2] = regexp(nme, regexptranslate('wildcard', S.import.select.blocks{idx_change}));
            if ~strcmp(nme(s1:s2),S.import.select.blocks_changenames{idx_change})
                nme = strrep(nme,nme(s1:s2),S.import.select.blocks_changenames{idx_change});
            end
        end
        savename = nme;
        EEG.setname = sprintf([S.(S.func).save.prefix{:} '%s' S.(S.func).save.suffix{:}],savename); % the output file is called: basename_orig
        EEG.filename = sprintf([S.(S.func).save.prefix{:} '%s' S.(S.func).save.suffix{:} '.set'],savename);
        EEG.filepath = fullfile(S.path.prep,'imported');
        
        if ~exist(EEG.filepath,'dir')
            mkdir(EEG.filepath)
        end

        fprintf('Saving %s%s.\n', EEG.filepath, EEG.filename);
        pop_saveset(EEG,'filename', EEG.filename, 'filepath', EEG.filepath);
    end
end

if exist('prev_filelist','var')
    S.(S.func).filelist = [S.(S.func).filelist prev_filelist];
    S.(S.func).dirlist = [S.(S.func).dirlist prev_dirlist];
    S.(S.func).subj_pdat_idx = [S.(S.func).subj_pdat_idx prev_subj_pdat_idx];
    S.(S.func).designtab = [S.(S.func).designtab; prev_designtab];
end


