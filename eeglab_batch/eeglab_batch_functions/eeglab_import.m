function S=eeglab_import(S)

S.func = 'import';
S = filehandler(S,'start');

% change to the input directory
eval(sprintf('%s', ['cd(''' S.path.import ''')']));

% report if there are no such files
if isempty(S.(S.func).filelist)
    error('No files found to import!\n');
end

% indices of S.filelist for each subject
% col_ind = find(ismember(S.(S.func).designtab(1,:),{'subjects'}));
% designtab = cell2table(S.(S.func).designtab(2:end,col_ind)); % convert to table because unique with rows does not work on cell arrays!
% [~,first_ind,file_ind]=unique(S.(S.func).designtab,'rows','stable');
% uni_ind = unique(file_ind);

for f = 1:length(S.(S.func).filelist)
    
%     % FIND THE FILES
%     subfiles = S.(S.func).filelist(sub);
%     subdirs = S.(S.func).designtab.dir(S.fn+sub);

    % run though subject files in a loop and convert
%     for f = 1:length(subfiles)
        filename = S.(S.func).filelist{f};
        dirname = S.(S.func).designtab.dir{ismember(S.(S.func).designtab.file,filename)};
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
                    % get impedances
                    imptable = import_impedance_vhdr(fullfile(dirname,filename),S);
                    S.(S.func).outtable.file{S.fn+f} = filename;
                    S.(S.func).outtable(S.fn+f,2:end)=array2table(NaN(1,width(S.(S.func).outtable)-1));
                    for ch=1:height(imptable)
                        S.(S.func).outtable.(imptable.chan(ch))(S.fn+f) = imptable.imp(ch);
                    end
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

        S.(S.func).file_processed=filename;
        S = filehandler(S,'update');
%     end
end

function outtable = import_impedance_vhdr(filename, S)

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = ["     ", ":"];

% Specify column names and types
opts.VariableNames = ["VarName1", "VarName2", "Var3"];
opts.SelectedVariableNames = ["VarName1", "VarName2"];
opts.VariableTypes = ["string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";

% Specify variable properties
opts = setvaropts(opts, ["VarName1", "Var3"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["VarName1", "Var3"], "EmptyFieldRule", "auto");
% opts = setvaropts(opts, "VarName2", "ThousandsSeparator", ",");

% Import the data
dat = readtable(filename, opts);

imline = find(contains(dat.VarName1,"Impedance [kOhm]"))+1;

outtable = table;
outtable.chan = dat.VarName1(imline:end);
outtable.imp = dat.VarName2(imline:end);

if isfield(S.import,'Ref')
    outtable(ismember(outtable.chan,"Ref"),:)=[];
    outtable.chan(ismember(outtable.chan,S.import.Ref))="Ref";
end
