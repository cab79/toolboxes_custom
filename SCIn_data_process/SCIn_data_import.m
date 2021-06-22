function [S,D] = SCIn_data_import(S)
% takes outputs from SCIn and puts them in a structure 'D', with field name
% being the S.loadprefixes

% get filename using only first input type (filename prefix)
S.load.prefix = S.load.prefixes(1);

% GET FILE LIST
S.path.file = S.path.raw;
S = getfilelist(S);

% change to the input directory
eval(sprintf('%s', ['cd(''' S.path.file ''')']));

% report if there are no such files
if isempty(S.filelist)
    error('Cannot import! See message above \n');
end

if ~isfield(S.load,'combine_blocks')
    S.load.combine_blocks = 0;
end

% run though all files in a loop and get data
for s = 1:length(S.select.subjects)
    
    % create output D
    D(s).subname = S.select.subjects{s};
        
    % FIND THE FILES FOR THIS SUBJECT
    subfiles = S.filelist(find(not(cellfun('isempty', strfind(S.filelist,D(s).subname)))));
    dirnames = S.dirlist(find(not(cellfun('isempty', strfind(S.filelist,D(s).subname)))));

    % CYCLE THROUGH EACH FILE FOR THIS SUBJECT
    for f = 1:length(subfiles)
        filename = subfiles{f};
        dirname = dirnames{f};
        S.load.prefix = S.load.prefixes{1};
        fprintf('\nProcessing %s.\n\n', filename);

        % assume only one field in loaded file
        temp(f) = load(fullfile(dirname,filename));
        fename = fieldnames(temp(f));
        if isstruct(temp(f).(fename{1}))
            outstruct = temp(f).(fename{1});
        else
            outstruct.Output = temp(f).(fename{1});
        end
        if S.load.combine_blocks && f>1
            fdnames = fieldnames(outstruct);
            for fd = 1:length(fdnames)
                if iscell(outstruct.(fdnames{fd})) || isnumeric(outstruct.(fdnames{fd}))
                    if strcmp(fdnames{fd},'stimtime')
                        % GET NUM TRIALS IN PREV BLOCK
                        prevtrials = length(D(s).(S.load.prefix).stimtime);
                    end
                    if strcmp(fdnames{fd},'presstrial')
                        %prevtrials = length(D(s).(S.load.prefix).stimtime);
                        dat = outstruct.(fdnames{fd})+prevtrials;
                    else
                        dat = outstruct.(fdnames{fd});
                    end
                    D(s).(S.load.prefix).(fdnames{fd}) = [D(s).(S.load.prefix).(fdnames{fd}) dat];
                end
            end
        else
            D(s).(S.load.prefix)(f) = outstruct;
        end

        % get other filetypes with same name
        for fn = 1:length(S.load.prefixes)-1
            S.load.prefix = S.load.prefixes{fn+1};
            filename = strrep(filename,S.load.prefixes{1},S.load.prefix);

            % assume only one field in loaded file
            temp2(f) = load(fullfile(dirname,filename));
            fename = fieldnames(temp2(f));
            if S.load.combine_blocks && f>1
                fdnames = fieldnames(temp2(f).(fename{1}));
                for fd = 1:length(fdnames)
                    if iscell(temp2(f).(fename{1}).(fdnames{fd})) || isnumeric(temp2(f).(fename{1}).(fdnames{fd}))
                        if strcmp(fdnames{fd},'blocks')
                            prevblock = max(D(s).(S.load.prefix).(fdnames{fd}));
                            dat = temp2(f).(fename{1}).(fdnames{fd})+prevblock;
                        elseif strcmp(fdnames{fd},'condnum')
                            prevcond = max(D(s).(S.load.prefix).(fdnames{fd}));
                            dat = temp2(f).(fename{1}).(fdnames{fd})+prevcond;
                        else
                            dat = temp2(f).(fename{1}).(fdnames{fd});
                        end
                        D(s).(S.load.prefix).(fdnames{fd}) = [D(s).(S.load.prefix).(fdnames{fd}) dat];
                    end
                end
            else
                D(s).(S.load.prefix)(f) = temp2(f).(fename{1});
            end
        end
    end
end