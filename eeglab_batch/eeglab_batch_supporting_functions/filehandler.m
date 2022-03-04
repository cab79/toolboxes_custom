function S = filehandler(S,step)

% save directory principle:
% 1. creates directory fullfile(S.path.(S.func),S.(S.func).save.suffix{:})
% - THIS SHOULD BE HANDLED BY THE CALLING FUNCTION INSTEAD?
% 2. by default, save S to S.path.(S.func), unless S.(S.func).save.dir{:}
    % is defined, in which case saves to that subdirectory within S.path.(S.func)
% 3. output table saves to fullfile(S.path.(S.func),S.(S.func).save.dir{:})
    % unless specified by S.(S.func).outtable_dir

% BETTER:
% loading: always specify loaddir and suffix - latter only relates to
% filenames. Can specific suffix first so it can be used in settings to
% specifiy load.dir.
% saving: always specify save suffix first, then save.dir for data, then
% savedirS for S, and also for table

% FUNCTION:
% creates save.dir for data and for S.

% now requires:
% S.(S.func).load.dir
% S.(S.func).save.dir
% S.(S.func).save.dir_S
% S.(S.func).save.dir_outtable

switch step
    case 'start'

%         % create directory for saving if doesn't exist
%         if isfield(S.(S.func),'save')
%             if ~exist(fullfile(S.path.(S.func),S.(S.func).save.suffix{:}),'dir')
%                 mkdir(fullfile(S.path.(S.func),S.(S.func).save.suffix{:}));
%             end
%         end
        % create directories for saving if doesn't exist
        if isfield(S.(S.func),'save') 
            if isfield(S.(S.func).save,'dir') && ~exist(fullfile(S.(S.func).save.dir),'dir')
                mkdir(fullfile(S.(S.func).save.dir));
            end
            if isfield(S.(S.func).save,'dir_S') && ~exist(fullfile(S.(S.func).save.dir_S),'dir')
                mkdir(fullfile(S.(S.func).save.dir_S));
            end
            if isfield(S.(S.func).save,'dir_outtable') && ~exist(fullfile(S.(S.func).save.dir_outtable),'dir')
                mkdir(fullfile(S.(S.func).save.dir_outtable));
            end
        end

        % load directory
        S.path.file = S.(S.func).load.dir;
        
        % get previously processed files
        if isfield(S.(S.func),'designtab')
            if S.(S.func).overwrite==0
                prev_designtab = S.(S.func).designtab;

                % if subjects have been specifically selected, re-process them even if they
                % were already processed
                if ~isempty(S.(S.func).select.subjects) && ~strcmp(S.(S.func).select.subjects,{''})
                    prev_designtab.processed(ismember(prev_designtab.subjects,S.(S.func).select.subjects))=0;
                end
            end
            S.(S.func) = rmfield(S.(S.func),'designtab');
        end
        
        % GET DESIGNTAB, INCLUDING FILE LIST
        S = getfilelist(S,S.(S.func).load.suffix);
        
        % update designtab with those already processed and generate
        % filelist to process
        if S.(S.func).overwrite==0 && exist('prev_designtab','var')
            idx_processed = ismember(S.(S.func).designtab.file,prev_designtab.file(prev_designtab.processed==1));
            S.(S.func).designtab.processed(idx_processed)=1;
            S.(S.func).filelist = S.(S.func).designtab.file(~idx_processed);
        else
            S.(S.func).filelist = S.(S.func).designtab.file;
        end
        
        %S.loadpath = S.path.file;
        
        % setup output table if needed
        if ~isfield(S.(S.func),'outtable') || S.(S.func).overwrite==1
            S.(S.func).outtable = table;
        else
            % remove entries from table if they will be re-processed
            S.(S.func).outtable(ismember(S.(S.func).outtable.file, S.(S.func).filelist),:) =[];
        end
        
        S.fn = height(S.(S.func).outtable);

    case 'update'
        idx_processed = ismember(S.(S.func).designtab.file,S.(S.func).file_processed);
        S.(S.func).designtab.processed(idx_processed)=1;

        % temporary - can remove this later as it's a patch
        % remove duplicate earlier entries from table if they have been re-processed
        [~,ia,~] = unique(S.(S.func).outtable.file,'last');
        newtable = S.(S.func).outtable;
        S.(S.func).outtable =newtable(ia,:);

        % save S
        save(fullfile(S.(S.func).save.dir_S,S.sname),'S'); % saves 'S' - will be overwritten each time

        % save table
        if isfield(S.(S.func),'outtable_name')
            writetable(S.(S.func).outtable,fullfile(S.(S.func).save.dir_outtable,S.(S.func).outtable_name))
        end

end
