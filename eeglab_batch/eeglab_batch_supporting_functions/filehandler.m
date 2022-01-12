function S = filehandler(S,step)

switch step
    case 'start'

        % create directory for saving if doesn't exist
        if ~exist(fullfile(S.path.(S.func),S.(S.func).save.suffix{:}),'dir')
            mkdir(fullfile(S.path.(S.func),S.(S.func).save.suffix{:}));
        end

        % load directory
        if isfield(S.(S.func),'S.path.file')
            S.path.file = fullfile(S.path.(S.func),S.(S.func).load.suffix{:});
        else
            S.path.file = fullfile(S.path.(S.func),S.(S.func).load.suffix{:});
        end
        
        % get previously processed files
        if isfield(S.(S.func),'designtab')
            if S.(S.func).overwrite==0
                prev_designtab = S.(S.func).designtab;
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
        end
        S.fn = height(S.(S.func).outtable);

    case 'update'
        idx_processed = ismember(S.(S.func).designtab.file,S.(S.func).file_processed);
        S.(S.func).designtab.processed(idx_processed)=1;

        % save S
        save(fullfile(S.path.prep,S.sname),'S'); % saves 'S' - will be overwritten each time

        % save table
        if isfield(S.prep,'outtable_name')
            writetable(S.prep.outtable,fullfile(fullfile(S.path.prep,S.prep.outtable_name)))
        end

%         if exist('prev_filelist','var')
%             S.(S.func).filelist = [S.(S.func).filelist prev_filelist];
%             S.(S.func).dirlist = [S.(S.func).dirlist prev_dirlist];
%             S.(S.func).subj_pdat_idx = [S.(S.func).subj_pdat_idx prev_subj_pdat_idx];
%             S.(S.func).designtab = [S.(S.func).designtab; prev_designtab];
%         end

end
