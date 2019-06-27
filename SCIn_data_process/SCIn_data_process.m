function [S,D] = SCIn_data_process(S,D)
% requires strcuture D containing all the data

if ~isfield(S,'trialmax')
    S.trialmax = {1000};
end
if ~isfield(S,'movingavg')
    S.movingavg = 20;
end
if ~isfield(S,'response_to_analyse')
    S.response_to_analyse = 'first'; % by default, analyse first response per trial only
end
if ~isfield(S,'fitsim')
    S.fitsim = 1;
end
if ~isfield(S,'meansim')
    S.meansim = 0; 
end

for d = 1:length(D)
    
    for op = 1:length(D(d).Output)

        % get settings if not in D already
%         if ~isfield(D(d).Output(op),'Settings')
%             A=load(fullfile('C:\Data\Matlab\Matlab_files\NTIP\SCIn\Sequences', D(d).Output(op).SeqName));
%             D(d).Output(op).Settings = A.settings;
%         end
        

        % specify stimtypes if needed
        if isfield(S.accuracy,'cond_stimtype') && ~isempty(S.accuracy.cond_stimtype)
            incl_stimtypes = ismember(D(d).Sequence.signal(S.signal.target,:),S.accuracy.cond_stimtype);
        else
            incl_stimtypes = ones(1,size(D(d).Sequence.signal,2));
        end
        
        
        % create condnum cell which either contains all trials or a
        % subset of trials according to S.trialmax
                % This allows testing of how many trials are needed to
                % produce robust condition differences.
        condnum={};
        if ~isempty(S.select.condtype)
            condind = D(d).Sequence.(S.select.condtype).condnum;
        else
            condind = D(d).Sequence.condnum;
        end
        conds = unique(condind);
        if ~isempty(S.trialmax)
            for tm = 1:length(S.trialmax)
                trialmax = S.trialmax{tm};
                % reduce the number of trials per condition-pair to S.trialmax
                for i = 1:ceil(length(conds)/2)
                    ind = (i-1)*2+1:(i-1)*2+2; % condition pair indices, e.g. 1 2, 3 4, 5 6
                    trialind = find(ismember(condind,conds(ind(ismember(ind,conds)))));
                    trialind_new = trialind(1:min(trialmax,length(trialind)));
                    condnum_orig{tm} = condind;
                    condnum{tm}(trialind) = nan;
                    condnum{tm}(trialind_new) = condnum_orig{tm}(trialind_new);
                end
            end
        else
            condnum{1} = condind;
        end
        
        % identify and index repeated responses on each trial, if they
        % exist, and remove if needed
        D(d).Processed(op).all_resp_ind = 1:length(D(d).Output.presstrial);
        [~,D(d).Processed(op).first_resp_ind] = unique(D(d).Output.presstrial);
        D(d).Processed(op).repeated_response_ind = find(diff(D(d).Output.presstrial)==0);
        D(d).Processed(op).last_resp_ind = D(d).Processed(op).all_resp_ind; D(d).Processed(op).last_resp_ind(D(d).Processed(op).repeated_response_ind)=[];
        % get condition values in each case
        D(d).Processed(op).first_resp_conds = D(d).Sequence.condnum(:,D(d).Output.presstrial(D(d).Processed(op).first_resp_ind));
        D(d).Processed(op).last_resp_conds = D(d).Sequence.condnum(:,D(d).Output.presstrial(D(d).Processed(op).last_resp_ind));
        D(d).Processed(op).repeated_response_conds = D(d).Sequence.condnum(:,D(d).Output.presstrial(D(d).Processed(op).repeated_response_ind));
        % choose which to analyse
        switch S.response_to_analyse
            case 'first'
                resp_index = D(d).Processed(op).first_resp_ind;
            case 'last'
                resp_index = D(d).Processed(op).last_resp_ind;
            case 'all'
                resp_index = D(d).Processed(op).all_resp_ind;
        end
        
         % process RTs first so that trials rejected can also be rejected
         % from choice data
         if S.RT.on
            
            % log response times
            D(d).Processed(op).logrt=nan(1,size(D(d).Sequence.signal,2));
            RT=D(d).Output.RT(resp_index);
            RT(RT<S.RT.min | RT>S.RT.max)=nan; % don't consider RTs less than e.g. 500ms and greater than e.g. 2s
            logRT=log(RT);
            logRT_SD_cutoffs = [nanmedian(logRT)-nanstd(logRT)*S.RT.SDreject, nanmedian(logRT)+nanstd(logRT)*S.RT.SDreject];
            logRT(logRT<logRT_SD_cutoffs(1) | logRT>logRT_SD_cutoffs(2))=nan;
            D(d).Processed(op).logrt(D(d).Output.presstrial(resp_index)) = logRT;
            
            % moving average RT over time
            ma=S.movingavg;
            for i = 1:length(D(d).Processed(op).logrt)
                D(d).Processed(op).malogrt(i) = nansum(D(d).Processed(op).logrt(max(1,i-ma+1):i))/length(max(1,i-ma+1):i);
            end
            
            % create condnum cell which either contains all trials or a
            % subset of trials according to S.trialmax
                    % This allows testing of how many trials are needed to
                    % produce robust condition differences.
%             condnum={};
%             conds = unique(condind);
%             if ~isempty(S.trialmax)
%                 for tm = 1:length(S.trialmax)
%                     trialmax = S.trialmax{tm};
%                     % reduce the number of trials per condition-pair to S.trialmax
%                     for i = 1:ceil(length(conds)/2)
%                         ind = (i-1)*2+1:(i-1)*2+2; % condition pair indices, e.g. 1 2, 3 4, 5 6
%                         trialind = find(ismember(condind,conds(ind(ismember(ind,conds)))));
%                         trialind_new = trialind(1:min(trialmax,length(trialind)));
%                         condnum_orig{tm} = condind;
%                         condnum{tm}(trialind) = nan;
%                         condnum{tm}(trialind_new) = condnum_orig{tm}(trialind_new);
%                     end
%                 end
%             else
%                 condnum{1} = condind;
%             end

            % for each version of condnum (diff number of trials in each),
            % get perc correct for each condition and block
            for cn = 1:length(condnum)
                
                % split into conditions
                condsuni = unique(condnum{cn});
                condsuni = condsuni(~isnan(condsuni));
                D(d).Processed(op).cond_logrt{cn} = nan(1,length(conds));
                for i = 1:length(condsuni)
                    D(d).Processed(op).cond_logrt{cn}(condsuni(i)) = nanmean(D(d).Processed(op).logrt(condnum{cn}==condsuni(i) & incl_stimtypes)); 
                end

                % split into stim intensity per condition
                if S.accuracy.cond_stimtype
                    stims = unique(D(d).Sequence.signal(S.signal.target,:));
                    for i = 1:length(condsuni)
                        for s = 1:length(stims)
                            D(d).Processed(op).stimcond_logrt{cn}{condsuni(i)}(s) = nanmean(D(d).Processed(op).logrt(condnum{cn}==condsuni(i) & D(d).Sequence.signal(S.signal.target,:)==stims(s))); 
                        end
                    end
                end
            end
        end
        
        if S.accuracy.on
            
            S.accuracy.buttons = D(d).Output(1).Settings.buttonopt;
            
            % create presssignal 
            D(d).Processed(op).presssignal = nan(1,size(D(d).Sequence.signal,2));
            presssignal = [];
            for i = 1:length(S.accuracy.buttons)
                % translate actual buttons to their "signal" meaning
                presssignal(strcmp(D(d).Output(op).pressbutton(resp_index),S.accuracy.buttons{i}))=S.accuracy.signal(i);
            end
            % exclude bad trials based on RT limits
            if exist('logRT','var')
                presssignal(isnan(logRT)) = nan;
            end
            D(d).Processed(op).presssignal(D(d).Output(op).presstrial(resp_index)) = presssignal;
            
            % for each target, re-label in terms of correct response
            correct_resp = D(d).Sequence.signal(S.signal.target,:);
            for i = 1:length(S.accuracy.target_resp{1})
                correct_resp(D(d).Sequence.signal(S.signal.target,:)==S.accuracy.target_resp{1}(i)) = S.accuracy.target_resp{2}(i);
            end
            D(d).Processed(op).correct_resp = correct_resp;
            
            % what are the correct responses?
            if isfield(D(d).Output(1).Settings,'adaptive') && isfield(D(d).Output(1).Settings.adaptive,'response_type') && strcmp(D(d).Output(1).Settings.adaptive.response_type,'samediff')
                if numel(S.signal.target)~=2
                    error('Same-Different responses require two targets')
                end
                D(d).Processed(op).correct_resp = abs(diff(correct_resp))+1;
            elseif isfield(D(d).Output(1).Settings,'adaptive') && isfield(D(d).Output(1).Settings.adaptive,'response_type') && strcmp(D(d).Output(1).Settings.adaptive.response_type,'2AFC')
                if numel(S.signal.target)~=2
                    error('2AFC responses require two targets and definition of S.signal.correct_resp')
                end
                [~,D(d).Processed(op).correct_resp] = max(correct_resp==S.signal.correct_resp, [],1);
            else
                % assume signal is the target
                D(d).Processed(op).correct_resp = correct_resp;
            end
            
            
            % compare actual responses to correct responses - which were
            % correct?
            D(d).Processed(op).correct = double(D(d).Processed(op).correct_resp==D(d).Processed(op).presssignal);
            D(d).Processed(op).correct(isnan(D(d).Processed(op).presssignal)) = nan;
            
            % moving average percent correct over time
            ma=S.movingavg;
            for i = 1:length(D(d).Processed(op).correct)
                D(d).Processed(op).macorrect(i) = 100*sum(D(d).Processed(op).correct(max(1,i-ma+1):i))/length(max(1,i-ma+1):i);
            end
            
            % for each version of condnum (diff number of trials in each),
            % get perc correct for each condition and block
            for cn = 1:length(condnum)
                
                % split into conditions
                condsuni = unique(condnum{cn});
                condsuni = condsuni(~isnan(condsuni));
                D(d).Processed(op).condcorrectfract{cn} = nan(1,length(conds));
                
                for i = 1:length(condsuni)
                    D(d).Processed(op).condcorrect{cn}{condsuni(i)} = D(d).Processed(op).correct(condnum{cn}==condsuni(i) & incl_stimtypes); 
                    D(d).Processed(op).numtrials{cn}{condsuni(i)} = sum(~isnan(D(d).Processed(op).condcorrect{cn}{condsuni(i)}));
                    D(d).Processed(op).condcorrectfract{cn}(condsuni(i)) = nansum(D(d).Processed(op).condcorrect{cn}{condsuni(i)})/D(d).Processed(op).numtrials{cn}{condsuni(i)};
                end

                % split into blocks
                blocks = unique(D(d).Sequence.blocks);
                for i = 1:length(condsuni)
                    for b = 1:length(blocks)
                        D(d).Processed(op).blockcondcorrect{cn}{condsuni(i)}{b} = D(d).Processed(op).correct(condnum{cn}==condsuni(i) & D(d).Sequence.blocks==blocks(b) & ismember(D(d).Sequence.signal(S.signal.target,:),S.accuracy.cond_stimtype)); 
                        D(d).Processed(op).blocknumtrials{cn}{condsuni(i)}{b} = sum(~isnan(D(d).Processed(op).blockcondcorrect{cn}{condsuni(i)}{b}));
                        D(d).Processed(op).blockcondcorrectfract{cn}{condsuni(i)}(b) = nansum(D(d).Processed(op).blockcondcorrect{cn}{condsuni(i)}{b})/D(d).Processed(op).blocknumtrials{cn}{condsuni(i)}{b};
                        for ii = 1:length(D(d).Processed(op).blockcondcorrect{cn}{condsuni(i)}{b})
                            D(d).Processed(op).blockcondcorrectmovavg{cn}{condsuni(i)}{b}(ii) = 100*sum(D(d).Processed(op).blockcondcorrect{cn}{condsuni(i)}{b}(max(1,ii-ma+1):ii))/length(max(1,ii-ma+1):ii);;
                        end
                    end
                end
                
                
                % split into blocks before conditions, to get moving average of condition ratio of
                % percent correct
                blocks = unique(D(d).Sequence.blocks);
                for b = 1:length(blocks)
                    trialidx=find(D(d).Sequence.blocks==blocks(b) & ismember(D(d).Sequence.signal(S.signal.target,:),S.accuracy.cond_stimtype));
                    D(d).Processed(op).blockcorrect{cn}{b} = D(d).Processed(op).correct(trialidx);
                    cond = condnum{cn}(trialidx);
                    % for each trial, get indices of last ma trials
                    for ii = 1:length(trialidx)
                        idx=max(1,ii-ma+1):ii;
                        correcttrials=D(d).Processed(op).blockcorrect{cn}{b}(idx);
                        condidx=cond(idx);
                        ucond = unique(condidx);
                        correct=[];
                        Ntrials=[];
                        for ci =1:length(ucond)
                            % proportion correct, from 0 to 1: correctN / totalN
                            Ntrials(ci)=length(idx(condidx==ucond(ci)));
                            correct(ci) = nansum(D(d).Processed(op).blockcorrect{cn}{b}(idx(condidx==ucond(ci))))/Ntrials(ci);
                        end
                        if isempty(correct) || length(correct)<2 || any(Ntrials<2)
                            D(d).Processed(op).blockcorrectmovavg{cn}{b}(ii) = NaN;
                        else
                            % difference in proportion correct ranging from -1 to 1
                            D(d).Processed(op).blockcorrectmovavg{cn}{b}(ii) = correct(1)-correct(2);
                        end
                    end
                end
                
                % split into stim intensity per condition
                stims = unique(D(d).Sequence.signal(S.signal.target,:));
                if S.accuracy.cond_stimtype
                    for i = 1:length(condsuni)
                        for s = 1:length(stims)
                            D(d).Processed(op).stimcondcorrect{cn}{condsuni(i)}{s} = D(d).Processed(op).correct(condnum{cn}==condsuni(i) & D(d).Sequence.signal(S.signal.target,:)==stims(s)); 
                            D(d).Processed(op).stimnumtrials{cn}{condsuni(i)}{s} = sum(~isnan(D(d).Processed(op).stimcondcorrect{cn}{condsuni(i)}{s}));
                            D(d).Processed(op).stimcondcorrectfract{cn}{condsuni(i)}(s) = nansum(D(d).Processed(op).stimcondcorrect{cn}{condsuni(i)}{s})/D(d).Processed(op).stimnumtrials{cn}{condsuni(i)}{s};
                        end
                    end
                end
                
                % split into stim intensity per cue
                if S.signal.cue
                    cues = unique(D(d).Sequence.signal(S.signal.cue,:));
                    for i = 1:length(cues)
                        for s = 1:length(stims)
                            D(d).Processed(op).stimcuecorrect{cn}{cues(i)}{s} = D(d).Processed(op).correct(D(d).Sequence.signal(S.signal.cue,:)==cues(i) & D(d).Sequence.signal(S.signal.target,:)==stims(s)); 
                            D(d).Processed(op).stimcuenumtrials{cn}{cues(i)}{s} = sum(~isnan(D(d).Processed(op).stimcuecorrect{cn}{cues(i)}{s}));
                            D(d).Processed(op).stimcuecorrectfract{cn}{cues(i)}(s) = nansum(D(d).Processed(op).stimcuecorrect{cn}{cues(i)}{s})/D(d).Processed(op).stimcuenumtrials{cn}{cues(i)}{s};
                        end
                    end
                end
            end
        end
       
    end
    
    % average
    if S.meansim && S.fitsim==2
        S.mean_fields = {};
        if S.accuracy.on
            S.mean_fields = [S.mean_fields {'macorrect','blockcorrectmovavg','condcorrectfract','blockcondcorrectfract','stimcondcorrectfract','stimcuecorrectfract'}];
        end
        if S.RT.on
            S.mean_fields = [S.mean_fields {'logrt','cond_logrt','stimcond_logrt'}];
        end
        for mf = 1:length(S.mean_fields)
            if iscell(D(d).Processed(1).(S.mean_fields{mf}))
                for i1 = 1:length(D(d).Processed(1).(S.mean_fields{mf}))
                    if iscell(D(d).Processed(1).(S.mean_fields{mf}){i1})
                        for i2 = 1:length(D(d).Processed(1).(S.mean_fields{mf}){i1})
                            temp_mf=[];
                            for op = 1:length(D(d).Processed)
                                temp_mf(op,:) = D(d).Processed(op).(S.mean_fields{mf}){i1}{i2};
                            end
                            D(d).Processed(1).(S.mean_fields{mf}){i1}{i2} = mean(temp_mf,1);
                        end
                    else
                        temp_mf=[];
                        for op = 1:length(D(d).Processed)
                            temp_mf(op,:) = D(d).Processed(op).(S.mean_fields{mf}){i1};
                        end
                        D(d).Processed(1).(S.mean_fields{mf}){i1} = mean(temp_mf,1);
                    end
                end
            else
                temp_mf=[];
                for op = 1:length(D(d).Processed)
                    temp_mf(op,:) = D(d).Processed(op).(S.mean_fields{mf});
                end
                D(d).Processed(1).(S.mean_fields{mf}) = mean(temp_mf,1);
            end
        end
        D(d).Processed = D(d).Processed(1);
    end
    
    % save individual tables for each subject's processed data trials
    if S.save.tables
        if ~isempty(op)
            trial_data_table.trial_number = D(d).Output(1).presstrial(resp_index)';
            trial_data_table.block_number = D(d).Sequence.blocks(D(d).Output(1).presstrial(resp_index))';
            trial_data_table.condition_number = condind(D(d).Output(1).presstrial(resp_index))';
            trial_data_table.response = D(d).Output(1).pressbutton(resp_index)';
            trial_data_table.correct_response = D(d).Processed(1).correct_resp(D(d).Output(1).presstrial(resp_index))';
            trial_data_table.correct = D(d).Processed(1).correct(D(d).Output(1).presstrial(resp_index))';
            trial_data_table.accuracy_moving_average = D(d).Processed(1).macorrect(D(d).Output(1).presstrial(resp_index))';
            trial_data_table.log_response_time = D(d).Processed(1).logrt(D(d).Output(1).presstrial(resp_index))';
            trial_data_table.log_response_time_moving_average = D(d).Processed(1).malogrt(D(d).Output(1).presstrial(resp_index))';
            writetable(struct2table(trial_data_table),fullfile(S.path.prep,[S.expt '_' D(d).subname '.xlsx']));

        end
    end
    
end

% save a single table for all subjects' condition-mean data
if S.save.tables
    % add data to group condition table
    for d = 1:length(D)
        fraction_correct_table.subject{d,1} = D(d).subname;
        log_response_time_table.subject{d,1} = D(d).subname;
        for c = 1:length(conds)
            try
                fraction_correct_table.(['condition_' num2str(conds(c))])(d,1) = D(d).Processed(1).condcorrectfract{1}(c);
                log_response_time_table.(['condition_' num2str(conds(c))])(d,1) = D(d).Processed(1).cond_logrt{1}(c);
            catch
                fraction_correct_table.(['condition_' num2str(conds(c))])(d,1) = nan;
                log_response_time_table.(['condition_' num2str(conds(c))])(d,1) = nan;
            end
        end
    end
    writetable(struct2table(fraction_correct_table),fullfile(S.path.prep,[S.expt '_fraction_correct_' datestr(now,30) '.xlsx']));
    writetable(struct2table(log_response_time_table),fullfile(S.path.prep,[S.expt '_log_response_time_' datestr(now,30) '.xlsx']));
end