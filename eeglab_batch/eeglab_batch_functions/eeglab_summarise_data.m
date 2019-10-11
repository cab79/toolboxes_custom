function S=eeglab_summarise_data(S,part)
S.func = 'summ';

S.summ.out = table;

switch part
    case 'chan'
        
        % load directory
        S.prep.loaddir = fullfile(S.path.prep,S.(S.func).load.suffix{:});

        % GET FILE LIST
        S.path.file = S.prep.loaddir;
        S = getfilelist(S,S.(S.func).load.suffix);

        % change to the input directory
        eval(sprintf('%s', ['cd(''' S.path.prep ''')']));

        % report if there are no such files
        if isempty(S.(S.func).filelist)
            error('No files found!\n');
        end

        if ~any(strcmp('subject',S.summ.out.Properties.VariableNames))
            S.summ.out.subject = S.(S.func).select.subjects';
        end

        for s = 1:length(S.(S.func).select.subjects)
            subfiles = S.(S.func).filelist(find(not(cellfun('isempty', strfind(S.(S.func).filelist,S.(S.func).select.subjects{s})))));

            for f = 1:length(subfiles)
                filename = subfiles{f};

                fprintf('\nProcessing %s.\n\n', filename);

                % GET FILENAME PARTS
                [pth nme ext] = fileparts(filename); 
                fparts = strsplit(nme,'_');

                % LOAD DATA
                EEG = pop_loadset('filename',filename,'filepath',S.path.file);

                % summarise
                if isfield(EEG.chanlocs,'interp')
                    interp(s,f) = sum(cell2mat({EEG.chanlocs.interp}));
                else
                    interp(s,f) = nan;
                end
            end
        end
        for f = 1:size(interp,2)
            S.summ.out.(['ninterp' num2str(f)]) = interp(:,f);
        end
        
        % plots
        drawboxplot(S,'ninterp')
    
    case 'IC'
        
        % load directory
        S.prep.loaddir = fullfile(S.path.prep,S.(S.func).load.suffix{:});

        % GET FILE LIST
        S.path.file = S.prep.loaddir;
        S = getfilelist(S,S.(S.func).load.suffix);

        % change to the input directory
        eval(sprintf('%s', ['cd(''' S.path.prep ''')']));

        % report if there are no such files
        if isempty(S.(S.func).filelist)
            error('No files found!\n');
        end

        if ~any(strcmp('subject',S.summ.out.Properties.VariableNames))
            S.summ.out.subject = S.(S.func).select.subjects';
        end

        % run though all files in a loop
        for s = 1:length(S.(S.func).select.subjects)
            subfiles = S.(S.func).filelist(find(not(cellfun('isempty', strfind(S.(S.func).filelist,S.(S.func).select.subjects{s})))));

            for f = 1:length(subfiles)
                filename = subfiles{f};

                fprintf('\nProcessing %s.\n\n', filename);

                % GET FILENAME PARTS
                [pth nme ext] = fileparts(filename); 
                fparts = strsplit(nme,'_');

                % LOAD DATA
                EEG = pop_loadset('filename',filename,'filepath',S.path.file);

                % get components
                comp=EEG.reject.gcompreject;
                ncomp(s,f) = length(comp);
                ncomp_rej(s,f) = sum(comp);
                frac_comp_rej(s,f) = sum(comp)/length(comp);
            end

        end
        for f = 1:size(ncomp,2)
            S.summ.out.(['ncomp' num2str(f)]) = ncomp(:,f);
            S.summ.out.(['ncomp_rej' num2str(f)]) = ncomp_rej(:,f);
            S.summ.out.(['frac_comp_rej' num2str(f)]) = frac_comp_rej(:,f);
        end
        
        % plots
        drawboxplot(S,'ncomp')
        drawboxplot(S,'ncomp_rej')
        drawboxplot(S,'frac_comp_rej')

    case 'trial'
        
        % load directory
        S.prep.loaddir = fullfile(S.path.prep,S.(S.func).load.suffix{:});

        % GET FILE LIST
        S.path.file = S.prep.loaddir;
        S = getfilelist(S,'epoched_cleaned');

        % change to the input directory
        eval(sprintf('%s', ['cd(''' S.path.prep ''')']));

        % report if there are no such files
        if isempty(S.(S.func).filelist)
            error('No files found!\n');
        end

        if ~any(strcmp('subject',S.summ.out.Properties.VariableNames))
            S.summ.out.subject = S.(S.func).select.subjects';
        end

        for s = 1:length(S.(S.func).select.subjects)
            if isempty(S.(S.func).select.sessions)
                S.(S.func).select.sessions = {''};
            end
            if isempty(S.(S.func).select.blocks)
                S.(S.func).select.blocks = {''};
            end
            for a = 1:length(S.(S.func).select.sessions)
                if length(S.(S.func).select.sessions)>1
                    sess = ['_sess' num2str(a)];
                else
                    sess = '';
                end

                for b = 1:length(S.(S.func).select.blocks)
                    if length(S.(S.func).select.blocks)>1
                        block = ['_block' S.(S.func).select.blocks{b}];
                    else
                        block = '';
                    end
                    blockname='';
                    if length(S.(S.func).select.blocks{b})==1
                        blockname = ['_' S.(S.func).select.blocks{b} '_']; % add on underscores
                    end

                    % FIND THE FILES FOR THIS SUBJECT
                    if ~isempty(S.(S.func).select.sessions{a})
                        if isempty(blockname)
                            subfiles = S.(S.func).filelist(find(not(cellfun('isempty', strfind(S.(S.func).filelist,S.(S.func).select.subjects{s})))...
                                & not(cellfun('isempty', strfind(S.(S.func).filelist,S.(S.func).select.sessions{a})))));
                        else
                            subfiles = S.(S.func).filelist(find(not(cellfun('isempty', strfind(S.(S.func).filelist,S.(S.func).select.subjects{s})))...
                                & not(cellfun('isempty', strfind(S.(S.func).filelist,S.(S.func).select.sessions{a})))...
                                & not(cellfun('isempty', strfind(S.(S.func).filelist,blockname)))));
                        end
                    else
                        if isempty(blockname)
                            subfiles = S.(S.func).filelist(find(not(cellfun('isempty', strfind(S.(S.func).filelist,S.(S.func).select.subjects{s})))));
                        else
                            subfiles = S.(S.func).filelist(find(not(cellfun('isempty', strfind(S.(S.func).filelist,S.(S.func).select.subjects{s})))...
                                & not(cellfun('isempty', strfind(S.(S.func).filelist,blockname)))));
                        end
                    end

                    % CYCLE THROUGH EACH FILE FOR THIS SUBJECT
                    for f = 1:length(subfiles)
                        
                        if length(subfiles)>1
                            subf = ['_' num2str(f)];
                        else
                            subf = '';
                        end

                        % load 
                        EEG = pop_loadset('filename',subfiles{f},'filepath',S.path.file);

                        % get trials
                        if ~isfield(S.summ,'events') || isempty(S.summ.events)
                            S.summ.events = {''};
                        end
                        for c = 1:length(S.summ.events)
                            if ~isempty(EEG.epoch)
                                for e = 1:length(EEG.epoch)
                                    if iscell(EEG.epoch(e).eventtype)
                                        EEG.epoch(e).eventtype = EEG.epoch(e).eventtype{S.summ.selectevent};
                                    end
                                end
                                cellind = strcmp(S.summ.events{c},{EEG.epoch.eventtype});
                                % ntrials
                                ID = [sess block '_event' S.summ.events{c}(find(~isspace(S.summ.events{c}))) subf];
                                S.summ.out.(['ntrials' ID])(s) = sum(cellind);
                                
                                % SNR
                                base=dsearchn(EEG.times',1000*S.summ.snr_win{1}');
                                sig=dsearchn(EEG.times',1000*S.summ.snr_win{2}');
                                dat{1}=EEG.data(:,base(1):base(2),cellind);
                                dat{2}=EEG.data(:,sig(1):sig(2),cellind);
                                S.summ.out.(['SNR' ID])(s) = mean(rms(std(dat{2},{},1),2),3) / mean(rms(std(dat{1},{},1),2),3); % mean over trials of the RMS over time of the GFP over channels
                            else
                                S.summ.out.(['ntrials' ID])(s) = 0;
                                S.summ.out.(['SNR' ID])(s) = nan;
                            end
                        end

                    end

                end
            end
        end
        
        % plot ntrials
        headerstr = 'ntrials';
        cols = find(~cellfun(@isempty,strfind(S.summ.out.Properties.VariableNames,headerstr)));
        headers = S.summ.out.Properties.VariableNames(cols);
        dat = table2array(S.summ.out(:,cols));
        [~,datorder] = sort(mean(dat,1));
        figure; boxplot(dat); xticklabels(headers); set(gca, 'XTickLabelRotation', 90)
        
        % plots
        drawboxplot(S,'ntrials')
        drawboxplot(S,'SNR')
        
end

function drawboxplot(S,headerstr)
cols = find(~cellfun(@isempty,strfind(S.summ.out.Properties.VariableNames,headerstr)));
headers = S.summ.out.Properties.VariableNames(cols);
dat = table2array(S.summ.out(:,cols));
%[~,datorder] = sort(mean(dat,1));
figure; 
boxplot(dat); 
xticklabels(headers); 
set(gca, 'XTickLabelRotation', 90)
title(headerstr)