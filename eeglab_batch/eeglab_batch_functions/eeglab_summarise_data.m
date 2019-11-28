function S=eeglab_summarise_data(S,part)
S.func = 'summ';

if ~isfield(S.summ,'out')
    S.summ.out = table;
end

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
        figure; 
        subplot(2,2,1)
        outliers=drawboxplot(S,'ninterp',['number of interpolated channels']);
        S.summ.out.('ninterp_outlier_lower') = double(any(cell2mat(outliers(:,1)'),2));        
        S.summ.out.('ninterp_outlier_upper') = double(any(cell2mat(outliers(:,2)'),2));
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
            S.summ.out.(['ncomp_all' num2str(f)]) = ncomp(:,f);
            S.summ.out.(['ncomp_rej' num2str(f)]) = ncomp_rej(:,f);
            S.summ.out.(['frac_comp_rej' num2str(f)]) = frac_comp_rej(:,f);
        end
        
        % plots
        %figure; 
        subplot(2,2,2); outliers=drawboxplot(S,'ncomp_all',['total number of estimated ICs']);
        S.summ.out.('ncomp_all_outlier_lower') = double(any(cell2mat(outliers(:,1)'),2));        
        S.summ.out.('ncomp_all_outlier_upper') = double(any(cell2mat(outliers(:,2)'),2));
        subplot(2,2,3); outliers=drawboxplot(S,'ncomp_rej',['number of rejected ICA components']);
        S.summ.out.('ncomp_rej_outlier_lower') = double(any(cell2mat(outliers(:,1)'),2));        
        S.summ.out.('ncomp_rej_outlier_upper') = double(any(cell2mat(outliers(:,2)'),2));
        subplot(2,2,4); outliers=drawboxplot(S,'frac_comp_rej',['fraction of rejected ICA components']);
        S.summ.out.('frac_comp_rej_outlier_lower') = double(any(cell2mat(outliers(:,1)'),2));        
        S.summ.out.('frac_comp_rej_outlier_upper') = double(any(cell2mat(outliers(:,2)'),2));
        
    case 'trial'
        
        % load directory
        S.prep.loaddir = fullfile(S.path.prep,S.(S.func).load.suffix{:});

        % GET FILE LIST
        S.path.file = S.prep.loaddir;
        S = getfilelist(S,'epoched*');

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
                                selectchans = find(ismember({EEG.chanlocs(:).labels},S.summ.snr_chan));
                                S.summ.out.(['signal_gfp' ID])(s) = mean(rms(std(dat{2},{},1),2),3); % mean over trials of the RMS over time of the GFP over channels
                                S.summ.out.(['noise_gfp' ID])(s) = mean(rms(std(dat{1},{},1),2),3); % mean over trials of the RMS over time of the GFP over channels
                                S.summ.out.(['SNR_gfp' ID])(s) = mean(rms(std(dat{2},{},1),2),3) / mean(rms(std(dat{1},{},1),2),3); % mean over trials of the RMS over time of the GFP over channels
                                S.summ.out.(['signal_chan' ID])(s) = mean(rms(mean(dat{2}(selectchans,:),1),2),3); % mean over trials of the RMS over time of the GFP over channels
                                S.summ.out.(['noise_chan' ID])(s) = mean(rms(mean(dat{1}(selectchans,:),1),2),3); % mean over trials of the RMS over time of the GFP over channels
                                S.summ.out.(['SNR_chan' ID])(s) = mean(rms(mean(dat{2}(selectchans,:),1),2),3) / mean(rms(mean(dat{1}(selectchans,:),1),2),3); % mean over trials of the RMS over time of the GFP over channels
                            else
                                S.summ.out.(['ntrials' ID])(s) = 0;
                                S.summ.out.(['signal_gfp' ID])(s) = nan;
                                S.summ.out.(['noise_gfp' ID])(s) = nan;
                                S.summ.out.(['SNR_gfp' ID])(s) = nan;
                                S.summ.out.(['signal_chan' ID])(s) = nan;
                                S.summ.out.(['noise_chan' ID])(s) = nan;
                                S.summ.out.(['SNR_chan' ID])(s) = nan;
                            end
                        end

                    end

                end
            end
        end
        
        % mean SNR over conditions
        cols = find(~cellfun(@isempty,strfind(S.summ.out.Properties.VariableNames,'SNR_gfp')));
        dat = table2array(S.summ.out(:,cols));
        S.summ.out.('SNR_gfp_mean') = nanmean(dat,2);
        cols = find(~cellfun(@isempty,strfind(S.summ.out.Properties.VariableNames,'SNR_chan')));
        dat = table2array(S.summ.out(:,cols));
        S.summ.out.('SNR_chan_mean') = nanmean(dat,2);
        
        % plots
        figure; 
        outliers=drawboxplot(S,'ntrials',['total number of trials']);
        S.summ.out.('ntrials_outlier_lower') = double(any(cell2mat(outliers(:,1)'),2));
        S.summ.out.('ntrials_outlier_upper') = double(any(cell2mat(outliers(:,2)'),2));
        figure; 
        subplot(2,3,1); outliers=drawboxplot(S,'signal_gfp',['signal gfp']);
        S.summ.out.('signal_gfp_outlier_lower') = double(any(cell2mat(outliers(:,1)'),2)); 
        S.summ.out.('signal_gfp_outlier_upper') = double(any(cell2mat(outliers(:,2)'),2));
        subplot(2,3,2); outliers=drawboxplot(S,'noise_gfp',['noise gfp']);
        S.summ.out.('noise_gfp_outlier_lower') = double(any(cell2mat(outliers(:,1)'),2));  
        S.summ.out.('noise_gfp_outlier_upper') = double(any(cell2mat(outliers(:,2)'),2));
        subplot(2,3,3); outliers=drawboxplot(S,'SNR_gfp',['SNR gfp']);
        S.summ.out.('SNR_gfp_outlier_lower') = double(any(cell2mat(outliers(:,1)'),2));        
        S.summ.out.('SNR_gfp_outlier_upper') = double(any(cell2mat(outliers(:,2)'),2));
        subplot(2,3,4); outliers=drawboxplot(S,'signal_chan',['signal chan']);
        S.summ.out.('signal_chan_outlier_lower') = double(any(cell2mat(outliers(:,1)'),2)); 
        S.summ.out.('signal_chan_outlier_upper') = double(any(cell2mat(outliers(:,2)'),2));
        subplot(2,3,5); outliers=drawboxplot(S,'noise_chan',['noise chan']);
        S.summ.out.('noise_chan_outlier_lower') = double(any(cell2mat(outliers(:,1)'),2));  
        S.summ.out.('noise_chan_outlier_upper') = double(any(cell2mat(outliers(:,2)'),2));
        subplot(2,3,6); outliers=drawboxplot(S,'SNR_chan',['SNR chan']);
        S.summ.out.('SNR_chan_outlier_lower') = double(any(cell2mat(outliers(:,1)'),2));        
        S.summ.out.('SNR_chan_outlier_upper') = double(any(cell2mat(outliers(:,2)'),2));
        
end

function outliers=drawboxplot(S,headerstr,ttl)

% get data
cols = find(~cellfun(@isempty,strfind(S.summ.out.Properties.VariableNames,headerstr)));
headers = S.summ.out.Properties.VariableNames(cols);
dat = table2array(S.summ.out(:,cols));

% boxplots
hb=boxplot(dat,'Symbol','rx'); 
% hb=boxplot(dat,'PlotStyle','compact','Symbol','rx'); 
if length(headers)>1
    xticklabels(headers); 
    set(gca, 'XTickLabelRotation', 90)
else
    xticklabels({''}); 
end
title(ttl)

% get outliers
hOut = findobj(hb,'Tag','Outliers');
yy = get(hOut,'YData');
med = nanmedian(dat);
for i = 1:length(med)
    if iscell(yy)
        outliers{i,1} = ismember(dat(:,i),unique(yy{i}(yy{i}<med(i))));
        outliers{i,2} = ismember(dat(:,i),unique(yy{i}(yy{i}>med(i))));
    elseif isnan(yy)
        outliers{i,1} = zeros(length(dat),1);
        outliers{i,2} = zeros(length(dat),1);
    else
        outliers{i,1} = ismember(dat(:,i),unique(yy(yy<med(i))));
        outliers{i,2} = ismember(dat(:,i),unique(yy(yy>med(i))));
    end
end