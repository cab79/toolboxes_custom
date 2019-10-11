function [EEG,S] = manualRejectEEG(EEG_in,S,z_thresh)
% This function allows GUI-based manual rejection of epochs.

% TRIAL MARKING: Select any trials to be marked for rejection and click 'update marks'
% TRIAL REJECTION: type 'manrej' into pop-up box to remove marked trials. Or type trial
% numbers to remove them.

% CHANNEL INTERPOLATION: Type channel name into pop-up box

% UNDOING CHANGES: Type 'previous' or 'original' into top answer pop-up box to revert
% changes by one step, or to go back to original EEG set

% FINISHING: to finish and break the loop, type 'manrej' or trial numbers to
% reject into bottom box

% NOTES 
% * When data is loaded, 'bad channels' will be automatically greyed out. Type
% 'badchans' into the appropriate box to view them after refreshing the
% data
% * Viewing the channels as interpolated will interpolate them for real!
% Best practise is to view the channels first and revert as necessary
% before 
close all
eeglab

%% Inside function
global EEG
EEG=EEG_in;
prevEEG = EEG; % previous EEG version
origEEG = EEG; % will be original EEG version
ALLEEG = EEG;

if isfield(S.(S.func),'inclchan')
    chan = S.(S.func).inclchan;
    exchan = 1:EEG.nbchan; exchan(chan)=[];
else
    chan = 1:EEG.nbchan;
    exchan=[];
end

% calculate bad channels
hideElec = [];
% [~, indElec] = pop_rejchan(EEG, 'elec',[1:EEG.nbchan],...
%     'threshold',5,'norm','on','measure','kurt');
if ~isempty(z_thresh) && iscell(z_thresh) 
    for i = 1:length(S.(S.func).clean.man.freq)
        
        if isempty(S.(S.func).clean.man.freq{i})
            filterset = [0 0];
        else
            filterset = S.(S.func).clean.man.freq{i};
        end
        
        if filterset(1)>0; EEGf = pop_eegfiltnew( EEG, filterset(1), 0, [], 0);end
        if filterset(2)>0; EEGf = pop_eegfiltnew( EEGf, 0, filterset(2), [], 0);end
        eeg_2d = reshape(EEGf.data,EEGf.nbchan,[]);
        eeg_2d = eeg_2d(chan,:);
        level = std(eeg_2d,[],2);
        minvar = 1; % uV

        % 
        flat = find(level<minvar);
        % use simelp threshold
        outliers = [];
        % use variance in std scores and set second threshold
        if S.(S.func).clean.man.interactive
            f(i)=figure;
        end
        if z_thresh{i}(2)
            thresh2=(mean(level(level>minvar))+z_thresh{i}(2)*std(level(level>minvar)));
            outliers = unique([outliers;find(level>thresh2)]);
            if S.(S.func).clean.man.interactive
                hold on; plot(level); scatter(1:length(level),level,25,'b','filled'); plot([0,length(level)],[thresh2,thresh2],'y')
                title(['freq: ' num2str(filterset) ', thresh1=red, thresh2=yellow'])
            end
        end
        if z_thresh{i}(1)
            level_without =level;
            level_without(outliers)=[];
            new_out=find(level>mean(level_without)*z_thresh{i}(1));
            outliers = [outliers; new_out];
            if S.(S.func).clean.man.interactive
                hold on; plot([0,length(level)],mean(level_without)*[z_thresh{i}(1),z_thresh{i}(1)],'r')
                title(['freq: ' num2str(filterset) ', thresh1=red, thresh2=yellow'])
            end
        end
        S.(S.func).clean.man.rejchan = chan(outliers)';

        if S.(S.func).clean.man.hidebad
            % hide the outliers, flat channels and excluded channels
            hideElec= unique([hideElec sort([chan(flat) chan(outliers) exchan])]);
        else
            % display the outliers
            temp= 1:EEGf.nbchan; 
            temp(chan(outliers))=[];
            if ~isempty(hideElec)
                hideElec = intersect(hideElec,temp);
            else
                hideElec = temp;
            end
        end
    end
    
end

if numel(hideElec)>=1
    % add to EEG structure as badchan
    bchan = zeros(1,EEG.nbchan);
    bchan(hideElec) = 1;
    bchan = num2cell(bchan);
    [EEG.chanlocs.badchan] = bchan{:};
else
    fprintf('No bad channels to plot...\n')
end

bchan = num2cell(zeros(1,EEG.nbchan));

b = []; % initialise variable

% view data and remove some trials
% breakloop = 0;
% while ~breakloop
    %% Trial selection GUI
    % Computes 'bad' channels using kurtosis and displays them in grey against
    % blue 'good' channels 
    prevEEG2 = prevEEG;
    prevEEG = EEG;
    
    if S.(S.func).clean.man.interactive
        pop_eegplot(EEG, 1, 1, 0); % view data first
        ud=get(gcf,'UserData');     %2
        ud.winlength=10;            %3
        set(gcf,'UserData',ud);     %4
        eegplot('draws',0)          %5

        h = gcf;
        waitfor(h)
    end
    
    tempRej = EEG(1).reject; % keep rejection info
    
    if isempty(tempRej.rejmanual); tempRej = ALLEEG(1).reject; end % try ALLEEG if there are no rejection markers
    if isempty(tempRej.rejmanual); tempRej = ALLEEG(end).reject; end % try other ALLEEG if there are no rejection markers

    manrejlist = find(tempRej.rejmanual); % trials manually marked for rejection
    allchan = {EEG.chanlocs.labels}; % list of channelnames
    
    S.(S.func).clean.man.rejtrial = manrejlist';
    
%     %% Identify channels/trials to remove
%     prompt = {'Enter channels to view removal via interpolation (note: type ''PREVIOUS'' or ''ORIGINAL'' to reset changes):',...
%         'Enter trials to view removal (type ''manrej'' for pre-selected trials):',...
%         'Highlight channel(s) (type ''badchans'' for pre-selected channels):',...
%         'Fade channel(s) (type ''badchans'' for pre-selected channels):',...
%         'Channels to remove via interpolation','Trials to remove (type ''manrej'' for pre-selected trials)'};
%     title = 'Manual Trial/Channel Rejection for ICA';
%     dims = [1 70];
%     definput = {'','','','','',''};
%     answer = inputdlg(prompt,title,dims,definput);
%     
%     %% no changes
%     if isempty(answer)
%         fprintf('No changes, so this will be stored as the current dataset.\n')
%         breakloop = 1;
%         continue
%     end
%     
%     %% changes (real)
%     % interpolate channels (real)
%     if  ~isempty(answer{5})
%         str = answer{5};
%         gap = find(str==' ');
%         if ~isempty(gap); gap = [1 gap length(str)]; end
%         rejchans = {str}; % initialise as single answer
%         for k = 1:length(gap)-1
%             temp = str(gap(k):(gap(k+1)));
%             rejchans{k} = strtrim(temp); % remove whitespace
%         end
%         a = upper(allchan);
%         b = upper(rejchans);
%         chanloc = find(ismember(a,b));
%         EEG= eeg_interp(EEG, chanloc, 'spherical');
%         breakloop = 1;
%     end
%     
%     if  ~isempty(answer{6})         % remove trials (real)
%         if ischar(answer{6}) % use marked trials from 'EEG' struct
%             rej = manrejlist;
%         else
%             rej = str2num(answer{6});
%         end
%         fprintf('removing %g trials\n',numel(rej))
%         EEG = pop_rejepoch(EEG, rej);
%         breakloop = 1; % end viewing process
%     end
%     
%     %% changes (view only)
%     % interpolate channels (view only)
%     if  ~isempty(answer{1}) && breakloop==0
%         str = answer{1};
%         gap = find(str==' ');
%         if ~isempty(gap); gap = [1 gap length(str)]; end
%         rejchans = {str}; % initialise as single answer
%         for k = 1:length(gap)-1
%             temp = str(gap(k):(gap(k+1)));
%             rejchans{k} = strtrim(temp); % remove whitespace
%         end
%         a = upper(allchan);
%         b = upper(rejchans);
%         chanloc = find(ismember(a,b));
%         
%         % either: revert to original dataset, interpolate channels
%         % (temporarily)
%         if strcmpi(b,'PREVIOUS')
%             fprintf('Restoring EEG structure to previous state...\n')
%             EEG = prevEEG2; ALLEEG = EEG; tempRej = EEG(1).reject; b = [];
%         elseif strcmpi(b,'ORIGINAL')
%             fprintf('Restoring EEG structure to original state...\n')
%             EEG = origEEG; ALLEEG = EEG; tempRej = EEG(1).reject; b = [];
%         elseif ~isempty(chanloc)
%             EEG= eeg_interp(EEG, chanloc, 'spherical');
%             
%             % check for removed channels, update EEG structure if so (UNUSED)
%             Match=cellfun(@(x) ismember(x, {EEG.chanlocs.labels}), allchan, 'UniformOutput', 0);
%             r=find(~cell2mat(Match)); % find indices of missing channel
%             if ~isempty(r)
%                 warnmsg = sprintf('Removing %g channels rather than interpolating...',numel(r));
%                 warning(warnmsg)
%                 EEG.chanlocs(r) = []; % remove unwanted channel
%             end
%         end
%     end
%     
%     % remove trials (view only)
%     if  ~isempty(answer{2}) && breakloop==0
%         if strcmpi(answer{2},'MANREJ') % use marked trials from 'EEG' struct
%             rej = manrejlist;
%         else
%             rej = str2num(answer{2});
%         end
%         EEG = pop_rejepoch(EEG, rej);
%     end
%     
%     
%     %% highlight/fade channels (view only)
%     
%     % reset badchans
%     bchan = num2cell(zeros(1,EEG.nbchan));
%     [EEG.chanlocs.badchan] = bchan{:};
%     
%     if  ~isempty(answer{3}) || ~isempty(answer{4}) && breakloop==0
%         
%         if ~isempty(answer{3}) % highlight
%             highFade = answer{3};
%             hf = 'high';
%         elseif ~isempty(answer{4}) % highlight
%             highFade = answer{4};
%             hf = 'fade';
%         end
%         
%         % locate channels
%         str = highFade;
%         gap = find(str==' ');
%         if ~isempty(gap); gap = [1 gap length(str)]; end
%         hfchans = {str}; % initialise as single answer
%         for k = 1:length(gap)-1
%             temp = str(gap(k):(gap(k+1)));
%             hfchans{k} = strtrim(temp); % remove whitespace
%         end
%         
%         a = upper(allchan);
%         b = upper(hfchans);
%         hideElec= find(ismember(a,b));
%         
%         % check if either say 'badchans' (this will overwrite previous hideElecvalue)
%         if strcmpi(answer{3},'BADCHANS') || strcmpi(answer{4},'BADCHANS') % use marked trials from 'EEG' struct
%             [~, hideElec] = pop_rejchan(EEG, 'elec',[1:EEG.nbchan],...
%                 'threshold',5,'norm','on','measure','kurt');
%             
%             if strcmpi(answer{3},'BADCHANS') % highlight bad channels
%                 hf = 'high';
%             elseif strcmpi(answer{4},'BADCHANS') % fade bad channels
%                 hf = 'fade';
%             end
%         end
%             
%         % update bad channels list based on structure
%             if strcmpi(hf,'high') % highlight bad channels
%                 bchan = ones(1,EEG.nbchan);
%                 bchan(hideElec) = 0;
%             elseif strcmpi(hf,'fade') % fade bad channels
%                     bchan = zeros(1,EEG.nbchan);
%                     bchan(hideElec) = 1;
%             end
%                 
%                 bchan = num2cell(bchan);
%                 [EEG.chanlocs.badchan] = bchan{:};
%     end
    
    %% Update rejection markers
    % keep rejection markers updated as long as user hasn't
    % selected 'revert' or 'original'
%     if any(strcmpi(b,'PREVIOUS')) || any(strcmpi(b,'ORIGINAL'))
%     elseif breakloop==0
%         EEG.reject = tempRej; % restore rejection marks
%     end
% end