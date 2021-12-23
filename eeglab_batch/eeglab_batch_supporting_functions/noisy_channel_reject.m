function S = noisy_channel_reject(S,EEG)
% function maximises the amount of data remaining after removing flat
% channels (on certain trials)

% Planned future updates:
% add topographic plots showing locations of channels to be interpolated.

% OLD
% std threshold - less variance than this per trial will be rejected
% varthresh = S.prep.clean.flatchan.varthresh;

results = [];
idx={};

% inverse variance over data points on each trial and channel
varv = squeeze(std(EEG.data,[],2));
% invvar = 1./squeeze(std(EEG.data,[],2));

% any need excluding?
if sum(varv<S.prep.clean.noisychan.varthresh, 'all')==0
    S.prep.clean.noisychan.rejchan = [];
    S.prep.clean.noisychan.rejtrial = [];
    return
end

% plot
f1=figure;
tiledlayout(1,2)
nexttile
imagesc(varv,[0 S.prep.clean.noisychan.varthresh]); title(['var up to ' num2str(S.prep.clean.noisychan.varthresh)])
nexttile
imagesc(varv>S.prep.clean.noisychan.varthresh,[0 1]); title(['over threshold of ' num2str(S.prep.clean.noisychan.varthresh)])

% channels with high variance as fraction of trials
trialfrac = sum(varv>S.prep.clean.noisychan.varthresh,2)/EEG.trials;
[~,si] = sort(trialfrac,'descend');
var_trialfrac = S.prep.clean.noisychan.trialfrac;
S.prep.clean.noisychan.rejchan = find(isoutlier(trialfrac,'median') & trialfrac>var_trialfrac & trialfrac>median(trialfrac));
chans= {EEG.chanlocs.labels};
S.prep.clean.noisychan.chans_reject_labels_ordered = chans(si(1:length(S.prep.clean.noisychan.rejchan)));
S.prep.clean.noisychan.rejtrial = [];
S.prep.clean.noisychan.median_frac_trials = median(trialfrac);

if 0

    
    % order variance over trials: trialind is an index of which trials have the
    % greatest variance
    [~,trialind]=sort(sum(varv,1),'descend');
    
    % cycle through every combination of trial and channel removal
    row=0;
    for t = 1:length(trialind)
        
        % remove trials
        rmtrial = varv;
        rmtrial(:,trialind(1:t)) = []; 
        
        % order variance over trials
        [~,chanind]=sort(sum(rmtrial,2),'descend');
        
        for c = 1:length(chanind)
            
            % remove chans
            rmchan = rmtrial;
            rmchan(chanind(1:c),:) = []; 
            
            % if all high var removed, calculate area left
            if sum(rmchan>S.prep.clean.noisychan.varthresh,'all')==0
                row=row+1;
                results(row,:) = [t,c,length(rmchan(:))];
                idx{row} = {trialind(1:t),chanind(1:c)};
                continue
            end
        end
    end
    
    trial_weight = S.prep.clean.noisychan.trial_weight;
    chan_weight = S.prep.clean.noisychan.chan_weights;
    % use area per lost channel/trial
    metric=[];
    for cw = 1:length(chan_weight)
        metric(:,cw) = results(:,3)./((results(:,2)*chan_weight(cw)) + (results(:,1)*trial_weight));
    end
    f2=figure('units','normalized','outerposition',[0 0 1 1])
    for m = 1:size(metric,2)
        [~,order] = sort(metric(:,m),'descend');
        ordered_results = results(order,:);
        ordered_idx=idx(order);
        rmvar = zeros(size(varv));
        rmvar(ordered_idx{1}{2},:) = 1;
        rmvar(:,ordered_idx{1}{1}) = 1;
        subplot(ceil(length(chan_weight)/3),ceil(length(chan_weight)/4),m);
        imagesc(rmvar)
        title(num2str(m));
        xlabel(['bad trials = ' num2str(ordered_results(1,1))]);
        ylabel(['bad chans = ' num2str(ordered_results(1,2))]);
    end
    
    if length(chan_weight)>1
        answer = str2double(inputdlg('Choose plot number','',1,{'7'}));
    else
        answer = 1;
    end
    if isnan(answer) || ~answer
        S.prep.clean.noisychan.rejchan = [];
        S.prep.clean.noisychan.rejtrial = [];
        close(f2)
    else
        [~,order] = sort(metric(:,answer),'descend');
        ordered_results = results(order,:);
        ordered_idx=idx(order);
        S.prep.clean.noisychan.rejchan = sort(ordered_idx{1}{2});
        S.prep.clean.noisychan.rejtrial = sort(ordered_idx{1}{1})';
        chans= {EEG.chanlocs.labels};
        chans_reject_labels = chans(S.prep.clean.noisychan.rejchan)
        close(f2)
    end
end
close(f1)