function S = autocorr_channel_reject(S,EEG)
% function maximises the amount of data remaining after removing autocorr
% channels (on certain trials)

% Planned future updates:
% add topographic plots showing locations of channels to be interpolated.

% OLD
% std threshold - less variance than this per trial will be rejected
% varthresh = S.prep.clean.flatchan.varthresh;

dropautocorr = S.prep.clean.autochan.dropautocorr; 
autocorrint = S.prep.clean.autochan.autocorrint;% will compute autocorrelation with this many milliseconds lag

results = [];
idx={};

Ncorrint=round(autocorrint/(1000/EEG.srate)); % number of samples for lag
rej = false(1,EEG.nbchan);
for k=1:EEG.nbchan
    for t = 1:EEG.trials
        yy=xcorr(EEG.data(k,:,t),Ncorrint,'coeff');
        autocorr(k,t) = yy(1);
    end
end

% any need excluding?
if sum(autocorr<dropautocorr, 'all')==0
    S.prep.clean.autochan.rejchan = [];
    S.prep.clean.autochan.rejtrial = [];
    return
end

% plot
f1=figure;
tiledlayout(1,2)
nexttile
imagesc(autocorr,[0 1]); title(['autocorr'])
nexttile
imagesc(autocorr<dropautocorr,[0 1]); title(['under threshold of ' num2str(dropautocorr)])

% channels with little autocorrelation as fraction of trials
trialfrac = sum(autocorr<dropautocorr,2)/EEG.trials;
[~,si] = sort(trialfrac,'descend');
autocorr_trialfrac = S.prep.clean.autochan.autocorr_trialfrac;
S.prep.clean.autochan.rejchan = find(isoutlier(trialfrac,'median') & trialfrac>autocorr_trialfrac & trialfrac>median(trialfrac));
chans= {EEG.chanlocs.labels};
S.prep.clean.autochan.chans_reject_labels_ordered = chans(si(1:length(S.prep.clean.autochan.rejchan)));
S.prep.clean.autochan.rejtrial = [];
S.prep.clean.autochan.median_frac_trials = median(trialfrac);

if 0
    % order over trials
    [~,trialind]=sort(sum(autocorr,1),'ascend');
    
    % cycle through every combination of trial and channel removal
    row=0;
    for t = 1:length(trialind)
        
        % remove trials
        rmtrial = autocorr;
        rmtrial(:,trialind(1:t)) = []; 
        
        % order over trials
        [~,chanind]=sort(sum(rmtrial,2),'ascend');
        
        for c = 1:length(chanind)
            
            % remove chans
            rmchan = rmtrial;
            rmchan(chanind(1:c),:) = []; 
            
            % if all removed, calculate area left
            if sum(rmchan<dropautocorr,'all')==0
                row=row+1;
                results(row,:) = [t,c,length(rmchan(:))];
                idx{row} = {trialind(1:t),chanind(1:c)};
                continue
            end
        end
    end
    
    trial_weight = S.prep.clean.autochan.trial_weight;
    chan_weight = S.prep.clean.autochan.chan_weights;
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
        rmvar = zeros(size(autocorr));
        rmvar(ordered_idx{1}{2},:) = 1;
        rmvar(:,ordered_idx{1}{1}) = 1;
        subplot(ceil(length(chan_weight)/3),ceil(length(chan_weight)/4),m);
        imagesc(1-rmvar)
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
        S.prep.clean.autochan.rejchan = [];
        S.prep.clean.autochan.rejtrial = [];
        close(f2)
    else
        [~,order] = sort(metric(:,answer),'descend');
        ordered_results = results(order,:);
        ordered_idx=idx(order);
        S.prep.clean.autochan.rejchan = sort(ordered_idx{1}{2});
        S.prep.clean.autochan.rejtrial = sort(ordered_idx{1}{1})';
        close(f2)
    end
end
close(f1)