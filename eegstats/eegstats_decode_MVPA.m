function D=eegstats_decode_MVPA(S,varargin)

%% load prep data
if isempty(varargin)
   % try loading the data if not supplied in varargin
   try
       disp('loading prepped data')
       load(fullfile(S.mvpa.path.inputs.prep,'D.mat'),'D') 
   catch
       error('Supply a D structure or path/file of the data')
   end
else
    D= varargin{1};
end

%% set paths
S.path.code = {
    };
set_paths(S)

for d = 1:length(D)
            
    % get predictor and events
    pred = table2array(D(d).prep.dtab.(S.mvpa.pred));
    events = table2array(D(d).prep.dtab.eventType);
    
    % get indices
    trainidx=find(double(D(d).prep.dtab.train));
    testidx=find(double(D(d).prep.dtab.test));
    
    % data for decoding: (trials x voxels/samples)
    for s = 1:length(D(d).prep.Y)
        if s==1
            data = table2array(D(d).prep.Y(s).dtab);
        else
            data = horzcat(data,table2array(D(d).prep.Y(s).dtab));
        end
    end

    % mask the data
    mask = ones(size(data,2),1);
    if ~isempty(S.mvpa.maskimg)
        maskimg = spm_read_vols(spm_vol(S.mvpa.maskimg));
        mask = ungrid(maskimg,S)>0;
        mask=mask(:);
        data(:,~mask)=0;
    end
        
    % select training and testing data
    train_data = data(trainidx,:);
    test_data = data(testidx,:);

    % conditions for cosmo dataset
    if strcmp(S.mvpa.type,'class') 
        % classification: 0s and 1s over trials 
        conds = pred(trainidx);
    elseif strcmp(S.mvpa.type,'regress')  
        % regression: condition types to include
        conds = events(trainidx);
    end
    
    % create cosmo dataset
    cos = eeglab2cosmo(train_data,S.mvpa.select_timewin,conds); % UPDATE THIS FOR BRAINVISION DATA GETMARKERS

    % set the targets: predictor values
    if strcmp(S.mvpa.type,'class') 
        % classification
        cos.sa.targets=cos.sa.trialinfo(:,1); 
    elseif strcmp(S.mvpa.type,'regress') 
        % regression
        cos.sa.targets=pred(trainidx);
    end
    
    % set the chunks: each trial is a chunk
    cos.sa.chunks=[1:size(cos.samples,1)]'; 

    % use cosmo balancing
    if S.mvpa.balance_dataset_and_partitions 
        S.mvpa.balance_idx=[];
    end

    % run analysis: UPDATE FUNCTION TO S.MVPA
    [out,S] = run_cosmo_machine(cos,S);
    disp('MVPA complete')

    % prediction of test sample: dot product of weights with testdata
    test_pred=pred(testidx);
    if strcmp(S.mvpa.type,'class') % classification

        % variable to predict
        predicted=test_pred;
        out.predicted=predicted;

        % weight predictions
        out.testdata_pred = (out.weights(1,:) * test_data);
        if length(unique(predicted))==2 % if predicted consists of two values
            out.testdata_predbin=nan(1,length(out.testdata_pred));
            out.testdata_predbin(out.testdata_pred>0)=1;
            out.testdata_predbin(out.testdata_pred<0)=2; % this is opposite to intuition for GPC
            out.testdata_corr = -corr(out.testdata_pred',predicted,'type','Spearman');
            disp(['correlation: ' num2str(out.testdata_corr)])
            out.testdata_correct = sum(out.testdata_predbin==predicted')/length(out.testdata_pred);
            disp(['correct: ' num2str(out.testdata_correct)])
        else
            out.testdata_corr = corr(out.testdata_pred',predicted,'type','Spearman');
            disp(['correlation: ' num2str(out.testdata_corr)])
        end
        
        % transweight predictions
        out.testdata_predtrans = (out.transweights(1,:) * test_data);
        if length(unique(predicted))==2 % if predicted consists of two values
            out.testdata_predbintrans=nan(1,length(out.testdata_predtrans));
            out.testdata_predbintrans(out.testdata_predtrans>0)=1;
            out.testdata_predbintrans(out.testdata_predtrans<0)=2; % this is opposite to intuition for GPC
            out.testdata_corrtrans = -corr(out.testdata_predtrans',predicted,'type','Spearman');
            disp(['correlation trans: ' num2str(out.testdata_corrtrans)])
            out.testdata_correcttrans = sum(out.testdata_predbintrans==predicted')/length(out.testdata_predtrans);
            disp(['correct trans: ' num2str(out.testdata_correcttrans)])
        else
            out.testdata_corrtrans = corr(out.testdata_predtrans',predicted,'type','Spearman');
            disp(['correlation trans: ' num2str(out.testdata_corrtrans)])
        end
    elseif strcmp(S.mvpa.type,'regress') % regression
        % weight predictions
        out.testdata_pred = (out.weights * test_data) + out.offsets;
        out.testdata_corr = corr(out.testdata_pred',test_pred,'type','Spearman');
        disp(['correlation weigh: ' num2str(out.testdata_corr)])
        % transweight predictions
        out.testdata_predtrans = (out.transweights * testdata) + out.offsets;
        out.testdata_corrtrans = corr(out.testdata_predtrans',test_pred,'type','Spearman');
        disp(['correlation trans: ' num2str(out.testdata_corrtrans)])
    end

    % reshape weights to 2D
    sizdat = size(data);
    out.weights = reshape(out.weights(1,:),sizdat(1),sizdat(2));
    out.transweights = reshape(out.transweights(1,:),sizdat(1),sizdat(2));
    
    % Store results
    D(d).MVPA.out = out;
end

save(fullfile(S.mvpa.path.outputs, 'D.mat'),'D');

function set_paths(S)
for p = 1:length(S.path.code)
    if S.path.code{p,1}
        addpath(genpath(S.path.code{p,2}));
    else
        addpath(S.path.code{p,2});
    end
end

function Y = ungrid(YY,S)
% inputs
SPMdata = spm_eeg_load(fullfile(S.mvpa.path.SPMdata,S.mvpa.file.SPMdata));
chanind = 1:size(SPMdata,1); % index of chans
ntime = size(YY,3);
n = size(YY,1); % dimension
Cel = spm_eeg_locate_channels(SPMdata, n, chanind);
for ti = 1:ntime
    for ch = chanind
        Y(ch,ti) = YY(Cel(ch,1),Cel(ch,2),ti);
    end
end