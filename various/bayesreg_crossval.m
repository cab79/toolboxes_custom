function [out,S,X] = bayesreg_crossval(X,Y,S,groups,varargin)
% Y can be a matrix with multiple columns, each is analysed separately.

if ~isempty(varargin)
    catvars = varargin{1};
else
    catvars=[];
end

for s = 1:size(Y,2)
    y = Y(:,s);
    
    % z-score non-categorical only
    contvars = 1:size(X,2);
    contvars(catvars)=[];
    if S.zscore
        [X(:,contvars), X_mean(:,contvars), X_stand_de(:,contvars)] = zscore(X(:,contvars), [], 1);
        if ~strcmp(S.brr.model,'binomial')
            [y, y_mean, Y_stand_de] = zscore(y);
        end
    end

    N=size(X,1);
    K=S.brr.folds;
    if K
        if N>K
            CVsample = cvpartition(N,'KFold',K);
        elseif N==K
            CVsample = cvpartition(N,'LeaveOut');
        end
%         x_val_in  = crossvalind('Kfold', N, K);

    else
        % if no folds, run all trials for both training and testing
        x_val_in = ones(N,1);
        K=1;
    end
    for i = K : -1 : 1
        
        if S.brr.folds
            val_in = test(CVsample,i);
        else
            val_in = x_val_in == i;
        end
        
        es_in     = ~val_in;
        if ~any(es_in); es_in=val_in; end

        if (any(isnan(y(es_in))))
            stt(i).muB = NaN*ones(size(X,2),1); 
            stt(i).muSigma2 = NaN; 
            logl(i) = NaN; 
            waic(i) = NaN; 
            r2(i) = NaN; 
            predstt(i).neglike = NaN; 
            predstt(i).r2 = NaN; 
        else
%             if S.brr.usegroups 
%                 grps = unique(groupvec);
%                 if length(grps)>1
%                     for g = grps
%                         groups{g} = find(groupvec==g);
%                     end
            if ~isempty(groups)
                if ~isempty(catvars)
                    [beta, beta0, stt(i)] = bayesreg(X(es_in,:),y(es_in),S.brr.model,S.brr.prior,'nsamples',S.brr.nsamples,'burnin',S.brr.burnin,'thin',S.brr.thin,'display',false,'waic', true, 'catvars',catvars, 'groups',groups);
                else
                    [beta, beta0, stt(i)] = bayesreg(X(es_in,:),y(es_in),S.brr.model,S.brr.prior,'nsamples',S.brr.nsamples,'burnin',S.brr.burnin,'thin',S.brr.thin,'display',false,'waic', true, 'groups',groups);
                end
            else
                if ~isempty(catvars)
                    [beta, beta0, stt(i)] = bayesreg(X(es_in,:),y(es_in),S.brr.model,S.brr.prior,'nsamples',S.brr.nsamples,'burnin',S.brr.burnin,'thin',S.brr.thin,'display',false,'waic', true, 'catvars',catvars);
                else
                    [beta, beta0, stt(i)] = bayesreg(X(es_in,:),y(es_in),S.brr.model,S.brr.prior,'nsamples',S.brr.nsamples,'burnin',S.brr.burnin,'thin',S.brr.thin,'display',false,'waic', true);
                end
            end
%             else
%                 [beta, beta0, stt(i)] = bayesreg(X(es_in,:),y(es_in),S.brr.model,S.brr.prior,'nsamples',S.brr.nsamples,'burnin',S.brr.burnin,'thin',S.brr.thin,'display',false,'waic', S.brr.waic,'catvars',catvars);
%             end
            try
                [pred, predstt(i)] = br_predict(X(val_in,:), beta, beta0, stt(i), 'ytest', y(val_in), 'CI', [2.5, 97.5], 'display', false);
            catch
                [pred, predstt(i)] = br_predict(X(val_in,:), beta, beta0, stt(i), 'ytest', y(val_in), 'display', false);
            end

            logl(i) = stt(i).modelstats.logl; 
            waic(i) = stt(i).modelstats.waic; 
            r2(i) = stt(i).modelstats.r2; 
        end
    end

    out(s).muB = nanmean([stt(:).muB],2);
    out(s).r2 = nanmean(r2);
    out(s).waic = nanmean(waic);
    out(s).logl = nanmean(logl);
    out(s).neglike = nanmean([predstt(:).neglike]);
    try
        out(s).muSigma2 = nanmean([stt(:).muSigma2]);
        out(s).r2test = nanmean([predstt(:).r2]);
    end
    try
        out(s).xmean = X_mean;
        out(s).xstd = X_stand_de;
        out(s).ymean = y_mean;
        out(s).ystd = Y_stand_de;
    end
end