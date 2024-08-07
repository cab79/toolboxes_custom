function [varargout] = eegstats_encode_LMM(varargin)
% Runs Linear Mixed Model on one or a number of data samples (s)

if isempty(varargin)
    % assume Condor
    load input.mat;
    Y=struct;
    for s = 1:size(tempY,2)
        temp = array2table(tempY(:,s),'VariableNames',{'data'});
        Y(s).dtab = horzcat(dtab,temp);
    end
else
    % assume S,Y inputs
    S=varargin{1};
    Y=varargin{2};
end

LM=length(S.encode.model);
MC=length(S.encode.model_compare);
disp(['number of models:' num2str(LM)]);
disp(['number of model comparisons:' num2str(MC)]);

M = struct;
failed_s = cell(1,LM);
for s = 1:length(Y)

    % % Display the progress
    % progress = (s / length(Y)) * 100;
    % fprintf('Progress: %3.1f%%\r', progress);
    % % Use carriage return to overwrite the same line
    % if s < length(Y)
    %     fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'); % Number of \b should match the length of the progress text
    % else
    %     fprintf('\n'); % Move to the next line on the final iteration
    % end
    
    % for each model
    model=struct;
    for i = 1:LM

        % fit model
        try
            if ismember('train',Y(s).dtab.Properties.VariableNames)
                train_data = Y(s).dtab(find(double(Y(s).dtab.train)),:);
            else
                train_data = Y(s).dtab;
            end
            if S.encode.zscore
                try
                    mustBeNumericOrLogical(train_data.data)
                    train_data.data = applyz(train_data.data);
                catch
                    disp(['EEG data not numeric or logical for LM ' num2str(LM) ', sample ' num2str(s)])
                end
            end
            lmm=fitlme(train_data,S.encode.model{i},'FitMethod',S.encode.lmm.fitmethod,'DummyVarCoding', S.encode.lmm.coding);
            last_s_worked=s;
        catch
            failed_s{i} = [failed_s{i} s];
            if isfield(Y(s).dtab,'train')
                train_data = Y(last_s_worked).dtab(find(double(Y(last_s_worked).dtab.train)),:);
            else
                train_data = Y(last_s_worked).dtab;
            end
            if S.encode.zscore
                try
                    mustBeNumericOrLogical(train_data.data)
                    train_data.data = applyz(train_data.data);
                catch
                    disp(['EEG data not numeric or logical for LM ' num2str(LM) ', sample ' num2str(s)])
                end
            end
            lmm=fitlme(train_data,S.encode.model{i},'FitMethod',S.encode.lmm.fitmethod,'DummyVarCoding', S.encode.lmm.coding);
        end
        model(i).lmm = lmm; % temporary var for model comparisons
        [R,Rn]=randomEffects(lmm);
        
        % outputs common to all samples
        if s==1
            M.model(i).samples(1).def = char(lmm.Formula);
            M.model(i).samples(1).pred = lmm.PredictorNames;
            M.model(i).samples(1).fixeddesign = designMatrix(lmm,'Fixed');
            M.model(i).samples(1).randomdesign = designMatrix(lmm,'Random');
            M.model(i).samples(1).CoefficientNames = lmm.CoefficientNames;
            M.model(i).samples(1).RandomNames = Rn;
            M.model(i).samples(1).failed_s = failed_s{i};
        end
        
        % outputs for each sample
        M.model(i).samples(s).b = double(lmm.Coefficients(:,2));
        M.model(i).samples(s).r = R;
        M.model(i).samples(s).mse = lmm.MSE;
        M.model(i).samples(s).logl = lmm.LogLikelihood;
        M.model(i).samples(s).r2_ord = lmm.Rsquared.Ordinary;
        M.model(i).samples(s).r2_adj = lmm.Rsquared.Adjusted;
        M.model(i).samples(s).coeffcov = lmm.CoefficientCovariance;
        M.model(i).samples(s).psi = covarianceParameters(lmm);

        % contrasts
        if strcmp(S.encode.lmm.contrasts(1).Term,'anova') 
            anovastats=anova(lmm,'DFMethod','satterthwaite');
            M.model(i).samples(s).Term = anovastats.Term;
            M.model(i).samples(s).DF = [anovastats.DF1,anovastats.DF2];
            M.model(i).samples(s).F = anovastats.FStat;
            M.model(i).samples(s).p = anovastats.pValue;
        else
            for h=1:length(S.encode.lmm.contrasts)
                [pValue,FStat,DF1,DF2] = coefTest(lmm,S.encode.lmm.contrasts(h).H,'DFMethod','satterthwaite');
                M.model(i).samples(s).Term = S.encode.lmm.contrasts(h).Term;
                M.model(i).samples(s).DF(h,:) = [DF1,DF2];
                M.model(i).samples(s).F(h,1) = FStat;
                M.model(i).samples(s).p(h,1) = pValue;
            end
        end

        % residuals
        resid = residuals(lmm,'ResidualType','Standardized');
%         try
%             M.model(i).samples(s).hnorm=kstest(resid);
%             M.model(i).samples(s).skew=skewness(resid);
%             M.model(i).samples(s).kurt=kurtosis(resid);
%         catch
%             M.model(i).samples(s).hnorm=NaN;
%             M.model(i).samples(s).skew=NaN;
%             M.model(i).samples(s).kurt=NaN;
%         end
        if S.encode.save_residuals
%             resid = bsxfun(@times,resid,Y(s).data_std); % don't re-scale
            M.model(i).samples(s).resid=resid;
        end
        % fitted
        if S.encode.save_fitted
            ftd=fitted(lmm,'Conditional',true);
%             ftd = bsxfun(@plus,bsxfun(@times,ftd,Y(s).data_std),Y(s).data_mean); % don't re-scale
            M.model(i).samples(s).fitted=ftd;
        end
        % input
        if S.encode.save_input
            input=Y(s).dtab.data;
%             input = bsxfun(@plus,bsxfun(@times,input,Y(s).data_std),Y(s).data_mean); % don't re-scale
            M.model(i).samples(s).input=input;
        end
    end

    % compare models
    for i = 1:MC
        modi = S.encode.model_compare{i};
        M.model_comp(i).samples(s).pair=modi;
        try
            results = compare(model(modi(1)).lmm,model(modi(2)).lmm,'CheckNesting',true);
            M.model_comp(i).samples(s).pval = results.pValue(2); 
            M.model_comp(i).samples(s).LRStat = results.LRStat(2); 
            M.model_comp(i).samples(s).deltaDF = results.deltaDF;
        catch
           M.model_comp(i).samples(s).pval = nan;
           M.model_comp(i).samples(s).LRStat = nan; 
           M.model_comp(i).samples(s).deltaDF = [nan nan];
        end
    end
end

if isempty(varargin)
    % Condor output
    save('output.mat','M','S','chunk_info');
else
    varargout = {M};
end

function data = applyz(data)

data_mean = nanmean(data);
data_std = nanstd(data);
if data_std>0
    data = (data - data_mean) / data_std;
end