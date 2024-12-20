function [varargout] = eegstats_encode_BRR(varargin)

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

M = struct;
for s = 1:length(Y)

    % Display the progress
    progress = (s / length(Y)) * 100;
    fprintf('Progress: %3.1f%%\r', progress);
    % Use carriage return to overwrite the same line
    if s < length(Y)
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'); % Number of \b should match the length of the progress text
    else
        fprintf('\n'); % Move to the next line on the final iteration
    end
    
    % for each model
    for i = 1:LM

        % training data
        if ismember('train',Y(s).dtab.Properties.VariableNames)
            train_idx = find(double(Y(s).dtab.train));
        else
            train_idx = 1:height(Y(s).dtab);
        end
        train_data = Y(s).dtab.data(train_idx,:);
        
        
        % predictors
        if s==1
            [terms, groups, train_pred, categ] = pred_from_formula(S.encode.model{i},Y(s).dtab(train_idx,:));
        end
        
        % fit model
        Sb.brr = S.encode.brr;  
        Sb.zscore = S.encode.brr.z_score;
        if S.encode.brr.use_groups==0
            groups=[];
        end
        brr=bayesreg_crossval(train_pred,train_data,Sb,groups,categ);

        % outputs common to all samples
        M.model(i).samples(1).def = S.encode.model{i};
        M.model(i).samples(1).CoefficientNames = terms;
        M.model(i).samples(1).fixeddesign = train_pred; % store predictors
%         M.model(i).samples(1).CoefficientNames = brr.CoefficientNames;

        % outputs for each sample
        M.model(i).samples(s).b = brr.muB;
        M.model(i).samples(s).t = brr.tStat;
        M.model(i).samples(s).mse = brr.muSigma2;
        M.model(i).samples(s).logl = brr.logl;
        M.model(i).samples(s).waic = brr.waic;
        M.model(i).samples(s).neglike = brr.neglike;
        M.model(i).samples(s).r2_ord = brr.r2test;

        % df and p values
        df = size(train_pred, 1) - size(train_pred, 2);
        M.model(i).samples(s).df = df; 
        M.model(i).samples(s).p = 2 * (1 - tcdf(abs(brr.tStat), df));
        
        ftd = train_pred*brr.muB; % Predicted responses at each data point.
        resid = train_data-ftd; % Residuals.
        M.model(i).samples(s).ktest_normality = kstest(resid);
        % residuals
        if S.encode.save_residuals
            M.model(i).samples(s).resid=resid;
        end
        % fitted
        if S.encode.save_fitted
            M.model(i).samples(s).fitted=ftd;
        end
        % input
        if S.encode.save_input
            M.model(i).samples(s).input=train_data;
        end
    end

end

if isempty(varargin)
    % Condor output
    save('output.mat','M','S','chunk_info');
else
    varargout = {M};
end
    
function [terms, groups, pred, categ] = pred_from_formula(modeldef,dtab)

% parse the modeldef term into predictors
def=regexprep(modeldef, '\s+', ''); % remove spaces
terms = strsplit(def,'+'); % identify terms occuring after a '+'
terms(1)=[];

% % expand
maineffects = {};
interactions = {};
for tm = 1:length(terms) % find interactions
    terms{tm} = regexprep(terms{tm},'+','');
    if any(regexp(terms{tm}, '.[*].'))
        maineffects = strsplit(terms{tm},'*');
        interactions={};
        for i = 2:length(maineffects)
            allcombs = nchoosek(maineffects,i);
            for ac = 1:size(allcombs,1)
                interactions = [interactions {strjoin(allcombs(ac,:),':')}];
            end
        end
        terms{tm}='';
    end
end
terms = [terms maineffects interactions];
terms = unique(terms,'stable');

% if there are interactions, obtain predictors and dummy coding from LMM
if ~isempty(interactions)
    lmm=fitlme(dtab,modeldef); % MUST use reference coding for Bayesreg to recognise cat variables
    pred=designMatrix(lmm,'Fixed');
    terms=lmm.CoefficientNames;
    pnames=lmm.PredictorNames;
else % otherwise don't use LMM, to save time
    pnames = terms;
    for p = 1:length(terms)
        pred(:,p) = double(dtab.(terms{p}));
    end
end
% remove intercept if present
rm = strcmp(terms,'(Intercept)');
pred(:,rm)=[]; terms(rm)=[];
% get grouped and categorical predictors
for pn=1:length(pnames)
    groups{pn} = find(contains(terms,pnames{pn}));
end
for tm = 1:length(terms)
    categ(tm)=all(ismember(unique(pred(:,tm)),[0 1]));
end
categ=find(categ);

% subterms={};
% for tm = 1:length(terms) % find interactions
%     if any(regexp(terms{tm}, '.:.'))
%         subterms{tm} = strsplit(terms{tm},':');
%     else
%         subterms{tm} = terms{tm};
%     end
% end

% generate new interaction variables
% pred=[];
% categ=zeros(1,length(terms));
% for tm = 1:length(subterms) 
%     if iscell(subterms{tm}) % interactions
%         temp=[];
%         subcateg=[];
%         for stm = 1:length(subterms{tm})
%             var=dtab.(subterms{tm}{stm});
%             subcateg(stm)=iscategorical(var);
%             if subcateg(stm)
%                 temp(:,stm)=double(var)-1;
%             else
%                 temp(:,stm)=var;
%             end
%         end
%         pred(:,tm) = prod(temp,2);
%         if all(subcateg)
%             categ(tm)=1;
%         end
%     else
%         var=dtab.(subterms{tm});
%         categ(tm)=iscategorical(var);
%         if categ(tm)
%             pred(:,tm) = double(var)-1;
%         else
%             pred(:,tm) = var;
%         end
%     end
% end