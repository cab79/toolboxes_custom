function [T,traj,param,rt] = HGF_extract_results(D,S)

if ~isfield(S,'summary_stats')
    S.summary_stats = {};
end

% fields to include/exclude
fitfields = {'p_prc','p_obs','traj'};
%fitfields = {'p_prc','traj'};
ignorefields = {'p','ptrans'};

% compile parameters and trajectories into  param and traj structures
for d = 1:length(D) 
    for ff = 1:length(fitfields)
        subfields = fieldnames(D(d).HGF.fit.(fitfields{ff}));
        subfields = subfields(~ismember(subfields,ignorefields));
        for sf = 1:length(subfields)
            if strcmp(fitfields{ff},'traj') && isstruct(D(d).HGF.fit.(fitfields{ff}).(subfields{sf}))
                trajfields = fieldnames(D(d).HGF.fit.(fitfields{ff}).(subfields{sf}));
                for tf = 1:length(trajfields)
                    traj(d).([subfields{sf} '_' trajfields{tf}]) = D(d).HGF.fit.(fitfields{ff}).(subfields{sf}).(trajfields{tf});
                end
            elseif strcmp(fitfields{ff},'traj')
                traj(d).(subfields{sf}) = D(d).HGF.fit.(fitfields{ff}).(subfields{sf});
            elseif strcmp(fitfields{ff},'p_prc') || strcmp(fitfields{ff},'p_obs')
                param(d).(subfields{sf}) = D(d).HGF.fit.(fitfields{ff}).(subfields{sf});
                if all(param(d).(subfields{sf}))==0
                    param(d).(subfields{sf})(:)=NaN;
                end
            end
        end
    end
end

% for param, add values to table
clear tablecols
tablecols.subjects = {D(:).subname}';
paramfields = fieldnames(param);
for pf = 1:length(paramfields)
    if length(param(1).(paramfields{pf}))>1
        for pfn = 1:length(param(1).(paramfields{pf}))
            for d = 1:length(param)
                tablecols.([paramfields{pf} '_' num2str(pfn)])(d,1) = param(d).(paramfields{pf})(pfn);
            end
        end
    else
        for d = 1:length(param)
            tablecols.(paramfields{pf})(d,1) = param(d).(paramfields{pf});
        end
    end
end

%% ---------------------------------------------------------------------
%  Convert all *_xi parameters into soft-max weights (*_phi)
% ----------------------------------------------------------------------
xiFields  = paramfields( endsWith(paramfields,'_xi') );   % e.g. AL_xi, PL_xi …
nXi       = numel(xiFields);
if nXi>0
    for d = 1:length(param)                               % loop subjects
        % 1) collect this subject's ξ values  (assumed scalar each)
        xiVec = zeros(1,nXi);
        for k = 1:nXi
            xiVal = param(d).(xiFields{k});
            xiVec(k) = xiVal(1);                          % take first elem if vector
        end
        xiVec( isnan(xiVec) ) = 0;                        % ← treat NaN as 0 ❶

        % 2) soft-max 
        expVec  = exp(xiVec);
        phiVec  = expVec ./ sum(expVec);                  % weights sum to 1

        % 3) add φ fields to tablecols
        for k = 1:nXi
            phiName = strrep(xiFields{k},'_xi','_phi');   % e.g. AL_phi
            tablecols.(phiName)(d,1) = phiVec(k);
        end
    end
end


% create results table
T = S.designtab;
Tp = struct2table(tablecols);
Tp(~ismember(Tp.subjects,T.subjects),:)=[];
T(~ismember(T.subjects,Tp.subjects),:)=[];
T = join(T,Tp,'Keys','subjects');
% T = [T,struct2table(tablecols)];
% if length(D)>1
%     T = cell2table(S.designmat(2:end,:),...
%         'VariableNames',S.designmat(1,:));
%     T = [T,struct2table(tablecols)];
% else
%     T=struct2table(tablecols);
% end


% for traj, first split into levels, then create means and stds, and then add values to table
clear tablecols
trajfields = fieldnames(traj);
% first split into levels if needed
for tf = 1:length(trajfields)
    if size(traj(1).(trajfields{tf}),2)>1
        for tfn = 1:size(traj(1).(trajfields{tf}),2)
            for d = 1:length(traj)
                traj(d).([trajfields{tf} '_' num2str(tfn)]) = traj(d).(trajfields{tf})(:,tfn);
            end
        end
        traj = rmfield(traj,trajfields{tf});
    end
end
trajfields = fieldnames(traj);
for tf = 1:length(trajfields)
    for d = 1:length(traj)
    % create means and stds
        if any(strcmp(S.summary_stats,'trial_mean'))
            tablecols.([trajfields{tf} '_trial_mean'])(d,1) = nanmean(traj(d).(trajfields{tf}));
        end
        if any(strcmp(S.summary_stats,'trial_std'))
            tablecols.([trajfields{tf} '_trial_std'])(d,1) = nanstd(traj(d).(trajfields{tf}));
        end
        if any(strcmp(S.summary_stats,'trial_absmean'))
            tablecols.([trajfields{tf} '_trial_absmean'])(d,1) = nanmean(abs(traj(d).(trajfields{tf})));
        end
        if any(strcmp(S.summary_stats,'trial_absstd'))
            tablecols.([trajfields{tf} '_trial_absstd'])(d,1) = nanstd(abs(traj(d).(trajfields{tf})));
        end
    end
end

if isfield(S,'condmean') && ~isempty(S.condmean)
    if size(S.condmean,1)>1
        S.condmean = S.condmean(:);
    end

    condfields = fieldnames(S.cond);
    for tf = 1:length(S.condmean) % for each traj
        for d = 1:length(traj) % for each subject
            conds = D(d).dt.design(2,:);
            thistraj = traj(d).(S.condmean{tf});
            uni_cond = unique(conds);
            uni_cond(uni_cond==0) = [];
            all_conds = [];
            for c = 1:length(uni_cond) % for each cond
                all_conds(uni_cond(c)) = nanmean(thistraj(ismember(conds,uni_cond(c))));
            end
            for cf = 1:length(condfields) % for each average
                if iscell(S.cond.(condfields{cf}))
                    if length(S.cond.(condfields{cf}))==2
                        tablecols.([S.condmean{tf} '_' condfields{cf}])(d,1) = nanmean(all_conds(S.cond.(condfields{cf}){1}) - nanmean(all_conds(S.cond.(condfields{cf}){2})));
                    elseif length(S.cond.(condfields{cf}))==4
                        tablecols.([S.condmean{tf} '_' condfields{cf}])(d,1) = (nanmean(all_conds(S.cond.(condfields{cf}){1})) - nanmean(all_conds(S.cond.(condfields{cf}){2}))) - (nanmean(all_conds(S.cond.(condfields{cf}){3})) - nanmean(all_conds(S.cond.(condfields{cf}){4})));
                    end
                else
                    tablecols.([S.condmean{tf} '_' condfields{cf}])(d,1) = nanmean(all_conds(S.cond.(condfields{cf})));
                end
            end
        end
    end
end

try
    T = [T,struct2table(tablecols)];
end

%get RTs too
if isfield(D(1).HGF,'y')
    for d = 1:length(D)
        rt(d).rt = D(d).HGF.y(find(D(d).HGF.u(:,1)),2);
    end
else
    for d = 1:length(D)
        rt(d).rt = [];
    end
end

% 
% rmv_out=100; % values that are a factor of rmv_out * median are removed from trajectories
% results = orderfields(results);
% nsub =length(subs);
% 
% num_stim = 1; % number of event markers per trial; e.g. anticipation cue AND laser stimulus = 2 stim
% ncond = max(ci); % two hands * two DC
% 
% var = {'al0','al1','rb','om','dau','da','ud','wt','psi','epsi','mu','sa','muhat','sahat','AIC','BIC','LME','irr','be0','be1','be2','be3','be4','be5','ze','be'};
% fol = {'p_prc','p_prc','p_prc','p_prc','traj','traj','traj','traj','traj','traj','traj','traj','traj','traj','optim','optim','optim','irr','p_obs','p_obs','p_obs','p_obs','p_obs','p_obs','p_obs','p_obs'};
% num_stimtypes = [num_stim*ones(1,10) ones(1,8) ones(1,8)];
% % 3lev
% if num_lev==3
%     numvar = [2, 2, 0, 3, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]; %3 level model
%     num_cond = [1, 1, 1, 1, 1, 1, 1, ncond, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]; %num conds to extract. 1= average across conditions, >1 = extract each condition
%     incl_abs = [0 0 0 0 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
%     split_input = [0 0 0 0 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
% elseif num_lev==2
% % 2lev, e.g. SDT
%     numvar = [2, 2, 1, 0, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0];
%     num_cond = [1, 1, 1, 1, 1, 1, 1, ncond, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]; %num conds to extract. 1= average across conditions, >1 = extract each condition
%     incl_abs = [0 0 0 0 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
%     split_input = [0 0 0 0 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
% elseif num_lev==1
% % 1lev, e.g. KF
%     numvar = [2, 2, 0, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]; %3 level model
%     num_cond = [1, 1, 1, 1, 1, 1, 1, ncond, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]; %num conds to extract. 1= average across conditions, >1 = extract each condition
%     incl_abs = [0 0 0 0 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
%     split_input = [0 0 0 0 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
% end
% 
% nme_list = {};
% sub_sort=1:nsub; % sort order of subjects
% grp = pdata(1+files_ana',grp_col);
% 
% % remove variables
% var(numvar==0)=[];
% fol(numvar==0)=[];
% num_stimtypes(numvar==0)=[];
% num_cond(numvar==0)=[];
% incl_abs(numvar==0)=[];
% split_input(numvar==0)=[];
% numvar(numvar==0)=[];
% 
% mu_traj = nan(nsub,max(numvar),length(results(1).u));
% 
% % create variables and list of variable names
% V=struct;
% for v=1:length(var)
%     for nv = 1:numvar(v)
%         for ns = 1:num_stimtypes(v)
%             if numvar(v)>1
%                 if num_stimtypes(v)>1
%                     nme = [var{v} num2str(ns) '_' num2str(nv)];
%                 else
%                     nme = [var{v} num2str(nv)];
%                 end
%             else
%                  if num_stimtypes(v)>1
%                     nme = [var{v} num2str(ns)];
%                  else
%                     nme = var{v};
%                  end
%             end
%             eval([nme ' = nan(nsub,1);']);
%             %V.(nme) = nan(nsub,1);
%             nme_list{length(nme_list)+1} = [nme '_av']; 
%             %V.(nme).name = [nme '_av'];
%             if incl_abs(v)==1
%                 nme_list{length(nme_list)+1} = [nme '_abs']; 
%                 %V.(nme).name = [nme '_abs'];
%             end
%             if split_input(v)==1
%                 inp = unique(results(1).u(:,1));
%                 inp(isnan(inp))=[];
%                 for i = 1:length(inp)
%                     nme_list{length(nme_list)+1} = [nme '_i' num2str(inp(i))]; 
%                 end
%             end
% 
%             if num_cond(v)>1
%                 for nc = 1:num_cond(v)
%                     nme_list{length(nme_list)+1} = [nme '_' num2str(nc)];
%                 end
%             end
%             if incl_abs(v)==1
%                 if num_cond(v)>1
%                     for nc = 1:num_cond(v)
%                         nme_list{length(nme_list)+1} = [nme '_' num2str(nc) '_abs']; 
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% exp = nan(nsub,3);
% for s = 1:nsub;
%     sub = subs{s};
%     cond=[];
%     if isfield(results,'conds')
%         for c = 1:ncond
%             cii = find(ci==c);
%             for ciii = 1:length(cii)
%                 ciiii = results(s).conds==cii(ciii);
%                 cond(ciiii) = c;
%             end
%         end
%     end
%     stimlen = length(results(s).traj.mu);
%     
%     % indices of each stimulus type (e.g. anticipation cue and pain
%     % stimulus)
%     inds = [];
%     for i = 1:num_stim
%         inds(i,:) = i:num_stim:(stimlen-num_stim+i);
%     end
%     
%     for nv = 1:numvar(strcmp(var,'mu'))
%         a = results(s).traj.mu(:,nv);
%         mu_traj(s,nv,1:length(a)) = a;
%     end
%     mu_traj=mu_traj(sub_sort,:,:);
%     
%     for v=1:length(var) % for each variable
%         for nv = 1:numvar(v) % for level/type of that variable
%             for ns = 1:num_stimtypes(v) % and each stimulus type
%                 if numvar(v)>1
%                     if num_stimtypes(v)>1
%                         nme = [var{v} num2str(ns) '_' num2str(nv)];
%                         if strcmp(fol{v},'p_obs') || strcmp(fol{v},'p_prc') || strcmp(fol{v},'optim')
%                             eval([nme '(s,1) = results(s).(fol{v}).' [var{v} num2str(ns)] '(nv);']);
%                         elseif strcmp(fol{v},'traj');
%                             eval(['a = results(s).(fol{v}).' var{v} '(inds(ns,:),nv);']);
%                             %a(a>abs(rmv_out*median(a)))=[];
%                             %a(a<-abs(rmv_out*median(a)))=[];
%                             if num_cond(v)>1
%                                 for nc = 1:num_cond(v)
%                                     eval([nme '(s,nc) = nanmean(a(cond==nc));']);
%                                 end
%                             elseif split_input(v)==1
%                                 for i = 1:length(inp)
%                                     eval([nme '(s,i) = nanmean(a(results(s).u(:,1)==inp(i)));']);
%                                 end
%                             else
%                                 eval([nme '(s,1) = nanmean(a,1);']);
%                             end
%                         end
%                     else
%                         nme = [var{v} num2str(nv)];
%                         if strcmp(fol{v},'p_obs') || strcmp(fol{v},'p_prc') || strcmp(fol{v},'optim')
%                             eval([nme '(s,1) = results(s).(fol{v}).' var{v} '(nv);']);
%                         elseif strcmp(fol{v},'traj');
%                             eval(['a = results(s).(fol{v}).' var{v} '(:,nv);']);
%                             %a(a>abs(rmv_out*median(a)))=[];
%                             %a(a<-abs(rmv_out*median(a)))=[];
%                             if num_cond(v)>1
%                                 for nc = 1:num_cond(v)
%                                     eval([nme '(s,nc) = nanmean(a(cond==nc));']);
%                                 end
%                             elseif split_input(v)==1
%                                 for i = 1:length(inp)
%                                     eval([nme '(s,i) = nanmean(a(results(s).u(:,1)==inp(i)));']);
%                                 end
%                             else
%                                 eval([nme '(s,1) = nanmean(a,1);']);
%                             end
%                         end
%                     end
%                 else
%                      if num_stimtypes(v)>1
%                         nme = [var{v} num2str(ns)];
%                         if strcmp(fol{v},'p_obs') || strcmp(fol{v},'p_prc') || strcmp(fol{v},'optim')
%                             eval([nme '(s,1) = results(s).(fol{v}).' [var{v} num2str(ns)] ';']);
%                         elseif strcmp(fol{v},'traj');
%                             eval(['a = results(s).(fol{v}).' var{v} '(inds(ns,:));']);
%                             %a(a>abs(rmv_out*median(a)))=[];
%                             %a(a<-abs(rmv_out*median(a)))=[];
%                             if num_cond(v)>1
%                                 for nc = 1:num_cond(v)
%                                     eval([nme '(s,nc) = nanmean(a(cond==nc));']);
%                                 end
%                             elseif split_input(v)==1
%                                 for i = 1:length(inp)
%                                     eval([nme '(s,i) = nanmean(a(results(s).u(:,1)==inp(i)));']);
%                                 end
%                             else
%                                 eval([nme '(s,1) = nanmean(a,1);']);
%                             end
%                         elseif strcmp(fol{v},'irr');
%                             eval(['a = results(s).' var{v} '(inds(ns,:));']);
%                             %a(a>abs(rmv_out*median(a)))=[];
%                             %a(a<-abs(rmv_out*median(a)))=[];
%                             eval([nme '(s,1) = length(a);']);
%                         end 
%                      else
%                         nme = var{v};
%                         if strcmp(fol{v},'p_obs') || strcmp(fol{v},'p_prc') || strcmp(fol{v},'optim')
%                             eval([nme '(s,1) = results(s).(fol{v}).' var{v} ';']);
%                         elseif strcmp(fol{v},'traj');
%                             eval(['a = results(s).(fol{v}).' var{v} ';']);
%                             %a(a>abs(rmv_out*median(a)))=[];
%                             %a(a<-abs(rmv_out*median(a)))=[];
%                             if num_cond(v)>1
%                                 for nc = 1:num_cond(v)
%                                     eval([nme '(s,nc) = nanmean(a(cond==nc));']);
%                                 end
%                             elseif split_input(v)==1
%                                 for i = 1:length(inp)
%                                     eval([nme '(s,i) = nanmean(a(results(s).u(:,1)==inp(i)));']);
%                                 end
%                             else
%                                 eval([nme '(s,1) = nanmean(a,1);']);
%                             end
%                         elseif strcmp(fol{v},'irr');
%                             eval(['a = results(s).' var{v} ';']);
%                             %a(a>abs(rmv_out*median(a)))=[];
%                             %a(a<-abs(rmv_out*median(a)))=[];
%                             eval([nme '(s,1) = length(a);']);
%                         end 
%                      end
%                 end
%             end
%         end
%     end
% 
% end
% 
% for v=1:length(var)
%     for nv = 1:numvar(v)
%         for ns = 1:num_stimtypes(v)
%             if numvar(v)>1
%                 if num_stimtypes(v)>1
%                     nme = [var{v} num2str(ns) '_' num2str(nv)];
%                 else
%                     nme = [var{v} num2str(nv)];
%                 end
%             else
%                  if num_stimtypes(v)>1
%                     nme = [var{v} num2str(ns)];
%                  else
%                     nme = var{v};
%                  end
%             end
%             
%             % average across conditions
%             eval([nme '_av = squeeze(nanmean(' nme ',2));']);
%             
%             % average within-condition
%             if num_cond(v)>1
%                 for nc = 1:num_cond(v)
%                     eval([nme '_' num2str(nc) '= squeeze(nanmean(' nme '(:,nc),2));']);
%                 end
%             end
%             
%             % average per input type
%             if split_input(v)==1
%                 for i = 1:length(inp)
%                     eval([nme '_i' num2str(inp(i)) '= squeeze(nanmean(' nme '(:,i),2));']);
%                 end
%             end
% 
%             % create absolute values
%             if incl_abs(v)==1
%                 eval([nme '_abs = abs(' nme '_av);']);
%                 if num_cond(v)>1
%                     for nc = 1:num_cond(v)
%                         eval([nme '_' num2str(nc) '_abs = abs(' nme '_' num2str(nc) ');']);
%                     end
%                 end
%             end
%             
%             %create tied ranks of data
%             if ranked == 1;
%                 % on averages
%                 eval(['vsize = size(' nme '_av);']);
%                 eval(['isort = floor(tiedrank(' nme '_av(:)));']);
%                 eval([nme '_av= reshape(isort,vsize(1),vsize(2));']);
%                 % on abs averaged data
%                 if incl_abs(v)==1
%                     eval(['vsize = size(' nme '_abs);']);
%                     eval(['isort = floor(tiedrank(' nme '_abs(:)));']);
%                     eval([nme '_abs= reshape(isort,vsize(1),vsize(2));']);
%                 end
%                 
%                 % on conditions
%                 if num_cond(v)>1
%                     if incl_abs(v)==1
%                         eval([nme '_cond_abs = abs(' nme ');']);
%                         eval(['vsize = size(' nme '_cond_abs);']);
%                         eval(['isort = floor(tiedrank(' nme '_cond_abs(:)));']);
%                         eval([nme '_cond_abs= reshape(isort,vsize(1),vsize(2));']);    
%                     end
%                     eval(['vsize = size(' nme ');']);
%                     eval(['isort = floor(tiedrank(' nme '(:)));']);
%                     eval([nme '= reshape(isort,vsize(1),vsize(2));']);   
%                     for nc = 1:num_cond(v) 
%                         if incl_abs(v)==1
%                             eval([nme '_' num2str(nc) '_abs= squeeze(nanmean(' nme '_cond_abs(:,nc),2));']);     
%                         end
%                         eval([nme '_' num2str(nc) '= squeeze(nanmean(' nme '(:,nc),2));']);
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% 