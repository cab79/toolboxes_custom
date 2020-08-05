function [coeff,Fcoeff,Rcoeff] = LMM_FReffects(lme,effectstr,rdm_levels,IVs)
% EXTRACT COEFFICIENTS

effects=strsplit(effectstr,':');
% get effects
[F,Fn]=fixedEffects(lme);
Fn=Fn.Name;
[R,Rn]=randomEffects(lme);
Rnames = Rn.Level(strcmp(Rn.Name,'(Intercept)'));
[~,Ridx] = ismember(rdm_levels, Rnames);
Rn=Rn.Name;
% is group a factor?
if any(strcmp(IVs,'group'))
    grprange=lme.VariableInfo('group','Range');
    grps=grprange{:,:};
    grps=grps{:};
    gcell=[{''} grps(2:end)];
else
    gcell={''};
    grps={''};
end
% build within-subject terms
wIVs = IVs; 
wIVs(strcmp(wIVs,'group'))=[];
effect_is_IV = find(ismember(wIVs,effects));
if ~isempty(effect_is_IV)
    wIVs(effect_is_IV)=[];
    effectstr = [effectstr '_1'];
end
if ~isempty(wIVs)
    [~, celltypes] = get_cells(lme.Variables, lme.Variables.Properties.VariableNames, wIVs);
    terms={};term={};
    for c = 1:height(celltypes)
        for w = 1:length(wIVs)
            if double(celltypes.(wIVs{w})(c))>1
                terms{c,w} = [wIVs{w} '_' num2str(double(celltypes.(wIVs{w})(c))-1) ':'];
            else
                terms{c,w} = '';
            end
        end
        term{c,1}=horzcat(terms{c,:});
    end
    term = regexprep(term,'^:','');
    term = regexprep(term,':$','');
else
    term={''};
end

% initialise tables
Rcoeff=table;
Fcoeff=table;
% for each group
for g = 1:length(gcell)

    % get terms for this group
    if g==1
        eff_exp = '^((?!group).)*$'; % must not contain group term
    else
        eff_exp = ['.*group_' gcell{g} '.*']; % must contain group term
    end
    eff_fun = @(s)~cellfun('isempty',regexp(Fn,eff_exp));
    eff_out = cellfun(eff_fun,effects,'UniformOutput',false);
    eff_idx = all(horzcat(eff_out{:}),2);
    Fg=F(eff_idx);
    Fng=Fn(eff_idx);
    
    % remove group term to make later processing simpler
    Fng= regexprep(Fng,['(:?)group_' gcell{g} '(:?)'],'');
    
    % current version: random and fixed effects must correspond (i.e. full random slopes and intercepts model)
    if g==1
        
        if ~isequal(Fng,unique(Rn,'stable'))
%             error('random and fixed effect terms do not correspond')
            
        end
    end
    
    % get term indices and intercepts
    if g==1
        Fng{strcmp(Fng,'(Intercept)')}='';
    end
    for t = 1:size(term,1)
        term_idx(t) = find(strcmp(Fng,term(t)));
    end
    
    % identify the slopes of interest
%     eff_exp = '.*'; % expr = no colon before the first effect term
%     for eff = 1:length(effects)
%         eff_exp = [eff_exp effects{eff} '.*']; % expr = no colon after last effect term
%     end
%     eff_fun = @(s)~cellfun('isempty',regexp(Fng,eff_exp));
%     eff_out = cellfun(eff_fun,effects,'UniformOutput',false);
%     eff_idx = find(all(horzcat(eff_out{:}),2));
%     eff_terms = Fng(eff_idx);
    
    % build effect terms)
    eff_terms={};
    slopes=1;
    for t = 1:size(term,1)
        eff_terms{t,1} = [term{t} ':' effectstr];
        eff_terms{t,1} = regexprep(eff_terms{t,1},'^:','');
        if any(strcmp(Fng,eff_terms(t)))
            slp_idx(t) = find(strcmp(Fng,eff_terms(t)));
        else
            slopes = 0;
        end
    end
   
    % fixed effects
    Fx_int = [Fg(term_idx(1)), Fg(term_idx(2:end))'+Fg(term_idx(1))];
    Fcoeff.intercept(g,:)=Fx_int;
    if slopes
        Fx_slp = [Fg(slp_idx(1)), Fg(slp_idx(2:end))'+Fg(slp_idx(1))];
        Fcoeff.slope(g,:)=Fx_slp;
    else
        Fcoeff.slope(g,:)=zeros(size(Fcoeff.intercept(g,:)));
    end
    
    % if g>1 add to first group's (reference) data
    if g>1
%         Fcoeff{g,2:end} = Fcoeff{g,2:end} + Fcoeff{1,2:end}; % don't include intercept
        Fcoeff.intercept(g,:)=Fcoeff.intercept(g,:)+Fcoeff.intercept(1,:);
        Fcoeff.slope(g,:)=Fcoeff.slope(g,:)+Fcoeff.slope(1,:);
    end
end
% random effects
Rn(strcmp(Rn,'(Intercept)'))={''};
for t = 1:length(term)
    if t==1
        Rcoeff.intercept(:,t)=R(strcmp(Rn,term{t}));
        if slopes
            Rcoeff.slope(:,t)=R(strcmp(Rn,eff_terms{t}));
        end
    else
        if slopes
            Rcoeff.intercept(:,t)=R(strcmp(Rn,term{t})) + R(strcmp(Rn,term{1}));
            Rcoeff.slope(:,t)=R(strcmp(Rn,eff_terms{t})) + R(strcmp(Rn,eff_terms{1}));
        end
    end
end
if ~slopes
    Rcoeff.slope=zeros(size(Rcoeff.intercept));
end
Rcoeff=Rcoeff(Ridx,:);
% create coeff table
coeff=table;
coeff.rdm_levels=rdm_levels;
for r = 1:height(Rcoeff)
    coeff.group(r) = unique(lme.Variables.group(strcmp(lme.Variables.ID,rdm_levels(r))));
end
Rcoeff = [coeff,Rcoeff];
coeff=Rcoeff;
% add in fixed effects
if slopes
    for r = 1:height(coeff)
        gi = 1;
        if any(strcmp(coeff.group(r),gcell))
            gi=find(strcmp(coeff.group(r),gcell));
        end
        coeff{r,3:end} = coeff{r,3:end} + Fcoeff{gi,:};
    end
end
gtab=table; gtab.group=grps';
Fcoeff=[gtab Fcoeff];


function [cells,celltypes] = get_cells(design,headernames,IVname)
for iv = 1:length(IVname)
    header_col(iv) = find(ismember(headernames,IVname(iv)));
end
temp=design(:,header_col);
[celltypes,~,cells] = unique(temp,'rows','stable');


