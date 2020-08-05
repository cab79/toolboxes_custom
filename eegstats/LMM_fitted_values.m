function [fitM_val,fitC_val] = LMM_fitted_values(lme,fitM,fitC)
% EXTRACT FITTED VALUES PER CELL/SUBJECT FOR PLOTTING

% define IVs
IVs = lme.PredictorNames;
catV=[];
contV=[];
for v = 1:length(IVs)
    dat=lme.Variables.(IVs{v});
    if strcmp(IVs{v},'ID')
        idV = v;
    elseif iscategorical(dat) || iscell(dat)
        catV=[catV v];
    else
        contV=[contV v];
    end
end
% obtain marginal fitted values per cell
[cells, fitM_val] = IV_cells_labels(lme.Variables,IVs(catV),[]);
[u,~,irows] = unique(cells);
for i = u'
    Ufitted(i,1)=mean(fitM(irows==i));
end
fitM_val.fitted=Ufitted;
% obtain conditional fitted values per cell/subject
[cells, fitC_val] = IV_cells_labels(lme.Variables,IVs([idV catV]),[]);
[u,~,irows] = unique(cells);
for i = u'
    Cfitted(i,1)=mean(fitC(irows==i));
end
fitC_val.fitted=Cfitted;
