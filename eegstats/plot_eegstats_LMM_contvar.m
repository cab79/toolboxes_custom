function plot_eegstats_LMM_contvar(S)

% plot continuous variables and their interactions, e.g. for interpreting LMM
% significant interaction effects. Plots are means/CIs over subjects for
% coefficient lines and means/SDs for descriptive plots.

% Based on: https://janhove.github.io/analysis/2017/06/26/continuous-interactions 

dbstop if error
cd(S.path.stats_load)

%% get variables 

% load model data
temp=load(fullfile(S.path.stats_load,S.file.stats_load),'D');
D=temp.D;
% load model data
temp=load(fullfile(S.path.dtab_inputs),'D');
Dp=temp.D;

% extract predictor variables data from D
pred={};
for m = 1:length(S.model.index)
    terms={D.model(S.model.index).con(:).term};
    termi = ismember(terms,S.model.contrast_term);
    pred{m} = D.model(S.model.index).fixeddesign(:,termi);
end

% extract outcome variables
for s = 1:length(S.model.datasample)
    outcome(:,s) = Dp.prep.Y(S.model.datasample(s)).dtab.data;
end

% calculate means over subjects, per trial
tnums = unique(Dp.prep.dtab.trial);
[U,~,iU] = unique(Dp.prep.dtab.ID);
for m = 1:length(S.model.index)
    predmean{m} = nan(max(tnums),length(S.model.contrast_term),length(U));
    for u = 1:length(U)
        predmean{m}(Dp.prep.dtab.trial(iU==u),:,u) = pred{m}(iU==u,:);
    end
    predmean{m} = nanmean(predmean{m},3);
end
outmean = nan(max(tnums),length(S.model.datasample),length(U));
for u = 1:length(U)
    outmean(Dp.prep.dtab.trial(iU==u),:,u) = outcome(iU==u,:);
end
outmean = nanmean(outmean,3);

%% PLOTS: 1 figure per model, rows per data sample, columns per pred/out combination

% exploratory lowess plots of each predictor combination and each predictor vs. outcome
for m = 1:length(S.model.index)
    figure
    t=tiledlayout(length(S.model.datasample),3);
    title(t,'exploratory plots')
    for s = 1:length(S.model.datasample)
        % pred1 v outcome
        nexttile; 
        fit_1D(predmean{m}(:,1),outmean(:,s));
        ax(1)=gca;
        xlabel(S.model.contrast_term{1})
        ylabel(['data sample ' num2str(s)])
        % pred2 v outcome
        nexttile; 
        fit_1D(predmean{m}(:,2),outmean(:,s));
        ax(2)=gca;
        xlabel(S.model.contrast_term{2})
        ylabel(['data sample ' num2str(s)])
        % pred1 v pred2
        nexttile; 
        fit_1D(predmean{m}(:,1),predmean{m}(:,2));
        ax(3)=gca;
        xlabel(S.model.contrast_term{1})
        ylabel(S.model.contrast_term{2})
        linkaxes([ax(1) ax(2)],'y')
    end
end

% conditioning plots of each variable vs. outcome
for m = 1:length(S.model.index)
    figure
    ax=[];
    t=tiledlayout(length(S.model.datasample),S.model.cond_nranges);
    title(t,['conditioning plots on ' S.model.cond_var])
    for s = 1:length(S.model.datasample)

        % conditioning and plotting variables
        vc = find(ismember(S.model.contrast_term,S.model.cond_var)); % conditioning
        vp = find(~ismember(S.model.contrast_term,S.model.cond_var)); % plotting
        condvar = predmean{m}(:,vc); 
        plotvar = predmean{m}(:,vp); 
        
        % overlapping quantile ranges of the conditioning variable
        condval = [min(condvar) quantile(condvar,S.model.cond_nranges) max(condvar)];
        
        % data indices for each range
        for r = 1:S.model.cond_nranges
            cvals = condval([r,r+2]);
            range = find(condvar>=cvals(1) & condvar<=cvals(2));
            
            % plotvar v outcome for each range of conditioning variable
            nexttile; 
            fit_1D(plotvar(range),outmean(range,s));
            ax(r)=gca;
            xlabel([S.model.contrast_term{1} ', range ' num2str(r)])
            ylabel(['data sample ' num2str(s)])
        end
        linkaxes(ax)

        
    end
    
end
% coefficients: each predictor vs. outcome

% coefficients: three-dimensional plot of interaction

% coefficients: contour plot of interaction

function fit_1D(x,z,S)
rmnan = isnan(x) | isnan(z);
x(rmnan)=[];z(rmnan)=[];
%span = span_lowess_CV(x,z);
span = 0.99; 
xx = linspace(min(x),max(x),length(x));

hold on
scatter(x,z,'k');
ll=lsline;
ll.Color = 'b';
ll.LineWidth = 2;
line(xx,lowessfit_1D([x,z],xx,span),'color','r','linestyle','-', 'linewidth',2)
hold off

function lowess_2D(x,y,z)
f = fit([x y],z,'lowess');
plot(f,[x y],z);

function span = span_lowess_CV(x,z)
% https://blogs.mathworks.com/loren/2011/01/13/data-driven-fitting/
disp('calculating cross-validated lowess span')
num = 9;
spans = linspace(.1,.9,num);
sse = zeros(size(spans));
cp = cvpartition(length(z),'k',10);

for j=1:length(spans)
    disp(['span ' num2str(j) '/' num2str(length(spans))])
    f = @(train,test) norm(test(:,2) - lowessfit_1D(train,test(:,1),spans(j)))^2;
    sse(j) = sum(crossval(f,[x,z],'partition',cp));
end

[minsse,minj] = min(sse);
span = spans(minj);

function ys=lowessfit_1D(xy,xs,span)
% https://blogs.mathworks.com/loren/2011/01/13/data-driven-fitting/
%MYLOWESS Lowess smoothing, preserving x values
%   YS=MYLOWESS(XY,XS) returns the smoothed version of the x/y data in the
%   two-column matrix XY, but evaluates the smooth at XS and returns the
%   smoothed values in YS.  Any values outside the range of XY are taken to
%   be equal to the closest values.

if nargin<3 || isempty(span)
    span = .3;
end

% Sort and get smoothed version of xy data
xy = sortrows(xy);
x1 = xy(:,1);
y1 = xy(:,2);
ys1 = smooth(x1,y1,span,'rloess');

% Remove repeats so we can interpolate
t = diff(x1)==0;
x1(t)=[]; ys1(t) = [];

% Interpolate to evaluate this at the xs values
ys = interp1(x1,ys1,xs,'linear',NaN);

% Some of the original points may have x values outside the range of the
% resampled data.  Those are now NaN because we could not interpolate them.
% Replace NaN by the closest smoothed value.  This amounts to extending the
% smooth curve using a horizontal line.
if any(isnan(ys))
    ys(xs<x1(1)) = ys1(1);
    ys(xs>x1(end)) = ys1(end);
end