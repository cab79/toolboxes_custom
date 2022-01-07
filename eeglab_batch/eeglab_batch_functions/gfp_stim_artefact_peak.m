function S=gfp_stim_artefact_peak(S)

dbstop if error
S.func = 'gfp';
if isfield(S.(S.func),'GFP')
    S.(S.func) = rmfield(S.(S.func),'GFP');
end

% select channels
S=select_chans(S);

switch S.(S.func).select.datatype
    case 'ERP'
        S.path.file = S.path.erp;
    case 'TF'
        S.path.file = S.path.freq;
    case 'Freq'
        error('GFP analysis not possible with Freq data, only ERP and TF data')
end

S = getfilelist(S);

figure; hold on
for f = 1:length(S.(S.func).filelist)
    file = S.(S.func).filelist{f};
    load(fullfile(S.path.file,file));
    % smooth gavg for peak identification
    for e = 1:size(tldata{1}.avg,1)
        tldata{1}.avg(e,:) = smooth(squeeze(tldata{1}.avg(e,:)),5,'lowess');
    end
    GFP(f,:) = std(tldata{1}.avg,0,1);
    plot(tldata{1}.time,GFP(f,:));
end
GFP_mean = mean(GFP,1);
plot(tldata{1}.time,GFP_mean,'LineWidth',3);
[peak,loc] = max(GFP_mean)
tldata{1}.time(loc(1))

