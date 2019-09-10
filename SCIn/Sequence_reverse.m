function Sequence_reverse(loadname)

[pth,fnme,ext]=fileparts(loadname);
if isempty(pth)
    pth = pwd;
end

load(fullfile(pth,[fnme ext]));

fields=fieldnames(seq);
for fn = 1:length(fields)
    seq.(fields{fn}) = fliplr(seq.(fields{fn}));
end

savename = fullfile(pth,[fnme '_reversed.' ext]);
save(savename,'seq','settings','filename');