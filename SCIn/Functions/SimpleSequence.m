function h = SimpleSequence(h)

disp('Creating sequence...');
h.Seq =struct;


%% create final sequences/blocks
if ~isfield(h.Seq,'signal')
    if length(h.Settings.stim.dur) > 1
        h.Seq.signal=1:h.Settings.ntrials;
        h.Seq.condnum=1:h.Settings.ntrials;
        h.Seq.blocks=ones(1,h.Settings.ntrials);
    else
        h.Seq.signal=ones(1,h.Settings.ntrials);
        h.Seq.condnum=ones(1,h.Settings.ntrials);
        h.Seq.blocks=ones(1,h.Settings.ntrials);
    end
end