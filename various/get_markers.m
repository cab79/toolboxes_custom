
function [conds, tnums, fnums, bnums] = get_markers(EEG,type)

conds = nan(1,length(EEG.epoch));
tnums = nan(1,length(EEG.epoch));
fnums = nan(1,length(EEG.epoch));
bnums = nan(1,length(EEG.epoch));

switch type
    case 'EGI'
        for ep = 1:length(EEG.epoch)
            stimevidx = find(strcmp('STIM',EEG.epoch(ep).eventtype));
            if ep<length(EEG.epoch); stimevidx1 = find(strcmp('STIM',EEG.epoch(ep+1).eventtype));end
            if ~isempty(stimevidx)
                stimcodes = EEG.epoch(ep).eventcodes{stimevidx(end)};
                if ~any(strcmp('CNUM',stimcodes(:,1)))
                    error('change CNUM to FNUM to analyse conditions');
                end
                conds(1,ep) = stimcodes{strcmp('CNUM',stimcodes(:,1)),2};
                tnums(1,ep) = stimcodes{strcmp('TNUM',stimcodes(:,1)),2};
                fnums(1,ep) = stimcodes{strcmp('FNUM',stimcodes(:,1)),2};
                bnums(1,ep) = stimcodes{strcmp('BNUM',stimcodes(:,1)),2};

            else
                if length(stimevidx1)==2
                    stimcodes = EEG.epoch(ep+1).eventcodes{stimevidx1(1)};
                    conds(1,ep) = stimcodes{strcmp('CNUM',stimcodes(:,1)),2};
                    tnums(1,ep) = stimcodes{strcmp('TNUM',stimcodes(:,1)),2};
                    fnums(1,ep) = stimcodes{strcmp('FNUM',stimcodes(:,1)),2};
                    bnums(1,ep) = stimcodes{strcmp('BNUM',stimcodes(:,1)),2};
                else
                    error(['too many / too few STIMs on trial ' num2str(ep+1)])
                end
            end
        end
    case 'BV'
        evtypes = unique(horzcat(EEG.epoch(:).eventtype));
        for ep = 1:length(EEG.epoch)
            stimevidx = find([EEG.epoch(ep).eventlatency{:}]==0);
            if ~isempty(stimevidx)
                str=EEG.epoch(ep).eventtype{stimevidx};
                conds(1,ep) = str2double(regexp(str,'.* (\d+)','tokens','once'));
                tnums(1,ep) = EEG.epoch(ep).eventurevent{stimevidx};
            end
        end
end