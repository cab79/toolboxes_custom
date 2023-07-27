function D = eegstats_dataprep_reformat(S,varargin)

%% find prepared data
if isempty(varargin)
   % try loading the data if not supplied in varargin
   try
       disp('loading data...');
       load(S.prep.path.inputs,'D') 
   catch
       error('Supply a D structure or path/file of the data')
   end
else
    D= varargin{1};
end

%% outputs: grouped data
switch S.prep.output.format

    case 'group'
        %% outputs: change subject to grouped data
        
        if length(D)==1
            error('data already grouped, or there is only one subject')
        end
        
        % create a single vertical concatenated design matrix
        dim=D(1).prep.dim;
        prep.dtab=D(1).prep.dtab;
        prep.grpdata{1} = D(1).prep.data;
        for d = 2:length(D)
            prep.dtab = vertcat(prep.dtab,D(d).prep.dtab);
            prep.grpdata{d} = D(d).prep.data;
        end
        dim(3)=height(prep.dtab);
        D=struct;
        D.prep=prep;
        D.prep.dim=dim;
        
        % vertically concatenate EEG data over subjects
        temp=prep.grpdata;
        D.prep = rmfield(D.prep,'grpdata');
        D.prep.grpdata{1} = vertcat(temp{:}); 
        clear grpdata temp
        
    case 'subject'
        %% outputs: change grouped to subject data
        
        if length(D)>1
            error('data already in subject format')
        end
        
        E=struct;
        [U,iu,iU]=unique(D.prep.dtab.ID,'stable');
        for d=1:length(U)
            E(d).prep.dtab = D.prep.dtab(ismember(D.prep.dtab.ID,U{d}),:);
            E(d).prep.dim = D.prep.dim;
            E(d).prep.dim(3) = height(E(d).prep.dtab);
            if isfield(D.prep,'Y')
                for y = 1:length(D.prep.Y)
                    E(d).prep.Y(y).dtab = D.prep.Y(y).dtab(ismember(D.prep.dtab.ID,U{d}),:);
                    E(d).prep.Y(y).data_mean = D.prep.Y(y).data_mean;
                    E(d).prep.Y(y).data_std = D.prep.Y(y).data_std;
                end
            end
        end
        D=E;
        clear E
        
    case 'Y'
        %% create Y structure using for encoding models
        for d=1:length(D)
            for ip = 1:length(D(d).prep)
                if iscell(D(d).prep(ip).grpdata)
                    D(d).prep(ip).grpdata = D(d).prep(ip).grpdata{:};
                end
                for s = 1:size(D(d).prep(ip).grpdata,2)
                    % duplicate dtab over data samples
            %         D(d).prep(ip).Y(s).dtab = D(d).prep(ip).dtab;
                    D(d).prep(ip).Y(s).dtab=table;

                    % add EEG after z-scoring over trials (unit variance) so that statistical coefficients are normalised
                    if S.prep.calc.eeg.zscore
                        D(d).prep(ip).Y(s).data_mean = nanmean(double(D(d).prep(ip).grpdata(:,s)));
                        D(d).prep(ip).Y(s).data_std = nanstd(double(D(d).prep(ip).grpdata(:,s)));
                        D(d).prep(ip).Y(s).dtab.data = (double(D(d).prep(ip).grpdata(:,s)) - D(d).prep(ip).Y(s).data_mean) / D(d).prep(ip).Y(s).data_std;
            %             [D(d).prep(ip).Y(s).dtab.data, D(d).prep(ip).Y(s).data_mean, D(d).prep(ip).Y(s).data_std] = zscore(double(D.prep.grpdata(:,s)));
                    else
                        D(d).prep(ip).Y(s).dtab.data = double(D(d).prep(ip).grpdata(:,s));
                    end
                end
            end
        end
        if isfield(D.prep,'grpdata') % don't need EEG data in both formats
            D.prep = rmfield(D.prep,'grpdata'); % save memory
        end
end

% save data to disk
if S.prep.output.save
    disp('saving data to disk...')
    save(fullfile(S.prep.path.outputs,[S.prep.sname '.mat']),'D','S','-v7.3')
    disp('...done')
end

