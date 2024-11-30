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

        % select subjects
        if isfield(S.prep.output,'select') && ~isempty(S.prep.output.select)
            idx = D.prep.dtab.(S.prep.output.select)==1;
            D.prep.dtab = D.prep.dtab(idx,:);
            D.prep.grpdata{1} = D.prep.grpdata{1}(idx,:);
            D.prep.dim(3) = length(D.prep.grpdata{1});
            D.prep.(S.prep.output.select) = idx;
        end
        
    case 'subject'
        %% outputs: change grouped to subject data
        
        if length(D)>1
            error('data already in subject format')
        end

        % E=struct;
        % [U,iu,iU]=unique(D.prep.dtab.ID,'stable');
        % for d=1:length(U)
        %     tic
        %     E(d).prep.dtab = D.prep.dtab(ismember(D.prep.dtab.ID,U{d}),:);
        %     E(d).prep.dim = D.prep.dim;
        %     E(d).prep.dim(3) = height(E(d).prep.dtab);
        %     if isfield(D.prep,'Y')
        %         for y = 1:length(D.prep.Y)
        %             E(d).prep.Y(y).dtab = D.prep.Y(y).dtab(ismember(D.prep.dtab.ID,U{d}),:);
        %             if isfield(D.prep.Y(y),'data_mean')
        %                 E(d).prep.Y(y).data_mean = D.prep.Y(y).data_mean;
        %                 E(d).prep.Y(y).data_std = D.prep.Y(y).data_std;
        %             end
        %         end
        %     end
        %     toc
        % end
        % D=E;
        % clear E
        
        % Extract unique IDs and their indices
        [U, ~, iU] = unique(D.prep.dtab.ID, 'stable');
        nUnique = length(U);
        
        % Preallocate the struct array E to improve speed
        E(nUnique).prep = struct();
        
        % Precompute dimensions and set flags for fields
        hasFieldY = isfield(D.prep, 'Y');
        if hasFieldY
            nY = length(D.prep.Y);
        end
        
        % Create a logical index matrix in advance to avoid repeated `ismember` calls
        ID_indices = arrayfun(@(x) iU == x, 1:nUnique, 'UniformOutput', false);
        
        % Extract frequently used deeply nested fields to local variables to reduce access time
        prep_dtab = D.prep.dtab;
        prep_dim = D.prep.dim;
        if hasFieldY
            prep_Y = D.prep.Y;
        end
        
        % Loop through unique IDs only once
        for d = 1:nUnique
            tic
            idx = ID_indices{d};  % Logical index for current ID
            E(d).prep.dtab = prep_dtab(idx, :);  % Get the rows with matching ID
            E(d).prep.dim = prep_dim;
            E(d).prep.dim(3) = sum(idx);  % More efficient than height()
        
            % Only proceed if field 'Y' exists in D.prep
            if hasFieldY
                % Preallocate Y struct for efficiency
                Eprep_Y = struct('dtab', [], 'data_mean', [], 'data_std', []);  % Create a local variable for E(d).prep.Y
                for y = 1:nY
                    % Extract the rows for dtab based on the logical index
                    Eprep_Y(y).dtab = prep_Y(y).dtab(idx, :);
                    % Copy fields data_mean and data_std if they exist
                    if isfield(prep_Y(y), 'data_mean')
                        Eprep_Y(y).data_mean = prep_Y(y).data_mean;
                        Eprep_Y(y).data_std = prep_Y(y).data_std;
                    end
                end
                % Assign the local variable to E(d).prep.Y
                E(d).prep.Y = Eprep_Y;
            end
            toc
        end
        
        % Replace the original structure D with the new struct array E
        D = E;
        clear E

        
    case 'Y'
        %% create Y structure using for encoding models
        for d=1:length(D)
            for ip = 1:length(D(d).prep)
                if ~isfield(D(d).prep(ip),'grpdata')
                    D(d).prep(ip).grpdata = D(d).prep(ip).data;
                end
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
            if isfield(D(d).prep,'grpdata') % don't need EEG data in both formats
                D(d).prep = rmfield(D(d).prep,'grpdata'); % save memory
            end
            if isfield(D(d).prep,'data') % don't need EEG data in both formats
                D(d).prep = rmfield(D(d).prep,'data'); % save memory
            end
        end
        
    case 'grpdata'
    %% Convert from Y structure back to grpdata matrix

    for d = 1:length(D)
        % Check if D(d).prep.Y exists
        if isfield(D(d).prep, 'Y')
            numY = length(D(d).prep.Y);

            % Initialize grpdata matrix
            numTrials = height(D(d).prep.dtab);
            grpdata = zeros(numTrials, numY);

            for s = 1:numY
                % Get data from Y(s).dtab.data
                grpdata(:, s) = D(d).prep.Y(s).dtab.data;

                % If data was z-scored, reverse the z-scoring
                if isfield(D(d).prep.Y(s), 'data_mean') && isfield(D(d).prep.Y(s), 'data_std')
                    grpdata(:, s) = grpdata(:, s) * D(d).prep.Y(s).data_std + D(d).prep.Y(s).data_mean;
                end
            end

            % Assign grpdata to D(d).prep.grpdata
            D(d).prep.grpdata = grpdata;

            % Optionally, remove Y structure to save memory
            if isfield(S.prep.output, 'removeY') && S.prep.output.removeY
                D(d).prep = rmfield(D(d).prep, 'Y');
            end
        else
            error('Y structure not found in D(%d).prep', d);
        end
    end

    case 'groupY'
        %% outputs: change subject Y data to grouped Y data

        if length(D)==1
            error('Data already grouped, or there is only one subject')
        end

        % Initialize new D structure
        prep = struct;

        % Concatenate dtab
        prep.dtab = D(1).prep.dtab;
        for d = 2:length(D)
            prep.dtab = vertcat(prep.dtab, D(d).prep.dtab);
        end

        % Set dim
        dim = D(1).prep.dim;
        dim(3) = height(prep.dtab);

        % Number of Y elements (samples)
        numY = length(D(1).prep.Y);

        % For each Y (sample)
        for s = 1:numY
            % Initialize D.prep.Y(s).dtab
            prep.Y(s).dtab = table;

            % Initialize data variables
            data = [];

            % For each subject
            for d = 1:length(D)
                % Append data
                data = [data; D(d).prep.Y(s).dtab.data];
            end

            % Assign concatenated data
            prep.Y(s).dtab.data = data;

            % Handle data_mean and data_std if they exist
            if isfield(D(1).prep.Y(s), 'data_mean') && isfield(D(1).prep.Y(s), 'data_std')
                % Recompute data_mean and data_std from the concatenated data
                prep.Y(s).data_mean = nanmean(data);
                prep.Y(s).data_std = nanstd(data);
            end
        end

        D = struct;
        D.prep = prep;
        D.prep.dim = dim;

        % select subjects if specified
        if isfield(S.prep.output,'select') && ~isempty(S.prep.output.select)
            idx = D.prep.dtab.(S.prep.output.select)==1;
            D.prep.dtab = D.prep.dtab(idx,:);
            for s = 1:numY
                D.prep.Y(s).dtab.data = D.prep.Y(s).dtab.data(idx,:);
            end
            D.prep.dim(3) = height(D.prep.dtab);
            D.prep.(S.prep.output.select) = idx;
        end

    case 'swapY'
        %% swap over Y variables with specified dtab variables

        % Get the number of subjects
        numSubjects = numel(D);
        
        % Loop through each subject
        for i = 1:numSubjects
            
            % Loop through each Y field instance
            for j = 1:numel(D(i).prep.Y)
                
                % Get the data from Y(y).dtab and add it to dtab
                if isfield(D(i).prep.Y(j), 'dtab')
                    D(i).prep.dtab.(S.prep.dtabHeadersToAdd{j}) = D(i).prep.Y(j).dtab.data;
                end
            end

            % Remove data from Y(y).dtab
            D(i).prep.Y = rmfield(D(i).prep.Y,'dtab');
            
            % Remove specified headers from dtab and place them within Y structure
            for m = 1:numel(S.prep.dtabHeadersToRemove)
                headerToRemove = S.prep.dtabHeadersToRemove{m};

                if contains(headerToRemove, D(i).prep.dtab.Properties.VariableNames)

                    % Add data to Y(y).dtab
                    D(i).prep.Y(m).dtab = table;
                    D(i).prep.Y(m).dtab.data = D(i).prep.dtab.(headerToRemove);
                    
                    % Remove data from dtab
                    D(i).prep.dtab = removevars(D(i).prep.dtab, headerToRemove);
                end
            end
            D(i).prep.dim(2) = numel(S.prep.dtabHeadersToRemove);
        end




end

% save data to disk
if S.prep.output.save
    disp('saving data to disk...')
    save(fullfile(S.prep.path.outputs,[S.prep.sname '.mat']),'D','S','-v7.3')
    disp('...done')
end

