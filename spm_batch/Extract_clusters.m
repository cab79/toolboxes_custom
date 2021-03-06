function Extract_clusters(S)

%% Requires input structure S containing fields as follows. See Cluster_processing script for examples.
%-------------------------------------------------------------
% S.spmstats_path: directory in which SPM analysis is saved 
% S.spm_dir: specific folder containing the SPM stats for this analysis
% S.contrasts: contrast name cell array - must match that in Matlabbatch (i.e. from design-batch script). Leave empty to proccess ALL contrasts in Matlabbatch
% S.clustab: stats headers to save in a table 
% S.batch: name of batch .mat file saved from design_batch.m and within same folder as SPM.mat

%OPTIONAL SETTINGS (advise not to change):
%Num=16;    %- number of maxima per cluster [3]
%Dis=8;    %- min distance among clusters {mm} [8]
%Str='';    %- header string
dbstop if error
%% RUN
if isempty(S.spm_dir)
    S.spm_paths = dir(fullfile(S.spmstats_path,'*spm*'));
    S.spm_paths = {S.spm_paths(:).name};
elseif ~iscell(S.spm_dir)
    S.spm_paths = {S.spm_dir};
else
    S.spm_paths = S.spm_dir;
end

for sp = 1:length(S.spm_paths)
    S.spm_path = fullfile(S.spmstats_path,S.spm_paths{sp});

    if isempty(S.contrasts)
        % load SPM design and list contrasts
        load(fullfile(S.spm_path,S.batch));
        S.contrasts = {};
        tf=[];
        for fc = 1:length(matlabbatch{3}.spm.stats.con.consess)
            try
                S.contrasts{fc} = matlabbatch{3}.spm.stats.con.consess{1,fc}.fcon.name;
                tf=[tf 1];
            catch
                S.contrasts{fc} = matlabbatch{3}.spm.stats.con.consess{1,fc}.tcon.name;
                tf=[tf 2];
            end
        end
    else
        tf=S.tf;
    end

    spm eeg
    S.clus_path = {};
    for cont = 1:length(S.contrasts)
        contrast = S.contrasts{cont};
        Nhead = size(S.clustab{tf(cont)},1);
        Nrhead=1;

        load(fullfile(S.spm_path,'SPM.mat'));
        allcon = {SPM.xCon.name};
        con = find(strcmp(allcon,contrast));
        
        % fudge to prevent a crash if two contrasts are named the same:
        if length(con)>1
            con = con(con>prev_con);
            con=con(1);
        end
        prev_con = con;

        % settings
        SPMset.swd=S.spm_path;
        SPMset.thresDesc = S.thresDesc; % no correction
        SPMset.u = S.clusformthresh; % uncorrected threshold
        SPMset.k = 0; % extent threshold
        SPMset.units = {'mm' 'mm' 'ms'};
        SPMset.Ic = con; % indices of contrasts (in SPM.xCon)
        SPMset.title = contrast;
        SPMset.n=1; %    - conjunction number <= number of contrasts
        if ~isempty(S.imgmask)
            SPMset.Im={[S.imgmask ',1']}; % masking
        else
            SPMset.Im=[]; % no masking
        end
        SPMset.pm=[];      % - p-value for masking (uncorrected)
        SPMset.Ex=[];%       - flag for exclusive or inclusive masking

        % load results and generate initial stats
        s_merged = rmfield(SPM, intersect(fieldnames(SPM), fieldnames(SPMset)));
        names = [fieldnames(s_merged); fieldnames(SPMset)];
        SPM = cell2struct([struct2cell(s_merged); struct2cell(SPMset)], names, 1);
        [hReg,xSPM,SPM] = spm_results_ui('setup',SPM);
        TabDat = spm_list('List',xSPM,hReg);

        % find column indices of table values to extract
        A=TabDat.hdr(1:2,:)';
        B=S.clustab{tf(cont)}';
        C = union(A,B) ;
        [dum,ai] = ismember(A,C) ;
        [dum,bi] = ismember(B,C) ;
        [tfn,loc] = ismember(ai,bi,'rows');
        colind = find(tfn);

        % add new column to save temporal extent of cluster
        clustable=S.clustab{tf(cont)};
        colwidth = size(clustable,2);
        clustable(2,colwidth+1)={'Min temporal extent (ms)'};
        clustable(2,colwidth+2)={'Max temporal extent (ms)'};

        % re-draw stats using cluster extent threshold
        if S.clusformthresh && ~S.clusformthresh==1
            SPMset.k=xSPM.uc(3)-1;
            s_merged = rmfield(SPM, intersect(fieldnames(SPM), fieldnames(SPMset)));
            names = [fieldnames(s_merged); fieldnames(SPMset)];
            SPM = cell2struct([struct2cell(s_merged); struct2cell(SPMset)], names, 1);
            [hReg,xSPM,SPM] = spm_results_ui('setup',SPM);
            TabDat = spm_list('List',xSPM,hReg);
        end

        % move on to the next contrast if there are no significant results
        if isempty(TabDat.dat)
            continue
        end

        % create folder to save cluster data
        conname = strrep(contrast,' ','');
        conname = strrep(conname,'*','_');
        S.clus_path{cont} = fullfile(S.spm_path,[conname '_clusters']); 
        if exist(S.clus_path{cont},'dir')
            S.clus_path{cont} = [S.clus_path{cont} '_new'];
        end
        if ~exist(S.clus_path{cont},'dir')
            mkdir(S.clus_path{cont});
        end

        % select each cluster in turn
        TabDatR = TabDat.dat(~cellfun('isempty',TabDat.dat(:,3)),:);
        nClus = TabDat.dat{1,2};
        for c=1:nClus
            cname=['c' num2str(c)];
            xyz=TabDatR{c,12};
            spm_results_ui('SetCoords',xyz);
            TabDat = spm_list('ListCluster',xSPM,hReg,100,1);
            clustable(Nhead+c,Nrhead)={cname};
            
            if isempty(TabDat.dat)
                continue
            end
            % extract table data
            clustable(Nhead+c,Nrhead+(1:length(colind)))=TabDat.dat(1,colind);
            clustable(Nhead+c,Nrhead+7:Nrhead+9) = num2cell(round([clustable{Nhead+c,Nrhead+6}])');
            % vec2str: https://uk.mathworks.com/matlabcentral/fileexchange/22937
            clustable{Nhead+c,Nrhead+6} = vec2str(round([clustable{Nhead+c,Nrhead+6}])',{},{},0);
            loc = [TabDat.dat{:,colind(end)}];
            locF = [TabDat.dat{:,9}];
            clustable(Nhead+c,colwidth+1)={min(loc(3,:))};
            clustable(Nhead+c,colwidth+2)={max(loc(3,:))};

            % save cluster volume image for masking
            cfname = fullfile(S.clus_path{cont},[cname '_mask']);
            mysavespm('current',cfname);  % mysavespm is a function from spm_results_ui

            % re-draw projection with cluster mask
            SPMset.Im={[cfname '.nii,1']}; % masking
            SPMset.pm=[];      % - p-value for masking (uncorrected)
            SPMset.Ex=[];%       - flag for exclusive or inclusive masking
            s_merged = rmfield(SPM, intersect(fieldnames(SPM), fieldnames(SPMset)));
            names = [fieldnames(s_merged); fieldnames(SPMset)];
            SPM = cell2struct([struct2cell(s_merged); struct2cell(SPMset)], names, 1);
            [hReg,xSPM,SPM] = spm_results_ui('setup',SPM);
            spm_results_ui('SetCoords',xyz);
            TabDat = spm_list('ListCluster',xSPM,hReg);

            % save thresholded SPM of the cluster
            cfname = fullfile(S.clus_path{cont},[cname '_spm']);
            mysavespm('thresh',cfname);  % mysavespm is a function from spm_results_ui

            %save figure 
            print(cfname,'-dpng') 

            % save VOI
            xY=struct;
            xY.xyz=xyz;
            xY.name=cname;
            xY.Ic=0;           %- contrast used to adjust data (0 - no adjustment)
            xY.def = 'cluster';%- VOI definition
            spm_regions_saveoptions(xSPM,SPM,hReg,xY,0); % modified version of spm_regions that can prevent saving xY
            pause(1)
            movefile(fullfile(S.spm_path,['VOI_' cname '.mat']),S.clus_path{cont});

            % unmask stats for the next cluster
            SPMset.k=xSPM.uc(3)-1;
            SPMset.Im=[]; % no masking
            s_merged = rmfield(SPM, intersect(fieldnames(SPM), fieldnames(SPMset)));
            names = [fieldnames(s_merged); fieldnames(SPMset)];
            SPM = cell2struct([struct2cell(s_merged); struct2cell(SPMset)], names, 1);
            [hReg,xSPM,SPM] = spm_results_ui('setup',SPM);
            TabDat = spm_list('List',xSPM,hReg);
        end

        save(fullfile(S.clus_path{cont},'cluster_table.mat'),'clustable');
        xlswrite(fullfile(S.clus_path{cont},'cluster_table.xlsx'),clustable);

        % save 
        %save(fullfile(S.clus_path{cont},'cluster_data.mat'),'S');
    end
end
