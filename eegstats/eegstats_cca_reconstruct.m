function eegstats_cca_reconstruct(comp)
dbstop if error
close all
pth='C:\Data\CORE\eeg\ana\prep\cleaned\part2\eegstats_dataprep';
fname = 'eegstats_dtab20201002T105819.mat';
cd(pth);

load(fullfile(pth,fname))
cca_comp_select = [comp];
centre_output = 0;
d=1;
[U,~,iU]=unique(D(d).prep.dtab.ID,'stable');
dat={};

% e.g.
% 12 temporal
% 12 x 5 spatial
% 17 both
for pt = length(D(d).prep.PCA):-1:1
    
    % unpack
    pcatype = D(d).prep.PCA(pt).type;
    O = D(d).prep.PCA(pt).O;
    oind = D(d).prep.PCA(pt).oind;
    currdim = D(d).prep.PCA(pt).currdim;
    options = D(d).prep.PCA(pt).options;
    nCCA=size(O(1).W,2);
    
    if pt==length(D(d).prep.PCA)
        CCAs = D(d).prep.PCA(pt).CCA;
        if isempty(cca_comp_select)
            compi=1:nCCA;
        else
            compi=cca_comp_select;
        end
    else
        CCAs = dat;
        compi=1:nCCA;
    end

    temp={};
    for u=1:length(U)
        
        dn = ndims(CCAs{u});
        
        % for multiple observations only
        if length(O)>1
            if dn==2
                sz=size(CCAs{u});
                CCAs{u} = permute(reshape(CCAs{u},sz(1),nCCA,[]),[1 3 2]);
            elseif dn==3 
                if strcmp(pcatype,'spatial')
                    CCAs{u} = permute(CCAs{u},[1 3 2]); % trials x chan x time
                elseif strcmp(pcatype,'temporal')
                elseif strcmp(pcatype,'both')
                end
                sz=size(CCAs{u});
%                 CCAs{u} = permute(reshape(CCAs{u},sz(1),sz(2),nCCA,[]),[1 2 4 3]);
            end
        end
        
        for o = 1:length(O)
            
            dn = ndims(CCAs{u});
            switch dn
                case 2
                    CCA = CCAs{u}(:,compi);
                case 3
                    if length(O)>1
                        CCA = reshape(CCAs{u}(:,o,compi),[],length(compi));
                    else
                        if strcmp(pcatype,'spatial')
                            CCA = reshape(permute(CCAs{u}(:,compi,:),[1 3 2]),[],length(compi));
                        else
                            CCA = reshape(CCAs{u}(:,:,compi),[],length(compi));
                        end
                    end
                case 4
                    CCA = reshape(CCAs{u}(:,o,o,compi),[],length(compi));
            end
            PCA = CCA*pinv(O(o).W(:,compi,u));
            
            if options.standardise_output
                tsigma = repmat(O(o).sigma{u},size(CCA,1),1);
                PCA = PCA.*tsigma;   
            end
            
            % add mean back
            if options.centre_output && centre_output
                tmu = repmat(O(o).mu{u},size(CCA,1),1);
                temp{u}(:,:,o) = PCA*pinv(O(o).COEFF{u}) + tmu;
            else
                temp{u}(:,:,o) = PCA*pinv(O(o).COEFF{u});
            end

        end
        
        dn = ndims(temp{u});
        if dn==2
            if strcmp(pcatype,'spatial')
                dat{u} = permute(reshape(temp{u},[],currdim(2),currdim(1)),[1 3 2]);
            elseif strcmp(pcatype,'temporal')
                dat{u} = reshape(temp{u},[],currdim(1),currdim(2));
            elseif strcmp(pcatype,'both')
                dat{u} = reshape(temp{u},[],currdim(1),currdim(2));
            end
        else
            dat{u}=temp{u};
        end
        
    end
end

grpdat = cat(1,dat{:});
dim = size(grpdat);
D.prep.dim = dim([2 3 1]);
grpdat = reshape(grpdat,size(grpdat,1),[]);
Y=struct;
for y=1:size(grpdat,2)
    Y(y).data_mean=mean(grpdat(:,y));
    Y(y).data_std=std(grpdat(:,y));
    % zscore
    Y(y).dtab=table;
    Y(y).dtab.data=(grpdat(:,y) - Y(y).data_mean) ./ Y(y).data_std;
end
D.prep.Y=Y;
sname = [strrep(num2str(cca_comp_select),'  ','')];
save(fullfile(pth,strrep(fname,'.mat',['_CCAcomp' sname '.mat'])),'D','S','-v7.3')

% check recontructions against original data
if 0
    E=load(fullfile(pth,'eegstats_dtab20200713T191743_orig.mat'));
    figure; scatter([D.prep.Y(:).data_mean],[E.D.prep.Y(:).data_mean])
    figure; scatter(D.prep.Y(1).dtab.data,E.D.prep.Y(1).dtab.data)
    figure; plot(D.prep.Y(1).dtab.data)
    figure; plot(E.D.prep.Y(1).dtab.data)
    rescale=0;
    do_boxplots = 0; % much slower!
    plotdiff(E.D,D,'ODD',{'group','CRPS'},rescale,do_boxplots,fname);
    plotdiff(E.D,D,'ODD',{'group','HC'},rescale,do_boxplots,fname);
    plotdiff(E.D,D,'ODD',{},rescale,do_boxplots,fname);
    plotdiff(E.D,D,'group',{},rescale,do_boxplots,fname);
    
end

function plotdiff(D1,D2,fac,selectfac,rescale,do_boxplots,fname)
% plots data WITHOUT RE-SCALING TO ORIGINAL MEAN/STD

S.img.file.coord = 'G:\Q_backup\Projects\CORE\eeg\egilocdata.mat';
S.img.file.coord_struct = 'c2d'; % name of structure within file that has coordinates
S.img.chansel = [2,3,4,5,6,7,9,10,11,12,13,15,16,18,19,20,22,23,24,26,27,28,29,30,31,33,34,35,36,37,39,40,41,42,45,46,47,50,51,52,53,54,55,57,58,59,60,61,62,65,66,67,70,71,72,75,76,77,78,79,80,83,84,85,86,87,90,91,92,93,96,97,98,100,101,102,103,104,105,106,108,109,110,111,112,115,116,117,118,122,123,124];
S.img.imgsize = 32;
S.img.interp_method = '';

% original
if ~isempty(selectfac)
    figure('name',[fac ', ' selectfac{1} num2str(selectfac{2})])
    [ulev,~,uind] = unique(D1.prep.dtab.(selectfac{1}));
    levind = find(ismember(ulev,selectfac{2}));
    dtab = D1.prep.dtab(uind==levind,:);
    for s=1:length(D1.prep.Y)
        D1.prep.Y(s).dtab=D1.prep.Y(s).dtab(uind==levind,:);
    end
else
    figure('name',fac)
    dtab = D1.prep.dtab;
end
[ulev,~,uind] = unique(dtab.(fac));
for s=1:length(D1.prep.Y)
    if rescale
        input(s,:) = bsxfun(@plus,bsxfun(@times,D1.prep.Y(s).dtab.data,D1.prep.Y(s).data_std),D1.prep.Y(s).data_mean);
    else
        input(s,:) = D1.prep.Y(s).dtab.data;
    end
    lev1(s,1)=mean(input(s, uind==1));
    lev2(s,1)=mean(input(s, uind==2));
end
if do_boxplots
    sname = strrep(fname,'.mat','_D1topo.mat');
    if ~exist(sname,'file')
        topo_input_2D = reshape(topotime_3D(reshape(input,D1.prep.dim(1),D1.prep.dim(2),[]),S),[],D1.prep.dim(2),size(input,2));
        save(sname,'topo_input_2D','-V7.3');
    else
        load(sname,'topo_input_2D');
    end
end
topo_lev1 = topotime_3D(reshape(lev1,D1.prep.dim(1),D1.prep.dim(2)),S);
topo_lev2 = topotime_3D(reshape(lev2,D1.prep.dim(1),D1.prep.dim(2)),S);
topo_lev1_2D=reshape(topo_lev1,[],size(topo_lev1,3));
topo_lev2_2D=reshape(topo_lev2,[],size(topo_lev2,3));
[~,maxdiff] = max(nanstd(topo_lev1_2D)-nanstd(topo_lev2_2D));
[~,maxvar] = max(topo_lev1_2D(:,maxdiff)-topo_lev2_2D(:,maxdiff));
subplot(3,5,1); pcolor(topo_lev1(:,:,maxdiff)), shading interp; axis off; title('lev1 orig')
subplot(3,5,2); pcolor(topo_lev2(:,:,maxdiff)), shading interp; axis off; title('lev2 orig')
subplot(3,5,3); plot(topo_lev1_2D(maxvar,:),'r'); hold on; plot(topo_lev2_2D(maxvar,:),'b'); 
if do_boxplots
    subplot(3,5,4); boxplot(squeeze(topo_input_2D(maxvar,maxdiff,:)),uind); title('trials') 
    % subject means
    [slev,~,sind] = unique(D1.prep.dtab.ID);
    subdat={};
    for u = 1:length(ulev)
        for sub = 1:length(slev)
            subdat{u}(sub,1) =mean(squeeze(topo_input_2D(maxvar,maxdiff, uind==u & sind==sub)));
        end
        subdat{u} = subdat{u}(~isnan(subdat{u}));
    end
    suind = [ones(length(subdat{1}),1);2*ones(length(subdat{2}),1)];
    subdat = vertcat(subdat{:});
    subplot(3,5,5); boxplot(subdat,suind); title('subjects') 
end
% CCA
if ~isempty(selectfac)
    [ulev,~,uind] = unique(D2.prep.dtab.(selectfac{1}));
    levind = find(ismember(ulev,selectfac{2}));
    dtab = D2.prep.dtab(uind==levind,:);
    for s=1:length(D2.prep.Y)
        D2.prep.Y(s).dtab=D2.prep.Y(s).dtab(uind==levind,:);
    end
else
    dtab = D2.prep.dtab;
end
[ulev,~,uind] = unique(dtab.(fac));
input=[];
for s=1:length(D2.prep.Y)
    if rescale
        input(s,:) = bsxfun(@plus,bsxfun(@times,D2.prep.Y(s).dtab.data,D2.prep.Y(s).data_std),D2.prep.Y(s).data_mean);
    else
        input(s,:) = D2.prep.Y(s).dtab.data;
    end
    lev1(s,1)=mean(input(s,uind==1));
    lev2(s,1)=mean(input(s,uind==2));
end

if do_boxplots
    sname = strrep(fname,'.mat','_D2topo.mat');
    if ~exist(sname,'file')
        topo_input_2D = reshape(topotime_3D(reshape(input,D2.prep.dim(1),D2.prep.dim(2),[]),S),[],D2.prep.dim(2),size(input,2));
        save(sname,'topo_input_2D','-V7.3');
    else
        load(sname,'topo_input_2D');
    end
end
topo_lev1 = topotime_3D(reshape(lev1,D2.prep.dim(1),D2.prep.dim(2)),S);
topo_lev2 = topotime_3D(reshape(lev2,D2.prep.dim(1),D2.prep.dim(2)),S);
topo_lev1_2D=reshape(topo_lev1,[],size(topo_lev1,3));
topo_lev2_2D=reshape(topo_lev2,[],size(topo_lev2,3));
subplot(3,5,6); pcolor(topo_lev1(:,:,maxdiff)), shading interp; axis off; title('lev1 CCA at orig')
subplot(3,5,7); pcolor(topo_lev2(:,:,maxdiff)), shading interp; axis off; title('lev2 CCA at orig')
subplot(3,5,8); plot(topo_lev1_2D(maxvar,:),'r'); hold on; plot(topo_lev2_2D(maxvar,:),'b'); 
if do_boxplots
    subplot(3,5,9); boxplot(squeeze(topo_input_2D(maxvar,maxdiff,:)),uind); title('trials') 
    % subject means
    [slev,~,sind] = unique(D2.prep.dtab.ID);
    subdat={};
    for u = 1:length(ulev)
        for sub = 1:length(slev)
            subdat{u}(sub,1) =mean(squeeze(topo_input_2D(maxvar,maxdiff, uind==u & sind==sub)));
        end
        subdat{u} = subdat{u}(~isnan(subdat{u}));
    end
    suind = [ones(length(subdat{1}),1);2*ones(length(subdat{2}),1)];
    subdat = vertcat(subdat{:});
    subplot(3,5,10); boxplot(subdat,suind); title('subjects') 
end
[~,maxdiff] = max(nanstd(topo_lev1_2D)-nanstd(topo_lev2_2D));
[~,maxvar] = max(topo_lev1_2D(:,maxdiff)-topo_lev2_2D(:,maxdiff));
subplot(3,5,11); pcolor(topo_lev1(:,:,maxdiff)), shading interp; axis off; title('lev1 CCA at CCA')
subplot(3,5,12); pcolor(topo_lev2(:,:,maxdiff)), shading interp; axis off; title('lev2 CCA at CCA')
subplot(3,5,13); plot(topo_lev1_2D(maxvar,:),'r'); hold on; plot(topo_lev2_2D(maxvar,:),'b'); 
if do_boxplots
    subplot(3,5,14); boxplot(squeeze(topo_input_2D(maxvar,maxdiff,:)),uind); title('trials') 
    % subject means
    [slev,~,sind] = unique(D2.prep.dtab.ID);
    subdat={};
    for u = 1:length(ulev)
        for sub = 1:length(slev)
            subdat{u}(sub,1) =mean(squeeze(topo_input_2D(maxvar,maxdiff, uind==u & sind==sub)));
        end
        subdat{u} = subdat{u}(~isnan(subdat{u}));
    end
    suind = [ones(length(subdat{1}),1);2*ones(length(subdat{2}),1)];
    subdat = vertcat(subdat{:});
    subplot(3,5,15); boxplot(subdat,suind); title('subjects') 
end