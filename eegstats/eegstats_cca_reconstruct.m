function eegstats_cca_reconstruct
dbstop if error
close all
pth='C:\Data\CORE\eeg\ana\prep\cleaned\part2\eegstats_dataprep';
fname = 'eegstats_dtab20200330T125514_temporalboth.mat';

load(fullfile(pth,fname))
cca_comp_select = [];
d=1;
[U,~,iU]=unique(D(d).prep.dtab.ID,'stable');
dat={};

% e.g.
% 12 temporal
% 12 x 5 spatial
% 17 both
for pt = length(D(d).prep.CCA):-1:1
    
    % unpack
    pcatype = D(d).prep.CCA(pt).pcatype;
    O = D(d).prep.CCA(pt).O;
    oind = D(d).prep.CCA(pt).oind;
    currdim = D(d).prep.CCA(pt).currdim;
    options = D(d).prep.CCA(pt).options;
    nCCA=size(O(1).W,2);
    
    if pt==length(D(d).prep.CCA)
        CCAs = D(d).prep.CCA(pt).CCA;
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
            
            if options.standardise
                tsigma = repmat(O(o).sigma{u},size(CCA,1),1);
                PCA = PCA.*tsigma;   
            end
            
            % add mean back
            tmu = repmat(O(o).mu{u},size(CCA,1),1);
            temp{u}(:,:,o) = PCA*pinv(O(o).COEFF{u}) + tmu;

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
save(fullfile(pth,strrep(fname,'.mat',['_CCAcomp' sname '.mat'])),'D','-v7.3')

% check recontructions against original data
if 1
    close all
    E=load('eegstats_dtab20200325T161941_noCCA.mat');
    figure; scatter([D.prep.Y(:).data_mean],[E.D.prep.Y(:).data_mean])
    figure; scatter(D.prep.Y(1).dtab.data,E.D.prep.Y(1).dtab.data)
    figure; plot(D.prep.Y(1).dtab.data)
    figure; plot(E.D.prep.Y(1).dtab.data)
    plotdiff(E.D,D,'ODD');
    plotdiff(E.D,D,'group');
    
end

function plotdiff(D1,D2,fac)

S.img.file.coord = 'G:\Q_backup\Projects\CORE\eeg\egilocdata.mat';
S.img.file.coord_struct = 'c2d'; % name of structure within file that has coordinates
S.img.chansel = [2,3,4,5,6,7,9,10,11,12,13,15,16,18,19,20,22,23,24,26,27,28,29,30,31,33,34,35,36,37,39,40,41,42,45,46,47,50,51,52,53,54,55,57,58,59,60,61,62,65,66,67,70,71,72,75,76,77,78,79,80,83,84,85,86,87,90,91,92,93,96,97,98,100,101,102,103,104,105,106,108,109,110,111,112,115,116,117,118,122,123,124];
S.img.imgsize = 32;
S.img.interp_method = '';


figure('name',fac)
% original
[ulev,~,uind] = unique(D1.prep.dtab.(fac));
for y=1:length(D1.prep.Y)
    lev1(y,1)=mean(D1.prep.Y(y).dtab.data(uind==1));
    lev2(y,1)=mean(D1.prep.Y(y).dtab.data(uind==2));
end
topo_lev1 = topotime_3D(reshape(lev1,D1.prep.dim(1),D1.prep.dim(2)),S);
topo_lev2 = topotime_3D(reshape(lev2,D1.prep.dim(1),D1.prep.dim(2)),S);
topo_lev1_2D=reshape(topo_lev1,[],size(topo_lev1,3));
topo_lev2_2D=reshape(topo_lev2,[],size(topo_lev2,3));
[~,maxdiff] = max(nanstd(topo_lev1_2D)-nanstd(topo_lev2_2D));
[~,maxvar] = max(topo_lev1_2D(:,maxdiff)-topo_lev2_2D(:,maxdiff));
subplot(3,3,1); pcolor(topo_lev1(:,:,maxdiff)), shading interp; axis off; title('lev1 orig')
subplot(3,3,2); pcolor(topo_lev2(:,:,maxdiff)), shading interp; axis off; title('lev2 orig')
subplot(3,3,3); plot(topo_lev1_2D(maxvar,:),'r'); hold on; plot(topo_lev2_2D(maxvar,:),'b'); 
% CCA
[ulev,~,uind] = unique(D2.prep.dtab.(fac));
for y=1:length(D2.prep.Y)
    lev1(y,1)=mean(D2.prep.Y(y).dtab.data(uind==1));
    lev2(y,1)=mean(D2.prep.Y(y).dtab.data(uind==2));
end
topo_lev1 = topotime_3D(reshape(lev1,D2.prep.dim(1),D2.prep.dim(2)),S);
topo_lev2 = topotime_3D(reshape(lev2,D2.prep.dim(1),D2.prep.dim(2)),S);
topo_lev1_2D=reshape(topo_lev1,[],size(topo_lev1,3));
topo_lev2_2D=reshape(topo_lev2,[],size(topo_lev2,3));
subplot(3,3,4); pcolor(topo_lev1(:,:,maxdiff)), shading interp; axis off; title('lev1 CCA at orig')
subplot(3,3,5); pcolor(topo_lev2(:,:,maxdiff)), shading interp; axis off; title('lev2 CCA at orig')
subplot(3,3,6); plot(topo_lev1_2D(maxvar,:),'r'); hold on; plot(topo_lev2_2D(maxvar,:),'b'); 
[~,maxdiff] = max(nanstd(topo_lev1_2D)-nanstd(topo_lev2_2D));
[~,maxvar] = max(topo_lev1_2D(:,maxdiff)-topo_lev2_2D(:,maxdiff));
subplot(3,3,7); pcolor(topo_lev1(:,:,maxdiff)), shading interp; axis off; title('lev1 CCA at CCA')
subplot(3,3,8); pcolor(topo_lev2(:,:,maxdiff)), shading interp; axis off; title('lev2 CCA at CCA')
subplot(3,3,9); plot(topo_lev1_2D(maxvar,:),'r'); hold on; plot(topo_lev2_2D(maxvar,:),'b'); 