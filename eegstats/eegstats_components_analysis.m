function [D,grpdata]=eegstats_components_analysis(D,grpdata,S,varargin)
%- select temporal, spatial, both
%- select PCA/CCA method

% 1st level functions:
% - PCA
% - sPCA
% - factor
% - eigs
% - PLS

% 2nd level functions:
% - CCA

% inputs and outputs:
    % 'D': data structure containing experimental design table dtab, data
    % dimensions, etc.
    % 'grpdata': vertically concatenated EEG data array (over subjects)
    % Inputs are EEG chans/time/trials, outputs are components over trials.
    

% optional Y for PLS
if ~isempty(varagarin)
    Y=varargin{1};
else
    Y=[];
end

%%
[U,~,iU]=unique(D.prep.dtab.ID,'stable');
origgrpdata=grpdata;
newdim = D.prep.dim;

for pt = 1:length(S.type)
    NUM_FAC = S.type_numfac{pt};
    newgrpdata{pt}={};
    currdim = newdim;

    subdata={};
    dat={};
    switch S.data_in
        case 'conds'
            % Average over trials within each condition, per subject.
            % Currently assumes there are no missing conditions - needs
            % updating to allow missing conditions.
            [C,~,iC]=unique(double(D.prep.dtab.eventTypes));
            nreps = C(end);
            for u=1:length(U)
                for c=1:length(C)
                    dat{u}(c,:,:) = reshape(mean(grpdata{d}(iU==u & iC==c,:),1),1,currdim(1),currdim(2));
                end
                repidx{u} = C;
            end
        case 'trials'
            trialidx = unique([D.prep.tnums{:}]); % repetitions
            nreps = trialidx(end);
            for u=1:length(U)
                dat{u} = reshape(grpdata{d}(iU==u,:),[],currdim(1),currdim(2));
                repidx{u} = D.prep.tnums{u};
            end
    end

%         % baseline correct
%         if any(S.prep.calc.eeg.cca.base)
%             for u=1:length(U)
%                 if exist('select_samples','var')
%                     select = dsearchn(select_samples',[S.prep.calc.eeg.cca.base(1),S.prep.calc.eeg.cca.base(end)]');
%                 else
%                     select = dsearchn(total_samples',[S.prep.calc.eeg.cca.base(1),S.prep.calc.eeg.cca.base(end)]');
%                 end
%                 subdata{u} = bsxfun(@minus,subdata{u},mean(subdata{u}(:,select(1):select(2)),2));
%             end
%         end

    if strcmp(S.type{pt},'spatial')
        obs_dim = 2; % observation dimension
        var_dim = 1; % variable dimension
        % split observations?
        if S.type_obs_sep(pt)
            for obs = 1:currdim(2)
                for u=1:length(U)
                    subdata{obs}{u} = reshape(permute(dat{u}(:,:,obs),[2 1 3]),currdim(1),[]); 
                    subdata{obs}{u} = permute(subdata{obs}{u},[2 1]);
                end
            end
%                     newdim(obs_dim)=1;
        else
            % reshape
            for u=1:length(U)
                subdata{1}{u} = reshape(permute(dat{u},[2 1 3]),currdim(1),[]); % chan x trials*time
                subdata{1}{u} = permute(subdata{1}{u},[2 1]); % trials*time x chan 
            end
        end
    elseif strcmp(S.type{pt},'temporal')
        obs_dim = 1; % observation dimension
        var_dim = 2; % variable dimension
        % split observations?
        if S.type_obs_sep(pt)
            for obs = 1:currdim(1)
                for u=1:length(U)
                    subdata{obs}{u} = reshape(dat{u}(:,obs,:),[],currdim(2));
                end
            end
%                     newdim(obs_dim)=1;
        else
            for u=1:length(U)
                subdata{1}{u} = reshape(dat{u},[],currdim(2)); % trials*chan x time
            end
        end
    elseif strcmp(S.type{pt},'both')
        obs_dim = 0; % observation dimension
        var_dim = 0; % variable dimension
        for u=1:length(U)
            subdata{1}{u} = reshape(dat{u},[],currdim(1)*currdim(2)); % trials x chan*time
        end
    end

    % repetitions (usually time or chan, multiplied by conds or trials)
    for u=1:length(U)
        rep{u} = repidx{u};
        if obs_dim && length(subdata)==1
            for di = 2:newdim(obs_dim)
                rep{u} = [rep{u}, repidx{u} + (di-1)*nreps];
            end
        end
    end

    %%% PCA/CCA functions %%%
    O=struct;
    for o = 1:length(subdata)
        if ~obs_dim || length(subdata)>1
            obsdim=1;
        else
            obsdim=newdim(obs_dim);
        end
        
        % PCA
        [O(o).COEFF,O(o).W,O(o).scores,O(o).mu,O(o).sigma,NUM_FAC,YL] = perform_PCA(subdata{o},NUM_FAC,S,rep,nreps*obsdim,Y);
        O(o).scores_avg = squeeze(mean(reshape(O(o).scores,sz(1),nreps,obsdim,sz(3)),2));
        
        if S.CCA
            % CCA % https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.23689
            reg=1; r= 0.2;
            NUM_FAC(2) = min(NUM_FAC(2),NUM_FAC(1));
            [O(o).W,O(o).mdata] = mccas(O(o).scores,NUM_FAC(2),reg,r,O(o).W,S);
            sz=size(O(o).mdata);
            O(o).mdata_avg = squeeze(mean(reshape(O(o).mdata,sz(1),nreps,obsdim,sz(3)),2));
        end
    end

%             
    % get PCA and CCA scores for single trials
    for u=1:length(U)
        if strcmp(S.type{pt},'spatial')
            if length(O)==1
                temp = permute(reshape(grpdata{d}(iU==u,:),[],currdim(1),currdim(2)),[1 3 2]);% trials x time x chan
                temp = reshape(temp,[],size(temp,3)); 
            else
                temp = reshape(grpdata{d}(iU==u,:),[],currdim(1),currdim(2));% trials x chan x comp 
            end
        elseif strcmp(S.type{pt},'temporal')
            if length(O)==1
                temp = reshape(grpdata{d}(iU==u,:),[],currdim(1),currdim(2));% trials x chan x time
                temp = reshape(temp,[],size(temp,3)); 
            else
                temp = permute(reshape(grpdata{d}(iU==u,:),[],currdim(1),currdim(2)),[1 3 2]);% trials x time x comp
            end
        elseif strcmp(S.type{pt},'both')
            if length(O)==1
                temp = reshape(grpdata{d}(iU==u,:),[],currdim(1),currdim(2));% trials x chan x time
                temp = reshape(temp,[],size(temp,2)*size(temp,3)); 
            else
                temp = permute(reshape(grpdata{d}(iU==u,:),[],currdim(1),currdim(2)),[1 3 2]);% trials x time x chan
            end
        end

        if var_dim
            newdim(var_dim) = NUM_FAC(2);
        end
        
        
        oind{u}=[];
        PCAs{u} = [];
        for o = 1:length(O)

            if S.centre_output
                tmu = repmat(O(o).mu{u},size(temp,1),1);
                O(o).PCAs{u} = (temp(:,:,o)-tmu)*O(o).COEFF{u};
            else
                O(o).PCAs{u} = temp(:,:,o)*O(o).COEFF{u}; 
            end

            if S.standardise_output
                tsigma = repmat(O(o).sigma{u},size(temp(:,:,o),1),1);
                O(o).PCAs{u} = O(o).PCAs{u}./tsigma;   
            end

            oind{u} = [oind{u}, o*ones(1,size(O(o).PCAs{u},2))];
            PCAs{u} = [PCAs{u}, O(o).PCAs{u}];

        end
        
        if S.CCA
            CCAs{u} = [];
            for o = 1:length(O)
                CCAs{u} = [CCAs{u}, O(o).PCAs{u}*O(o).W(:,:,u)];
            end
            data_out=CCAs{u};
        else
            data_out=PCAs{u};
        end

        if length(O)>1
            temp = permute(reshape(data_out,length(repidx{u}),NUM_FAC(2),[]),[1 3 2]); % trials x time/chan x cc
        else
            temp = reshape(data_out,length(repidx{u}),[],NUM_FAC(2)); % trials x time/chan x cc
        end

        tempsz = [size(temp,2),size(temp,3)];

        if strcmp(S.type{pt},'spatial')
            temp = reshape(permute(temp,[1 3 2]),[],tempsz(1)*tempsz(2));
        elseif strcmp(S.type{pt},'temporal')
            temp = reshape(temp,[],tempsz(1)*tempsz(2));
        elseif strcmp(S.type{pt},'both')
            temp = squeeze(temp);
        end

        newgrpdata{pt}{d}(iU==u,:) = temp;
    end
    grpdata=newgrpdata{pt};
    if ~obs_dim
        newdim([1 2]) = tempsz;
    end

    % store outputs
    D.prep.PCA(pt).type = S.type{pt};
    D.prep.PCA(pt).O = O;
    D.prep.PCA(pt).oind = oind;
    D.prep.PCA(pt).currdim = currdim;
    D.prep.PCA(pt).options = S.prep.calc.eeg.cca;
    D.prep.PCA(pt).PCA = PCAs;
    if S.CCA
        D.prep.PCA(pt).CCA = CCAs;
    end

    % PLOT
    for o = 1:length(O)
        if S.CCA
            plot_mdata{o}=O(o).mdata;
            plot_mdata_avg{o}=O(o).mdata_avg;
        else
            plot_mdata{o}=O(o).scores;
            plot_mdata_avg{o}=O(o).scores_avg;
        end
    end
    
    cidx=floor(linspace(1,NUM_FAC(2),min(NUM_FAC(2),16)));
    for o = 1:length(O)
        figure('name',[S.type{pt} ', obs ' num2str(o)]); clf;
        cci=0; chandatacorr=[];
        for cc = 1:NUM_FAC(2)
            if strcmp(S.type{pt},'spatial')

                % topo data
                O(o).chandata=[];
                for u=1:length(U)
                    O(o).chandata(cc,:,u) = O(o).COEFF{u}*O(o).W(:,cc,u);
                end

                % PLOTS
                if ismember(cc,cidx)
                    % waveform plot
                    cci=cci+1;
                    mCCA = mean(plot_mdata_avg{o}(cc,:,:),3);
                    if length(mCCA) == length(S.prep.select.samples)
                        subplot(4,4*2,cci)
                        hold on
                        for u=1:length(U)
                            % plot(CCA{u}(:,cc));
                            plot(plot_mdata_avg{o}(cc,:,u));
                        end
                        % mean
                        % ccall = cellfun(@(X) X(:,cc), CCA,'UniformOutput',0);
                        % mCCA = mean(horzcat(ccall{:}),2);
                        plot(mCCA,'k','LineWidth',2);
                        hold off
                    end

                    % topoplot
                    cci=cci+1;
                    subplot(4,4*2,cci)
                    mchandata = squeeze(mean(O(o).chandata(cc,:,:),3))';
                    topo = topotime_3D(mchandata,S);
                    pcolor(topo), shading interp; axis off
                    title(['comp ' num2str(cc)])

%                             % topoplot from correlation
%                             cci=cci+1;
%                             subplot(4,4*2,cci)
%                             temp = [];
%                             for u=1:length(U)
%                                 s_subdata = reshape(newgrpdata{1}{d}(iU==u,:),[],40,19);
%                                 temp(1,:,u) = corr(plot_mdata{o}(cc,repidx{u},u)',s_subdata(:,:,o)); % 3x60x20 x trials(57)xComp(1083=57x19)
%                             end
%                             chandatacorr(cc,:) = mean(temp,3); % mean over subjects
%                             topo = topotime_3D(chandatacorr(cc,:)',S);
%                             pcolor(topo), shading interp; axis off
%                             title(['comp corr ' num2str(cc)])
                end

            elseif strcmp(S.type{pt},'temporal')

                % waveform data
                for u=1:length(U)
                    O(o).timedata(cc,:,u) = O(o).COEFF{u}(:,1:NUM_FAC(1))*O(o).W(1:NUM_FAC(1),cc,u);  
                end

                % PLOTS
                if ismember(cc,cidx)
                    % waveform plot
                    cci=cci+1;
                    subplot(4,4*2,cci)
                    hold on
                    for u=1:length(U)
                        plot(O(o).timedata(cc,:,u));
                    end
                    mtimedata = squeeze(mean(O(o).timedata(cc,:,:),2));
                    plot(mtimedata,'k','LineWidth',2);
                    hold off

                    % topoplot
                    cci=cci+1;
                    subplot(4,4*2,cci)
                    mCCA = nanmean(plot_mdata_avg{o}(cc,:,:),3)';
                    if length(mCCA) == length(S.img.chansel)
                        topo = topotime_3D(mCCA,S);
                        pcolor(topo), shading interp; axis off
                    end
                    title(['comp ' num2str(cc)])
                end

            end

        end
        D.prep.PCA(pt).O(o).chandatacorr = chandatacorr;
    end


end

% final chan-time data
for o = 1:length(O)
    for cc = 1:NUM_FAC(2)
        temp = [];
        for u=1:length(U)
            temp(1,:,u) = corr(plot_mdata{o}(cc,repidx{u},u)',origgrpdata{d}(iU==u,:),'type','Spearman'); % 3x60x20 x trials(57)xComp(1083=57x19)
        end
        temp = nanmean(temp,3);
        st_subdata = reshape(temp,D.prep.dim(1),D.prep.dim(2));
        chantimecorr(cc,:,:,o) = st_subdata; % mean over subjects
    end
end
D.prep.PCA(pt).chantimecorr = chantimecorr;

% plot - spatiotemporal maxima
%close all
for cc = 1:NUM_FAC(2)
    dat=squeeze(chantimecorr(cc,:,:));
    topo = topotime_3D(dat,S);
    thresh = [nanmean(topo(:))-2*nanstd(topo(:)),nanmean(topo(:))+2*nanstd(topo(:))];

    dat=topo;
    dat(dat>thresh(1) & dat<thresh(2)) = nan;
    dat(isnan(dat))=0;
    % connected components
    ccon = bwconncomp(dat,26);
    % remove small clusters
    ClusExtent = cellfun(@length,ccon.PixelIdxList);
    ccon.PixelIdxList(ClusExtent<max(ClusExtent)/10)=[];
    nc = length(ccon.PixelIdxList);

    % plot
    figure('name',['component ' num2str(cc)]); 
    ax(1) = subplot(nc+1,2,2); hold on
    plot(squeeze(max(topo,[],[1 2])),'b')
    plot([1,size(topo,3)],[thresh(2),thresh(2)],'b--')
    plot(squeeze(min(topo,[],[1 2])),'r')
    plot([1,size(topo,3)],[thresh(1),thresh(1)],'r--')

    % plot maps
    pi=2;
    for cci = 1:nc
        % subscripts
        [~,i] = max(abs(topo(ccon.PixelIdxList{cci})));
        [x,y,z]=ind2sub(size(topo),ccon.PixelIdxList{cci}(i));

        % line
        plot(ax(1),[z,z],[0,topo(x,y,z)],'k')

        % topo
        pi=pi+1;
        subplot(nc+1,2,pi); hold on
        pcolor(topo(:,:,z)), shading interp; axis off
        scatter(y,x,'k','filled')
        hold off

        % waveform
        pi=pi+1;
        subplot(nc+1,2,pi); hold on
        plot(squeeze(topo(x,y,:)),'k');
        if topo(x,y,z)>0
            col='b';
        else
            col='r';
        end
        plot([z,z],[0,topo(x,y,z)],col)
        hold off
        ylabel(['z=' num2str(z)])
    end
end


% plot - spatial and temporal variance maxima
% based on the idea that CCAs are multivariate and maximise
% variance, except this is variance over topotime rather than over
% trials.
%close all
for cc = 1:NUM_FAC(2)
    dat=squeeze(chantimecorr(cc,:,:));
    topo = topotime_3D(dat,S);
    sz=size(topo);

    stmask=[];
    sd=2;
    while isempty(find(stmask))
        sd = sd-0.1;
        % temporal regions of high variance over space
        temp = reshape(topo,[],sz(3));
        t_gfp = nanstd(temp,[],1);
        tthresh = mean(t_gfp)+sd*std(t_gfp);

        % spatial regions of high variance over time
        s_gfp = nanstd(topo,[],3);
        s_gfp(isnan(topo(:,:,1))) = nan;
        sthresh = nanmean(s_gfp(:))+sd*nanstd(s_gfp(:));

        % spatiotemporal mask
        tmask = permute(repmat(t_gfp>=tthresh,sz(1),1,sz(2)),[1 3 2]);
        smask = repmat(s_gfp>=sthresh,1,1,sz(3));
        stmask = smask.*tmask;
    end

    dat=stmask;
    dat(dat>thresh(1) & dat<thresh(2)) = nan;
    dat(isnan(dat))=0;
    % connected components
    ccon = bwconncomp(dat,26);
    % remove small clusters
    ClusExtent = cellfun(@length,ccon.PixelIdxList);
    ccon.PixelIdxList(ClusExtent<max(ClusExtent)/10)=[];
    nc = length(ccon.PixelIdxList);

    % plot
    figure('name',['component ' num2str(cc)]); 
    pi=0;
    % s_gfp
    pi=pi+1;
    subplot(nc+1,2,pi); hold on
    pcolor(s_gfp), shading interp; axis off
    hold off
    title('temporal std')
    % t_gfp
    pi=pi+1;
    subplot(nc+1,2,pi); hold on
    plot(t_gfp,'k');
    plot([1,length(t_gfp)],[tthresh,tthresh],'b--')
    hold off
    title('spatial std')
    for cci = 1:nc
        % subscripts
        [x,y,z]=ind2sub(sz,ccon.PixelIdxList{cci});
        max_zi = find(t_gfp==max(t_gfp(z)));
        max_xyi = find(s_gfp==max(s_gfp(x,y),[],'all'));

        % topo
        pi=pi+1;
        subplot(nc+1,2,pi); hold on
        pcolor(topo(:,:,max_zi)), shading interp; axis off
        scatter(y,x,'k')
        hold off
        title('component topography')

        % waveform
        pi=pi+1;
        subplot(nc+1,2,pi); hold on
        [xw,yw] = ind2sub(sz([1 2]),max_xyi);
        plot(squeeze(topo(xw,yw,:)),'k');
        plot([max_zi,max_zi],[0,topo(xw,yw,max_zi)],'b')
        hold off
        ylabel(['z=' num2str(max_zi)])
        title('component waveform')
    end
end

D.prep.dim = newdim;

        

function [COEFF,W,PCdata,mu,sigma,NUM_FAC,YL] = perform_PCA(dataAll,NUM_FAC,S,rep,nobs,Y)
% feed in index of 2nd dimension (trials or time)
               
TEMP_FAC=repmat(size(dataAll{1},2),1,2);
if isempty(NUM_FAC)
    NUM_FAC([1,2]) = TEMP_FAC;
else
    NUM_FAC([1,2]) = min([NUM_FAC;TEMP_FAC]);
end

maxmat = floor((S.maxGig*1e9)^(1/2));
if size(dataAll,2)>maxmat
    error(['downsample the data by ' num2str(size(dataAll,2)/maxmat) ' to enable CCA'])
end

% centre data
cdata={};
for i = 1:length(dataAll)
    tmpscore = squeeze(dataAll{i});
    mu{i} = mean(tmpscore);
    cdata{i} = tmpscore - repmat(mu{i},size(tmpscore,1),1);
end

YL=[];
% apply MCCA
switch S.PCAmethod

    case 'FA'
        % find nfac for each subject
        for i = 1:length(cdata)
            FactorResults = erp_pca(cdata{i},NUM_FAC(1));
            randResults = erp_pca(randn(size(cdata{i})),NUM_FAC(1));
            randResults.scree = randResults.scree*(sum(FactorResults.scree)/sum(randResults.scree)); %scale to same size as real data
            nfac_temp = find(FactorResults.scree<randResults.scree);
            nfac(i) = min(NUM_FAC(1),nfac_temp(1)-1);
        end
        % take median nfac
        NUM_FAC(1) = floor(median(nfac));
        for i = 1:length(cdata)
            FactorResults = erp_pca(cdata{i},NUM_FAC(1));
            COEFF{i} = FactorResults.FacCof;
            W(:,:,i) = COEFF{i}';
            score{i} = FactorResults.FacScr;
        end
    case 'PCA'
       for i = 1:length(cdata)
            [COEFF{i}, score{i},~,~,explained] = pca(cdata{i},'NumComponents',NUM_FAC(1),'Centered',false,'Algorithm','eig');
            [COEFFrand{i}, scorerand{i},~,~,explainedrand] = pca(randn(size(cdata{i})),'NumComponents',NUM_FAC(1));
            nfac_temp = find(explained<explainedrand);
            nfac(i) = min(NUM_FAC(1),nfac_temp(1)-1);
       end
        % take median nfac
        NUM_FAC(1) = floor(median(nfac));
        for i = 1:length(cdata)
            COEFF{i} = COEFF{i}(:,1:NUM_FAC(1));
            W(:,:,i) = COEFF{i}';
            score{i} = score{i}(:,1:NUM_FAC(1));
        end
    case 'eigs'
       for i = 1:length(cdata)
            [COEFF{i},values]=eigs(cdata{i},NUM_FAC(1));
            explained = diag(values)/sum(diag(values));
            [COEFFrand{i},valuesrand]=eigs(randn(size(cdata{i})),NUM_FAC(1));
            explainedrand = diag(valuesrand)/sum(diag(valuesrand));
            nfac_temp = find(explained<explainedrand);
            nfac(i) = min(NUM_FAC(1),nfac_temp(1)-1);
       end
        % take median nfac
        NUM_FAC(1) = floor(median(nfac));
        for i = 1:length(cdata)
            COEFF{i} = COEFF{i}(:,1:NUM_FAC(1));
            W(:,:,i) = COEFF{i}';
            score{i} = cdata{i}*COEFF{i}(:,1:NUM_FAC(1));
        end
    case 'sPCA'
        
        delta = inf;
        maxiter = 3000;
        convergenceCriterion = 1e-9;
        verbose = false;
        
        for i = 1:length(cdata)
            [COEFF{i},values] = spca(cdata{i}, [], NUM_FAC(1), delta, -NUM_FAC(1), maxiter, convergenceCriterion, verbose);
            W(:,:,i) = COEFF{i}';
            score{i} = cdata{i}*COEFF{i};
        end
    case 'PLS'
        [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(cdata{i},Y,NUM_FAC(1));
end
% obtain subject-specific PCs
PCdata = nan(NUM_FAC(1),nobs);
sigma = {};
for i = 1:length(cdata)
    if S.CCA % standardise scores if they will be used for CCA
       sigma{i} = sqrt(var(score{i})); 
       score{i} = score{i}./repmat(sqrt(var(score{i})),size(score{i},1),1); % standardise so that each PC from each subject contributes equally to CCA solution
    end
    PCdata(:,rep{i},i) = score{i}';
end


function [W, mdata, mWeights] = mccas(data,K,reg,r,weights,Sx)
%% regularized-multiset CCA
%   Date           Programmers               Description of change
%   ====        =================            =====================
%  09/10/2016     Qiong Zhang                 Original code
%% Citation
%  Zhang,Q., Borst,J., Kass, R.E., & Anderson, J.A. (2016) Between-Subject
%  Alignment of MEG Datasets at the Neural Representational Space. 

%% INPUT
% data - input training data to CCA (samples * sensors * subjects)
% K - number of CCAs to retain
% reg - add regularization (1) or not (0)
% r - degree of regularization added
% weights - PCA weight maps for each subject

%% OUTPUT
% W - CCA weights that transform PCAs to CCAs for each subject (PCAs * CCAs * subjects)
% mdata - transformed averaged data (CCAs * sapmles * subjects)
% mWeights - projection weights from sensor to CCAs (CCAs * sensors * subjects)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim = size(data,1);
num = size(data,2);
sub = size(data,3);
dim2 = size(weights,1);
num2 = size(weights,2);
sub2 = size(weights,3);
if(K>dim)
    K = dim;
end
temp=[];
for i=1:sub
    temp=[temp;data(:,:,i)];
end
R=cov(temp.','partialrows'); % cab: added partialrows to allow missing data (nan)
S = zeros(size(R));
for i = 1:dim*sub
    tmp = ceil(i/dim);
    S((tmp-1)*dim+1:tmp*dim,i) = R((tmp-1)*dim+1:tmp*dim,i); 
end

% add regularization 
if(reg==1)    
    if(K>dim2)
        K = dim2;
    end
    temp=[];
    for i=1:sub2
        temp=[temp;weights(:,:,i)];
    end
    R2=cov(temp.');
    S2 = zeros(size(R2));
    for i = 1:dim2*sub2
        tmp = ceil(i/dim2);
        S2((tmp-1)*dim2+1:tmp*dim2,i) = R2((tmp-1)*dim2+1:tmp*dim2,i); 
    end
    R = R + r*R2;
    S = S + r*S2;
end

% obtain CCA solutions 
[tempW,values]=eigs((R-S),S,K);

W = zeros(dim,K,sub);
for i=1:sub
    if Sx.normalise_cca_weights 
        W(:,:,i)=tempW((i-1)*dim+1:i*dim,:)./norm(tempW((i-1)*dim+1:i*dim,:));
    else
        W(:,:,i)=tempW((i-1)*dim+1:i*dim,:);
    end
end

% projected data
mdata = zeros(K,num,sub);
for i = 1:sub
    mdata(:,:,i) = W(:,:,i)'*data(:,:,i);
    for j = 1:K
        if Sx.normalise_cca_weights
            mdata(j,:,i) = mdata(j,:,i)/norm(mdata(j,:,i));
        else
            mdata(j,:,i) = mdata(j,:,i);
        end
    end
end

% projection weights
mWeights = 0;
if(reg==1)
    mWeights = zeros(K,num2,sub2);
    for i = 1:sub2
        mWeights(:,:,i) = W(:,:,i)'*weights(:,:,i);
    end
end