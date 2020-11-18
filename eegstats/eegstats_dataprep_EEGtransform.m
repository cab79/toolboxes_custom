function D = eegstats_dataprep_EEGtransform(S,varargin)
% inputs: S is a settings structure (see example script)
% outputs: D is the prepared EEG data 

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

%% transformations: must take place on grouped data (over subjects) if it is grouped
for d = 1:length(D)
    disp(['data transformations: ' num2str(d) '/' num2str(length(D))])
    
    if S.prep.calc.eeg.pca.on
        
        %% use existing data or calculate instead?
        if isfield(S.prep.calc.eeg.pca,'use_existing_data') && ~isempty(S.prep.calc.eeg.pca.use_existing_data)
            Dy = load(S.prep.calc.eeg.pca.use_existing_data,'D');
            D.prep.PCA = Dy.D.prep.PCA;
            D.prep.Y = Dy.D.prep.Y;
            D.prep.dim = Dy.D.prep.dim;
        else
        
            %% inputs for PLS/CCA
            Sf.PCAmethod = S.prep.calc.eeg.pca.PCAmethod;
            Sf.type = S.prep.calc.eeg.pca.type;
            Sf.type_obs_sep = S.prep.calc.eeg.pca.type_sep; % analyse "observation" separately. E.g. if temporal PCA, apply separately to each spatial channel?
            Sf.type_numfac = S.prep.calc.eeg.pca.numfac; % max number of factors
            Sf.data_in = S.prep.calc.eeg.pca.data_in; % type of data to analyse: trials or averaged conditions
            Sf.centre_output = S.prep.calc.eeg.pca.centre;
            Sf.standardise_output = S.prep.calc.eeg.pca.standardise;
            Sf.pca_FA_type = S.prep.calc.eeg.pca.FA_type;
            Sf.maxGig = S.prep.calc.eeg.pca.maxGig;
            Sf.normalise_cca_weights = S.prep.calc.eeg.cca.normalise_weights;
            Sf.img=S.img;
            Sf.Y_use4output = S.prep.calc.eeg.pls.Y_use4output;
            Sf.select_ncomp = S.prep.calc.eeg.pca.select_ncomp;
            Sf.cca_reg_weights = S.prep.calc.eeg.cca.reg_weights;
            Sf.cca_test_nPCAcomp = S.prep.calc.eeg.cca.test_nPCAcomp;
            Sf.cca_method = S.prep.calc.eeg.cca.method;
            Sf.cca_reduce_by_scree = S.prep.calc.eeg.cca.reduce_by_scree;
            Sf.cca_select_nPCA_per_group = S.prep.calc.eeg.cca.select_nPCA_per_group;
            Sf.cca_FA_type = S.prep.calc.eeg.cca.FA_type;

            % PLS
            Y={};
            if any(strcmp(Sf.PCAmethod,'PLS'))
                if S.prep.calc.eeg.pca.on
                    Y = D(d).prep.pred_PCA.Yscore;
                    col_names = D(d).prep.pred_PCA.col_names;
                else
                    error('run PCA on predictors prior to PLS')
                end
            end

            % RUN function
            if any(strcmp(Sf.PCAmethod,'PLS')) && S.prep.calc.eeg.pls.Y_select % select a model and output that single model
                [D(d)]=eegstats_components_analysis(D(d),Sf,{Y,col_names});
            elseif any(strcmp(Sf.PCAmethod,'PLS')) % output all models
                D_orig=D;
                for ym = 1:length(Y)
                    [Dtemp]=eegstats_components_analysis(D_orig(d),Sf,{Y(ym),col_names(ym)});
                    D(d).prep(ym) = D(d).prep(1);
                    D(d).prep(ym).grpdata=Dtemp.prep.grpdata;
                    D(d).prep(ym).PCA=Dtemp.prep.PCA;
                    D(d).prep(ym).dim=Dtemp.prep.dim;
                    close all
                    mse_temp=mean(cat(3,D(d).prep(ym).PCA.O.PLS.MSE{:}),3); 
                    msecv_temp=mean(cat(3,D(d).prep(ym).PCA.O.PLS.MSE_CV{:}),3); 
                    pvar_temp=mean(cat(3,D(d).prep(ym).PCA.O.PLS.explained{:}),3); 
                    mse_plot(ym,:) = mse_temp(:,end);
                    msecv_plot(ym,:) = msecv_temp(:,end);
                    pvar_plot(ym,:) = sum(pvar_temp,2);
                end
                figure
                subplot(2,3,1); bar(mse_plot(:,1)); title('X MSE per model')
                subplot(2,3,2); bar(msecv_plot(:,1)); title('X CV MSE per model')
                subplot(2,3,3); bar(pvar_plot(:,1)); title('X % var explained per model')
                subplot(2,3,4); bar(mse_plot(:,2)); title('Y MSE per model')
                subplot(2,3,5); bar(msecv_plot(:,2)); title('Y CV MSE per model')
                subplot(2,3,6); bar(pvar_plot(:,2)); title('Y % var explained per model')
                clear Dtemp D_orig
            else
                [D(d)]=eegstats_components_analysis(D(d),Sf);
            end
        end
    else
        for ip = 1:length(D(d).prep)
            D(d).prep(ip).grpdata = D(d).prep(ip).grpdata{1};
        end
    end
    
end

% save data to disk
if S.prep.output.save
    disp('saving data to disk...')
    save(fullfile(S.prep.path.outputs,[S.prep.sname '.mat']),'D','S','-v7.3')
    disp('...done')
end

