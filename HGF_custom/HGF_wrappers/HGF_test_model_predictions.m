function varargout=HGF_test_model_predictions(D_act,D_sim,S, type) 
% outputs summary stats of model predictions of behaviour over trials

% Bayesian regularised regression (BRR) settings
S.brr.folds = 0;            % number of folds in the traindata. Set to 0 to not conduct predictive cross-validation.
S.brr.model = 'gaussian';   % error distribution - string, one of {'gaussian','laplace','t','binomial'}
S.brr.prior = 'g';        %- string, one of {'g','ridge','lasso','horseshoe','horseshoe+'}
S.brr.nsamples = 1000;   %- number of posterior MCMC samples (Default: 1000)  
S.brr.burnin = 1000;     %- number of burnin MCMC samples (Default: 1000)
S.brr.thin = 5;       %- level of thinning (Default: 5)
S.brr.catvars = [];    %- vector of variables (as column numbers of X) to treat
%                       as categorical variables, with appropriate expansion. 
%                       See examples\br_example5 (Default: none)
S.brr.nogrouping = false; %- stop automatic grouping of categorical predictors
%                       that is enabled with the 'catvars' options. (Default: false)
S.brr.usegroups = 0;     % ****Specified by S.traj cells**** - create groups of predictors. Grouping of variables
%                       works only with HS, HS+ and lasso prior
%                       distributions. The same variable can appear in
%                       multiple groups. See examples\br_example[9,10,11,12]  (Default: { [] } )  
S.brr.waic = true;       %- whether to calculate the WAIC -- disabling can lead
%                       to large speed-ups, especially for Gaussian models with large n
%                       (default: true)
S.zscore = 1;

% do this in case of missing conditions in some ppts
allcond=1;
for d = 1:length(D_act)
    allcond = unique([allcond, unique(D_act(d).Processed.condnum{1})]);
end
actual_condmeans=nan(length(D_act),length(allcond));
fitted_condmeans=nan(length(D_act),length(allcond));

switch type
    case 'Ch' % choices
        actual_macorrect=[];
        fitted_macorrect=[];
        for d = 1:length(D_act)
            try
                actual_choices = D_act(d).HGF(1).yALL(:,1);
            catch
                actual_choices = D_act(d).HGF(1).y(:,1);
            end

            % condition fraction correct
            condnum = D_act(d).Processed.condnum{1};
            ii = ~isnan(actual_choices); 
            unicond = unique(condnum(ii));
            actual_condmeans(d,unicond) = mean_by_condition(actual_choices(ii)==D_act(d).HGF(1).u(ii,1),condnum(ii)); 

            % split into stim intensity per block
            blocknum = D_act(d).Sequence.blocks;
            stimnum = D_act(d).Sequence.signal(S.signal.target,:);
            [~,~,blockstimnum] = unique([blocknum;stimnum]','rows');
            uniblockstim = unique(blockstimnum);
            actual_blockstimmeans(d,uniblockstim) = mean_by_condition(actual_choices(ii)==D_act(d).HGF(1).u(ii,1),blockstimnum(ii)); 

            % MA correct
            % moving average percent correct over time
            ma=S.movingavg;
            if ~isempty(ma) 
                if length(ma) ~=2
                    error('specify min and max values for moving average of % correct')
                end
                for i = 1:length(D_act(d).Processed.correct)
                    ind = i-(ma(2)+1):i;
                    if sum(~isnan(D_act(d).Processed.correct(ind(ind>0)))) < ma(1)
                        actual_macorrect(d,i) = NaN;
                    else
                        actual_macorrect(d,i) = 100*nansum(D_act(d).Processed.correct(ind(ind>0)))/length(ind(ind>0));
                    end
                end
            end

            fitted_choices_rep = nan(1,length(actual_choices));
            fitted_condmeans_rep = [];
            fitted_blockstimmeans_rep = [];
            fitted_correct = [];
            rand_correct = [];
            temp_macorrect=[];
            for rep = 1:S.numsimrep
                fitted_choices = D_sim(d).HGF(rep).sim.y(:,1);
                fitted_choices_rep(rep,ii) = fitted_choices(ii);
                fitted_condmeans_rep(rep,unicond) = mean_by_condition(fitted_choices(ii)==D_act(d).HGF(1).u(ii,1),condnum(ii)); 
                fitted_correct(rep) = sum(fitted_choices(ii) == actual_choices(ii))/length(actual_choices(ii));
                fitted_blockstimmeans_rep(rep,uniblockstim) = mean_by_condition(fitted_choices(ii)==D_act(d).HGF(1).u(ii,1),blockstimnum(ii)); 

                % to create benchmark for predictions, randomise simulated
                % responses over trials 

%                 %separately for each stimulus type (0 and 1)
%                 rand_choices = nan(length(fitted_choices),1);
%                 u0i = find(D_sim(d).HGF(rep).sim.u(:,1)==0);
%                 u1i = find(D_sim(d).HGF(rep).sim.u(:,1)==1);
%                 rand_choices(u0i) = D_sim(d).HGF(rep).sim.y(u0i(randperm(length(u0i))),1);
%                 rand_choices(u1i) = D_sim(d).HGF(rep).sim.y(u1i(randperm(length(u1i))),1);

                %over both stimulus types (0 and 1)
                rand_choices = fitted_choices(randperm(length(fitted_choices)),1);
                rand_correct(rep) = sum(rand_choices(ii) == actual_choices(ii))/length(actual_choices(ii));

                % MA correct
                correct = nan(length(fitted_choices),1);
                correct(ii) = fitted_choices(ii)==D_act(d).HGF(1).u(ii,1);
                if ~isempty(ma) 
                    if length(ma) ~=2
                        error('specify min and max values for moving average of % correct')
                    end
                    for i = 1:length(D_act(d).Processed.correct)
                        ind = i-(ma(2)+1):i;
                        if sum(~isnan(correct(ind(ind>0)))) < ma(1)
                            temp_macorrect(rep,i) = NaN;
                        else
                            temp_macorrect(rep,i) = 100*nansum(correct(ind(ind>0)))/length(ind(ind>0));
                        end
                    end
                end
            end
            actual_choices_all(d,:) = actual_choices;
            fitted_choices_mean(d,:) = mean(fitted_choices_rep,1);
            means_fitted(d,1)=mean(fitted_correct);
            stds_fitted(d,1)=std(fitted_correct);
            fitted_macorrect(d,:)=mean(temp_macorrect,1);
            % random-corrected
            means_randcorr_fitted(d,1)=mean(fitted_correct-rand_correct);
            stds_randcorr_fitted(d,1)=std(fitted_correct-rand_correct);
            fitted_condmeans(d,:) = mean(fitted_condmeans_rep,1);
            fitted_blockstimmeans(d,:) = mean(fitted_blockstimmeans_rep,1);
        end
        varargout = {means_fitted,stds_fitted,means_randcorr_fitted,stds_randcorr_fitted,actual_condmeans,fitted_condmeans,actual_blockstimmeans,fitted_blockstimmeans,actual_macorrect,fitted_macorrect,actual_choices_all,fitted_choices_mean};

    case 'RT' % response times
        for d = 1:length(D_act)

            % subject's responses
            try
                actual_rt = D_act(d).HGF(1).yALL(:,2);
            catch
                actual_rt = D_act(d).HGF(1).y(:,2);
            end

            % outliers
            [~,actualOL] = rmoutliers(actual_rt);

            % condition means
            condnum = D_act(d).Processed.condnum{1};
            ii = ~isnan(actual_rt);
            unicond = unique(condnum(ii));
            actual_condmeans(d,unicond) = mean_by_condition(actual_rt(ii),condnum(ii)); 

            % split into stim intensity per block
            blocknum = D_act(d).Sequence.blocks;
            stimnum = D_act(d).Sequence.signal(S.signal.target,:);
            [~,~,blockstimnum] = unique([blocknum;stimnum]','rows');
            uniblockstim = unique(blockstimnum);
            actual_blockstimmeans(d,uniblockstim) = mean_by_condition(actual_rt(ii),blockstimnum(ii)); 

            fitted_condmeans_rep = [];
            fitted_blockstimmeans_rep = [];
            fitted_rt_all=[];
            fitted_choices_rep = [];
            cc=[];
            clear brr
            disp(['BRR on RTs, ppt ' num2str(d)])
            tic
            HGF = D_sim(d).HGF;
            parfor rep = 1:S.numsimrep
                fitted_rt = HGF(rep).sim.y(:,2);
                fitted_choices = HGF(rep).sim.y(:,1);

                if S.run_BRR && S.BRR_on_each_rep

                    % outliers
                    [~,fittedOL] = rmoutliers(fitted_rt);
                   
                    % NaNs
                    ii = ~isnan(actual_rt) & ~isnan(fitted_rt) & ~actualOL & ~fittedOL;
    %                 ii = ~isnan(actual_rt) & ~isnan(fitted_rt);
    
                    fitted_rt(~ii)=nan;

                    % remove choice effects (maybe only needed if jointly
                    % modelling RT and choice in the HGF?)
                    mdlF=fitlm(fitted_choices(ii),fitted_rt(ii));
                    fitted_rt(ii) = mdlF.Residuals.Raw;

                    % bayes reg
                    brr(d,rep) = bayesreg_crossval(fitted_rt(ii),actual_rt(ii),S,{});
                    %brr(d,rep).N = sum(ii);
                else
                    ii = ~isnan(actual_rt) & ~isnan(fitted_rt);
                    fitted_rt(~ii)=nan;
                end
                fitted_rt_all(:,rep) = fitted_rt;
                fitted_choices_rep(:,rep) = fitted_choices;

                %correlate
                cc(rep) = corr(fitted_rt(ii),actual_rt(ii),'type','Spearman');
                
                % condition means
                fitted_condmeans_rep(rep,:) = mean_by_condition(fitted_rt ,condnum); 
                
                % blockstim means
                fitted_blockstimmeans_rep(rep,:) = mean_by_condition(fitted_rt ,blockstimnum); 
               
                
            end
            toc
            cc_mean(d,1)=mean(cc);

            % RTs
            fitted_rt_mean = nanmean(fitted_rt_all,2);
            [~,fittedOL] = rmoutliers(fitted_rt_mean);

            if S.run_BRR
                if S.BRR_on_each_rep
                    fdnames = fieldnames(brr);
                    for fd = 1:length(fdnames)
                        brr_mean.(fdnames{fd})(d,1)=mean([brr(:).(fdnames{fd})]);
                    end
                else
    
    
                    % remove choice effects (maybe only needed if jointly
                    % modelling RT and choice in the HGF?)
    
                    fitted_choices_mean = nanmean(fitted_choices_rep,2);
                    ii = ~isnan(fitted_choices_mean) & ~isnan(actual_rt) & ~isnan(fitted_rt_mean) & ~actualOL & ~fittedOL;
    %                 ii = ~isnan(fitted_choices_mean) & ~isnan(actual_rt) & ~isnan(fitted_rt_mean);
                    mdlF=fitlm(fitted_choices_mean(ii),fitted_rt_mean(ii));
                    fitted_rt_mean(~ii)=nan;
                    fitted_rt_mean(ii) = mdlF.Residuals.Raw;
    
                    % bayes reg
                    brr_mean(d,1) = bayesreg_crossval(fitted_rt_mean(ii),actual_rt(ii),S,{});
                    %brr_mean(d,1).N = sum(ii);
                end
            end

            %[~,fittedOL] = rmoutliers(fitted_rt_mean);
            %fitted_rt_mean(fittedOL)=nan;

            fitted_condmeans(d,:) = mean(fitted_condmeans_rep,1);
            fitted_blockstimmeans(d,:) = mean(fitted_blockstimmeans_rep,1);
            actual_rt_allsub(d,:) = actual_rt';
            fitted_rt_allsub(d,:) = fitted_rt_mean';
        end
        if ~exist("brr_mean","var") brr_mean=struct; end
        varargout = {cc_mean,brr_mean,actual_condmeans,fitted_condmeans,actual_blockstimmeans,fitted_blockstimmeans,actual_rt_allsub,fitted_rt_allsub};
        
    case {'HGFvar','CCA_FApred'} % decoded HGF trajectories
        for d = 1:length(D_act)
            actual_var(d,:,:) = D_act(d).HGF(1).yALL(:,3:end);
            for rep = 1:S.numsimrep
                sim_var = D_sim(d).HGF(rep).sim.y(:,3:end);
                
                nvar = size(sim_var,2);
                
                for nv=1:nvar

                    %only consider trials in which both fitted and actual
                    %responses both occur
                    ii = ~isnan(squeeze(actual_var(d,:,nv))') .* ~isnan(sim_var(:,nv));

                    %correlate
                    cc(nv,rep) = corr(sim_var(ii,nv),squeeze(actual_var(d,ii,nv))','type','Spearman');

                    % bayes reg
                    brr(d,nv,rep) = bayesreg_crossval(sim_var(ii,nv),squeeze(actual_var(d,ii,nv))',S,{});
                    waic(d,rep,nv) = brr(d,nv,rep).waic;
                end
            end
            cc_mean(d,:)=mean(cc,2);
        end
        cc_varmean=mean(cc_mean,2);
        waic = mean(waic,3);
        varargout = {actual_var,cc_varmean,brr,waic};
        
        
end

function condmeans = mean_by_condition(data,condnum) % condnum = D_fit(d).Processed.condnum{1};  

% split into conditions
condsuni = unique(condnum);
condsuni = condsuni(~isnan(condsuni));
condmeans = nan(1,length(condsuni));
for i = 1:length(condsuni)
    condmeans(i) = nanmean(data(condnum==condsuni(i))); 
end
