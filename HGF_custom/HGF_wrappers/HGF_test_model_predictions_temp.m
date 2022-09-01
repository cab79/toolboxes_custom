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
fitted_condmeans=nan(S.numsimrep,length(allcond));

switch type
    case 'Ch' % choices
        actual_macorrect=[];
        fitted_macorrect=[];
        for d = 1:length(D_act)
            actual_choices = D_act(d).HGF(1).y(:,1);

            % condition fraction correct
            condnum = D_act(d).Processed.condnum{1};
            ii = find(~isnan(actual_choices));
            unicond = unique(condnum(ii));
            actual_condmeans(d,unicond) = mean_by_condition(actual_choices(ii)==D_act(d).HGF(1).u(ii,1),condnum(ii)); 

            % split into stim intensity per block
            blocknum = D(d).Sequence.blocks;
            stimnum = D(d).Sequence.signal(S.signal.target,:);
            [uniblockstim,~,blockstimnum] = unique([blocknum;stimnum]','rows');
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

            fitted_condmeans_rep = [];
            fitted_blockstimmeans_rep = [];
            fitted_correct = [];
            rand_correct = [];
            for rep = 1:S.numsimrep
                fitted_choices = D_sim(d).HGF(rep).sim.y(:,1);
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
                temp_macorrect=[];
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
            means_fitted(d,1)=mean(fitted_correct);
            stds_fitted(d,1)=std(fitted_correct);
            fitted_macorrect(d,:)=mean(temp_macorrect,1);
            % random-corrected
            means_randcorr_fitted(d,1)=mean(fitted_correct-rand_correct);
            stds_randcorr_fitted(d,1)=std(fitted_correct-rand_correct);
            fitted_condmeans(d,:) = mean(fitted_condmeans_rep,1);
            fitted_blockstimmeans(d,:) = mean(fitted_blockstimmeans_rep,1);
        end
        varargout = {means_fitted,stds_fitted,means_randcorr_fitted,stds_randcorr_fitted,actual_condmeans,fitted_condmeans,actual_macorrect,fitted_macorrect};

    case 'RT' % response times
        for d = 1:length(D_act)

            % subject's responses
            actual_rt = D_act(d).HGF(1).y(:,2);

            % condition means
            condnum = D_act(d).Processed.condnum{1};
            ii = find(~isnan(actual_rt));
            unicond = unique(condnum(ii));
            actual_condmeans(d,unicond) = mean_by_condition(actual_rt(ii),condnum(ii)); 


            fitted_condmeans_rep = [];
            for rep = 1:S.numsimrep
                fitted_rt = D_sim(d).HGF(rep).sim.y(:,2);

                % NaNs
                ii = find(~isnan(actual_rt) & ~isnan(fitted_rt));

                % condition means
                fitted_condmeans_rep(rep,:) = mean_by_condition(fitted_rt ,condnum); 
                
                %correlate
                cc(rep) = corr(fitted_rt(ii),actual_rt(ii),'type','Spearman');
                
                % bayes reg
                brr(d,rep) = bayesreg_crossval(fitted_rt(ii),actual_rt(ii),S,{});
                
            end
            cc_mean(d,1)=mean(cc);
            fdnames = fieldnames(brr);
            for fd = 1:length(fdnames)
                brr_mean.(fdnames{fd})(d,1)=mean([brr(:).(fdnames{fd})]);
            end
            fitted_condmeans(d,:) = mean(fitted_condmeans_rep,1);
        end
        varargout = {cc_mean,brr_mean,actual_condmeans,fitted_condmeans};
        
    case {'HGFvar','CCA_FApred'} % decoded HGF trajectories
        for d = 1:length(D_act)
            actual_var(d,:,:) = D_act(d).HGF(1).y(:,3:end);
            for rep = 1:S.numsimrep
                sim_var = D_sim(d).HGF(rep).sim.y(:,3:end);
                
                nvar = size(sim_var,2);
                
                for nv=1:nvar

                    %only consider trials in which both fitted and actual
                    %responses both occur
                    ii = find(~isnan(squeeze(actual_var(d,:,nv))') .* ~isnan(sim_var(:,nv)));

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
