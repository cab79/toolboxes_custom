function c = response_model_config(r, S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Contains the configuration for the linear log-reaction time response model according to as
% developed with Louise Marshall and Sven Bestmann
% http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1002575
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The Gaussian noise observation model assumes that responses have a Gaussian distribution around
% the inferred mean of the relevant state. The only parameter of the model is the noise variance
% (NOT standard deviation) zeta.

% Config structure
c = struct;

% Perceptual model used
c.responses = S.resp_modelspec.responses;
if iscell(r.c_prc.response.priormodel)
    if strcmp(c.responses,'Ch')
        c.model = r.c_prc.response.priormodel{1};
    elseif strcmp(c.responses,'RT') || any(strcmp(c.responses, 'HGFvar'))
        c.model = r.c_prc.response.priormodel{2};
    end
else
    c.model = r.c_prc.response.priormodel;
end

% CAB: Number of levels
try
    l = r.c_prc.(c.model).n_priorlevels+1;
catch
    l=1;
end

c.nparams =[];
c.priormus=[];
c.priorsas=[];
c.st = [];
c.pn=0;

% for prediction error, use absolute values?
try
    c.abs = S.resp_modelspec.PE_abs;
catch
    c.abs = 1;
end

% for simulation
try 
    c.deterministic = S.deterministic;
catch
    c.deterministic=0;
end

if any(strcmp(c.responses , 'Ch'))
    % USE SOFTMAX MODEL
    
    % Decision based on which representation?
    c.rep = r.c_prc.response.rep; 
    
    % Beta
    c.soft.logbemu = log(48);
    c.soft.logbesa = 1;
    
    % Gather prior settings in vectors
    type='soft';
    c = paramvec(c,type);
end


if any(strcmp(c.responses, 'RT')) || any(strcmp(c.responses, 'EEG')) || any(strcmp(c.responses, 'HGFvar')) || any(strcmp(c.responses, 'CCAcomp'))  || any(strcmp(c.responses, 'CCA_FApred'))
    % USE REGRESSION MODEL
    
    c.params = S.resp_modelspec.params;
    
    % Sufficient statistics of Gaussian parameter priors
    
    if any(strcmp(c.responses, 'RT'))
        
        % Beta_0
        c.reg.be0mu = log(0.5); 
        c.reg.be0sa = 4;
    elseif any(strcmp(c.responses, 'EEG')) 
        
        % Beta_0
        c.reg.be0mu = 0; 
        c.reg.be0sa = 4;
    elseif any(strcmp(c.responses, 'HGFvar'))

        % variable names
        c.ynames = S.resp_modelspec.ynames;
        c.ynum = S.resp_modelspec.responsenum;
        if ~isempty(c.ynum)
            c.ynames=c.ynames(c.ynum);
        end

        for yn = 1:length(c.ynames)
            % Beta_0
            c.reg.(['be0' num2str(c.params(yn)) 'mu']) = 0; 
            c.reg.(['be0' num2str(c.params(yn)) 'sa']) = 4;
        end
    elseif any(strcmp(c.responses, 'CCAcomp'))

        if ~iscell(c.params)
            error("parameter numbers should be in a cell array for CCAcomp")
        end
        
        for yn = 1:length(c.params)
            % Beta_0
            c.reg.(['be0' num2str(yn) 'mu']) = 0; 
            c.reg.(['be0' num2str(yn) 'sa']) = 4;
        end
        
        c.params_cell = c.params;
        c.params = unique([c.params{:}]);
    
    elseif any(strcmp(c.responses, 'CCA_FApred'))

        if ~iscell(c.params)
            error("parameter numbers should be in a cell array for CCAcomp")
        end
        
        % variable names and weights
        c.ynames = S.resp_modelspec.ynames;
        c.weights = S.resp_modelspec.weights;
        c.factors = S.resp_modelspec.factors;
%         c.ynum = S.resp_modelspec.responsenum;
%         if ~isempty(c.ynum)
%             c.ynames=c.ynames(c.ynum);
%             c.weights=c.weights(c.ynum);
%         end
        
        for y = 1:length(c.factors)
            % Beta_0
            c.reg.(['be0' num2str(y) 'mu']) = 0; 
            c.reg.(['be0' num2str(y) 'sa']) = 4;
            % Zeta
            c.reg.(['ze' num2str(y) 'mu']) = log(log(20)); 
            c.reg.(['ze' num2str(y) 'sa']) = log(2);
            c.reg.(['ze' num2str(y) 'var']) = true;
        end
        
        c.params_cell = c.params;
        c.params = unique([c.params{:}]);
        
        for yn = 1:length(c.params)
            c.reg.(['be' num2str(c.params(yn)) 'mu']) = 0; 
            c.reg.(['be' num2str(c.params(yn)) 'sa']) = 4;
        end
    end
    
    if ~any(strcmp(c.responses, 'CCA_FApred'))

        % Beta_1
        if any(c.params==1)
            c.reg.be1mu = 0;
            c.reg.be1sa = 4;
        else
            c.reg.be1mu = 0;
            c.reg.be1sa = 0;
        end

        % Beta_2
        if any(c.params==2)
            c.reg.be2mu = 0;
            c.reg.be2sa = 4;
        else
            c.reg.be2mu = 0;
            c.reg.be2sa = 0;
        end

        % Beta_3
        if any(c.params==3)
            if l>1
                c.reg.be3mu = 0;
                c.reg.be3sa = 4;
            else
                c.reg.be3mu = 0;
                c.reg.be3sa = 0;
            end
        else
            c.reg.be3mu = 0;
            c.reg.be3sa = 0;
        end

        % Beta_4
        if any(c.params==4)
            if l>2
                c.reg.be4mu = 0;
                c.reg.be4sa = 4;
            else
                c.reg.be4mu = 0;
                c.reg.be4sa = 0;
            end
        else
            c.reg.be4mu = 0;
            c.reg.be4sa = 0;
        end

        % Beta_5
        if any(c.params==5)
            c.reg.be5mu = 0;
            c.reg.be5sa = 4;
        else
            c.reg.be5mu = 0;
            c.reg.be5sa = 0;
        end

        % Beta_6
        if any(c.params==6)
            c.reg.be6mu = 0;
            c.reg.be6sa = 4;
        else
            c.reg.be6mu = 0;
            c.reg.be6sa = 0;
        end

        % Beta_7
        if any(c.params==7)
            c.reg.be7mu = 0;
            c.reg.be7sa = 4;
        else
            c.reg.be7mu = 0;
            c.reg.be7sa = 0;
        end

        % Beta_8
        if any(c.params==8)
            c.reg.be8mu = 0;
            c.reg.be8sa = 4;
        else
            c.reg.be8mu = 0;
            c.reg.be8sa = 0;
        end

        % Beta_9
        if any(c.params==9)
            if l>2
                c.reg.be9mu = 0;
                c.reg.be9sa = 4;
            else
                c.reg.be9mu = 0;
                c.reg.be9sa = 0;
            end
        else
            c.reg.be9mu = 0;
            c.reg.be9sa = 0;
        end

        % Beta_10
        if any(c.params==10)
            if l>2
                c.reg.be10mu = 0;
                c.reg.be10sa = 4;
            else
                c.reg.be10mu = 0;
                c.reg.be10sa = 0;
            end
        else
            c.reg.be10mu = 0;
            c.reg.be10sa = 0;
        end


        % Beta_11
        if any(c.params==11)
            if l>2
                c.reg.be11mu = 0;
                c.reg.be11sa = 4;
            else
                c.reg.be11mu = 0;
                c.reg.be11sa = 0;
            end
        else
            c.reg.be11mu = 0;
            c.reg.be11sa = 0;
        end

        % Beta_12
        if any(c.params==12)
            c.reg.be12mu = 0;
            c.reg.be12sa = 4;
        else
            c.reg.be12mu = 0;
            c.reg.be12sa = 0;
        end

        % Beta_13
        if any(c.params==13)
            if l>1
                c.reg.be13mu = 0;
                c.reg.be13sa = 4;
            else
                c.reg.be13mu = 0;
                c.reg.be13sa = 0;
            end
        else
            c.reg.be13mu = 0;
            c.reg.be13sa = 0;
        end

        % Beta_14
        if any(c.params==14)
            if l>2
                c.reg.be14mu = 0;
                c.reg.be14sa = 4;
            else
                c.reg.be14mu = 0;
                c.reg.be14sa = 0;
            end
        else
            c.reg.be14mu = 0;
            c.reg.be14sa = 0;
        end

        % Beta_15
        if any(c.params==15)
            c.reg.be15mu = 0;
            c.reg.be15sa = 4;
        else
            c.reg.be15mu = 0;
            c.reg.be15sa = 0;
        end


        % Beta_16
        if any(c.params==16)
            if l>1
                c.reg.be16mu = 0;
                c.reg.be16sa = 4;
            else
                c.reg.be16mu = 0;
                c.reg.be16sa = 0;
            end
        else
            c.reg.be16mu = 0;
            c.reg.be16sa = 0;
        end


        % Beta_17
        if any(c.params==17)
            if l>2
                c.reg.be17mu = 0;
                c.reg.be17sa = 4;
            else
                c.reg.be17mu = 0;
                c.reg.be17sa = 0;
            end
        else
            c.reg.be17mu = 0;
            c.reg.be17sa = 0;
        end
        
        % Zeta
        c.reg.logzemu = log(log(20));
        c.reg.logzesa = log(2);
        c.reg.logzevar = true; % this is a variance parameter
        
    end
    
    
    % Gather prior settings in vectors
    type='reg';
    c = paramvec(c,type);
end

% Model filehandle
c.obs_fun = @response_model;

% Handle to function that transforms observation parameters to their native space
% from the space they are estimated in
c.transp_obs_fun = @response_model_transp;

return;

function c = paramvec(c,type)
fn=fieldnames(c.(type));
for i = 1:length(fn)
    if strcmp(fn{i}(end-1:end),'mu')
        c.pn=c.pn+1;
        c.pnames{c.pn,1} = [type '_' fn{i}(1:end-2)];
        nme_gen = strsplit(fn{i}(1:end-2),'log');
        c.pnames_gen{c.pn,1} = nme_gen{end};
        c.pnames_mod{c.pn,1} = [type '_' nme_gen{end}];
        eval(['c.priormus = [c.priormus c.(type).' fn{i} '];']);
        eval(['c.nparams(c.pn) = length(c.(type).' fn{i} ');']);
        if isfield(c.(type),[fn{i}(1:end-2) 'var'])
            c.varparam(c.pn)=1;
        else
            c.varparam(c.pn)=0;
        end
    elseif strcmp(fn{i}(end-1:end),'sa')
        eval(['c.priorsas = [c.priorsas c.(type).' fn{i} '];']);
    else
        continue
    end
    if isempty(c.st)
        c.st = 0;
    else
        c.st=sum(c.nparams(1:c.pn-1));
    end
    c.priormusi{c.pn} = c.st+1:sum(c.nparams(1:c.pn));
end
