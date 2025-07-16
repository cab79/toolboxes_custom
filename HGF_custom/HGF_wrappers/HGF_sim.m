function [D,S] = HGF_sim(D,S)

if S.parallel
    checkp = gcp('nocreate');
    if isempty(checkp)
        myPool = parpool;
    end
    parforArg = Inf;
else
    parforArg = 0;
end

for d = 1:length(D) 

    S.HGF.selectrep=1;
    HGF=D(d).HGF;
    D_in=D(d);

    tic
    if S.parallel
        parfor (ns = 1:S.numsimrep,parforArg)
        % for ns = 1:S.numsimrep
            %
            % if S.parallel
            %    D_out=HGF_run(D(d),S,1);
            % else
                D_out=HGF_run_nopara(D_in,S,1); % takes longer using parallel
            % end
            HGF(ns).sim = D_out.HGF(1).sim;
        end
    else
        for ns = 1:S.numsimrep
            %
            %if S.parallel
            %    D_out=HGF_run(D(d),S,1);
            %else
                D_out=HGF_run_nopara(D_in,S,1); % takes longer using parallel
            %end
            HGF(ns).sim = D_out.HGF(1).sim;
        end
    end
    D(d).HGF = HGF;
    toc
    % for binary choices
    for ns = 1:S.numsimrep
        if length(D(d).HGF(ns).sim.y)==length(D(d).HGF(1).u)
            D(d).Output(ns).presstrial = 1:length(D(d).HGF(ns).sim.y);
            try 
                buttonopt = D(d).Output(1).Settings.buttonopt;
            catch
                buttonopt={};
            end
            if length(unique(D(d).HGF(ns).sim.y+1))==2 % should be binary, otherwise may be RT
                if ~isempty(buttonopt)
                    D(d).Output(ns).pressbutton = buttonopt(D(d).HGF(ns).sim.y+1);
                else
                    D(d).Output(ns).pressbutton = D(d).HGF(ns).sim.y;
                end
            end
        else
            error('simulated data has the wrong number of trials')
        end
    end
end