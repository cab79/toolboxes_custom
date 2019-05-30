% zest_scratchpad

clear init

% params for initial p.d.f.
init.zestA = 1;
init.zestB = 2.5; % = .025; % B= 2.5 for Marvit et al., 2003
init.zestC = 2.5; % = .1; %C = 2.5 for Marvit et al., 2003;
% CHANGE BACK TO LOG
init.zestmaxrange = log10(20); % log of highest threshold value possible; 
init.zestminrange = log10(.01); % log of lowest threshold value possible; 

% Parameters for Response function (Weibull function) (P(yes) given stimulus)
init.zestfa = 0.1;   %gamma in the text, false alarm rate (guess rate for 2AFC)
init.zestmiss = 0.02; %delta in the text, miss rate (1/2 inattention rate for 2AFC)
init.zestbeta = 6; %10;    %beta in the text, slope of response function
init.zesteta = 0;     %eta in the text, "sweat factor" or response criterion parameter
	
% Starting params
init.zestinit_diffLvl = 3; %10; % initial difference level used for Fig 1A of Marvit et al.	was 3 db

% UNCOMMENT IF USING LOG
init.zestconvert = {'delta_L', 'sd_pdf'};

max_trials = 20;
thresh_tol = .01;

% known_DL = 1;

% % read in data table
% [testSub_numData testSub_txtData testSub_rawData] = xlsread('/data/attmap/pilot/attmap_pilot_testsub.xls');
% 
% % get column indices
% subResp_col_idx = find(ismember(testSub_rawData(1,:),'sub_resp'));
% obsNum_col_idx = find(ismember(testSub_rawData(1,:),'obs_num'));
% meanPDF_col_idx = find(ismember(testSub_rawData(1,:),'mean_pdf'));
% converged_col_idx = find(ismember(testSub_rawData(1,:),'converged'));
% 
% % get data vectors
% subResp_data = cell2mat(testSub_rawData(2:end,subResp_col_idx));
% obsNum_data = cell2mat(testSub_rawData(2:end,obsNum_col_idx));
% meanPDF_data = cell2mat(testSub_rawData(2:end,meanPDF_col_idx));
% converged_data = cell2mat(testSub_rawData(2:end,converged_col_idx));
% 
clear resps
clear thresh

subResp_data = [1 0 0 1 1 0 1 0 1 0];

% thresh = nan(size(subResp_data));
% initialize p.d.f.
ZEST_marvit(NaN,init);

for iresp = 1:length(subResp_data)
    thresh(iresp) = ZEST_marvit(subResp_data(iresp));
%     thresh(iresp) = ZEST_marvit(subResp_data(iresp));
%     if thresh(iresp)<known_DL
%         next_resp = 0;
%     else
%         next_resp = 1;
%     end
%     resps(iresp+1) = next_resp;
    
    if iresp > 1 && abs(thresh(iresp)-thresh(iresp-1)) <= thresh_tol
        final_thresh = thresh(iresp);
        sprintf('ZEST converged at %2.4f dB in %d trials',final_thresh,iresp)
        break
    elseif iresp == max_trials
        final_thresh = thresh(iresp);
        sprintf('ZEST reached maximum number of trials prior to convergence')
    end
end

sprintf('end of script')
