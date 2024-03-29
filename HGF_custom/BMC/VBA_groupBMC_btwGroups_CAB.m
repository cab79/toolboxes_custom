function [h,p,varargout] = VBA_groupBMC_btwGroups_CAB(Ls,options,pep_flag)
% test for between-groups difference in model frequencies
% function [h,p] = VBA_groupBMC_btwGroups(Ls,options)
% IN:
%   - Ls: {nmXns_1, nmXns_2} array of log-model evidences matrices of each group (nm models; ns_g subjects in the group).
%   - options: a structure containing the following fields:
%       .DisplayWin: flag for display window
%       .verbose: flag for summary statistics display
%       .families: a cell array of size nf, which contains the indices of
%       the models that belong to each of the nf families.
% OUT:
%   - h: test decision about a difference between the group (rejection of
%        the null hypothesis of equality)
%   - p: the posterior probability that the two groups have the same model
%        frequencies

if nargin < 2
    options = {};
end
% one frequency for all
L = horzcat(Ls{:});
options.figName = 'RFX-BMS: all subjects';
[posterior,out] = VBA_groupBMC_cab(L,options,pep_flag);
Fe = out.F(end);
varargout = {out};

% separate frequencies
Fd = 0;
for g = 1:length(Ls)
    options.figName = ['RFX-BMS: group ' options.grpNames{g}];
    [~,out] = VBA_groupBMC_cab(Ls{g},options,pep_flag);
    Fd = Fd+out.F(end);
    varargout = [varargout,{out}];
end

if length(varargout)>1
    varargout = {varargout};
end
% if length(Ls)==2
%     [posterior1,out1] = VBA_groupBMC_cab(Ls{1},options,pep_flag);
%     [posterior2,out2] = VBA_groupBMC_cab(Ls{2},options,pep_flag);
%     Fd = out1.F(end) + out2.F(end);
%     varargout = {[varargout,{out1},{out2}]};
% elseif length(Ls)==3
%     [posterior1,out1] = VBA_groupBMC_cab(Ls{1},options,pep_flag);
%     [posterior2,out2] = VBA_groupBMC_cab(Ls{2},options,pep_flag);
%     [posterior3,out3] = VBA_groupBMC_cab(Ls{3},options,pep_flag);
%     Fd = out1.F(end) + out2.F(end) + out3.F(end);
%     varargout = {[varargout,{out1},{out2},{out3}]};
% end

p = 1./(1+exp(Fd-Fe));
h = p<.05;

end