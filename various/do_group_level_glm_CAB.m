function correlationMats = do_group_level_glm_CAB(correlationMats, Settings)
%DO_GROUP_LEVEL_GLM run group comparison 
% 


%	Copyright 2014 OHBA
%	This program is free software: you can redistribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation, either version 3 of the License, or
%	(at your option) any later version.
%	
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%	
%	You should have received a copy of the GNU General Public License
%	along with this program.  If not, see <http://www.gnu.org/licenses/>.


%	$LastChangedBy: giles.colclough@gmail.com $
%	$Revision: 231 $
%	$LastChangedDate: 2014-08-07 20:53:06 +0100 (Thu, 07 Aug 2014) $
%	Contact: giles.colclough@eng.ox.ac.uk
%	Originally written on: MACI64 by Giles Colclough, 10-Apr-2014 13:33:21

groupLevel = struct;
for iFreq = Settings.nFreqBands:-1:1,
if strcmpi(Settings.paradigm, 'rest'),
    % univariate edge testing in turn for correlation, partial correlation
    % and regularized partials.
    [T, p, fweptmp, fdrptmp, COPEmat] = univariate_edge_test(correlationMats{iFreq}.envCorrelation_z, ...
                                                       Settings.GroupLevel.designMatrix,        ...
                                                       Settings.GroupLevel.contrasts,[],Settings);
    %for iContrast = size(p,3):-1:1,
    %    h(:,:,iContrast) = false_discovery_rate(p_to_z_two_tailed(p(:,:,iContrast)), ...
    %                                                    Settings.FDRalpha);
    %end%for
    
    groupLevel.correlation.T    = T;
    groupLevel.correlation.p    = p;
    groupLevel.correlation.FWEp = fweptmp;
    groupLevel.correlation.FDRp = fdrptmp;
    groupLevel.correlation.COPE = COPEmat;
    
    [T, p, fweptmp, fdrptmp, COPEmat] = univariate_edge_test(correlationMats{iFreq}.envPartialCorrelation_z, ...
                                                 Settings.GroupLevel.designMatrix,        ...
                                                 Settings.GroupLevel.contrasts,[],Settings);
    %for iContrast = size(p,3):-1:1,
    %    h(:,:,iContrast) = false_discovery_rate(p_to_z_two_tailed(p(:,:,iContrast)), ...
    %                                                    Settings.FDRalpha);
    %end%for
    
    groupLevel.partialCorrelation.T    = T;
    groupLevel.partialCorrelation.p    = p;
    groupLevel.partialCorrelation.FWEp = fweptmp;
    groupLevel.partialCorrelation.FDRp = fdrptmp;
    groupLevel.partialCorrelation.COPE = COPEmat;
    
    [T, p, fweptmp, fdrptmp, COPEmat] = univariate_edge_test(correlationMats{iFreq}.envPartialCorrelationRegularized_z, ...
                                                 Settings.GroupLevel.designMatrix,        ...
                                                 Settings.GroupLevel.contrasts,[],Settings);
    %for iContrast = size(p,3):-1:1,
    %    h(:,:,iContrast) = false_discovery_rate(p_to_z_two_tailed(p(:,:,iContrast)), ...
    %                                                    Settings.FDRalpha);
    %end%for
    
    groupLevel.partialCorrelationRegularized.T    = T;
    groupLevel.partialCorrelationRegularized.p    = p;
    groupLevel.partialCorrelationRegularized.FWEp = fweptmp;
    groupLevel.partialCorrelationRegularized.FDRp = fdrptmp;
    groupLevel.partialCorrelationRegularized.COPE = COPEmat;
    
    % add in design for reference
    groupLevel.designMatrix = Settings.GroupLevel.designMatrix;
    groupLevel.contrasts    = Settings.GroupLevel.contrasts;
    
    
elseif strcmpi(Settings.paradigm, 'task'),
    % we need to run a separate GLM for each first level contrast. 
    % And for each group-level contrast.
    % What a mare.
    for iContrast = length(Settings.SubjectLevel.contrasts):-1:1,
		% we use parameter estimates from the level below for each of
		% correlation, partial correlation and regularised partial
		% correlation
        Settings.GroupLevel.iInter=0; % NOT interaction on within-subject factors
		if isfield(correlationMats{iFreq}, 'subjectLevel'),
			COPE = correlationMats{iFreq}.subjectLevel(iContrast).cope;
		elseif isfield(correlationMats{iFreq}, 'firstLevel'),
			COPE = correlationMats{iFreq}.firstLevel(iContrast).cope;
		else
			error([mfilename ':WhereIsTheData'], ...
				  'Expected input to have either first level or subject level results. \n');
		end%if
        
        groupLevel = correlation_stats(groupLevel, Settings, COPE, iContrast);
		
    end%for
    
    % within-subject interactions
    for iInter = size(Settings.SubjectLevel.interaction,1):-1:1,
        Settings.GroupLevel.iInter=iInter; % interaction on within-subject factors
        iContrast = length(Settings.SubjectLevel.contrasts)+iInter:length(Settings.SubjectLevel.contrasts)+iInter+1;
		% we use parameter estimates from the level below for each of
		% correlation, partial correlation and regularised partial
		% correlation
        COPE = cell(2,1);
        for ic = 1:length(iContrast)
            if isfield(correlationMats{iFreq}, 'subjectLevel'),
                COPE{ic} = correlationMats{iFreq}.subjectLevel(iContrast(ic)).cope;
            elseif isfield(correlationMats{iFreq}, 'firstLevel'),
                COPE{ic} = correlationMats{iFreq}.firstLevel(iContrast(ic)).cope;
            else
                error([mfilename ':WhereIsTheData'], ...
                      'Expected input to have either first level or subject level results. \n');
            end%if
        end
        
        groupLevel = correlation_stats(groupLevel, Settings, COPE, iContrast);
		
    end%for
    
else
    error([mfilename 'BadParadigm'], ...
          'Unrecognised paradigm %s. \n', Settings.paradigm);
end%if

correlationMats{iFreq}.groupLevel = groupLevel;
end%for loop over frequencies

end % for function

function groupLevel = correlation_stats(groupLevel, Settings, COPE, iContrast)

if iscell(COPE) && length(iContrast)>1 % interaction
    COPEcell = COPE;
    COPE=struct;
    for i = 1:length(COPEcell)
        COPE.correlation(:,:,:,i) = COPEcell{i}.correlation;
        COPE.partialCorrelation(:,:,:,i) = COPEcell{i}.partialCorrelation;
        COPE.partialCorrelationRegularized(:,:,:,i) = COPEcell{i}.partialCorrelationRegularized;
    end
    iContrast=iContrast(1);
end
%-- correlation
[T, p, fweptmp, fdrptmp, COPEmat] = univariate_edge_test(COPE.correlation,                 ...
                                                   Settings.GroupLevel.designMatrix, ...
                                                   Settings.GroupLevel.contrasts,[],Settings);
%for iConG = size(p,3):-1:1,
%    h(:,:,iConG) = false_discovery_rate(p_to_z_two_tailed(p(:,:,iConG)), ...
%                                                Settings.FDRalpha);
%end%for

groupLevel(iContrast).correlation.T    = T;
groupLevel(iContrast).correlation.p    = p;
groupLevel(iContrast).correlation.FWEp = fweptmp;
groupLevel(iContrast).correlation.FDRp = fdrptmp;
groupLevel(iContrast).correlation.COPE = COPEmat;

%-- partial correlation
[T, p, fweptmp, fdrptmp, COPEmat] = univariate_edge_test(COPE.partialCorrelation,          ...
                                                   Settings.GroupLevel.designMatrix, ...
                                                   Settings.GroupLevel.contrasts,[],Settings);
%for iConG = size(p,3):-1:1,
%    h(:,:,iConG) = false_discovery_rate(p_to_z_two_tailed(p(:,:,iConG)), ...
%                                                Settings.FDRalpha);
%end%for

groupLevel(iContrast).partialCorrelation.T    = T;
groupLevel(iContrast).partialCorrelation.p    = p;
groupLevel(iContrast).partialCorrelation.FWEp = fweptmp;
groupLevel(iContrast).partialCorrelation.FDRp = fdrptmp;
groupLevel(iContrast).partialCorrelation.COPE = COPEmat;

%-- regularised partial correlation
if Settings.Regularize.do,
    [T, p, fweptmp, fdrptmp, COPEmat] = univariate_edge_test(COPE.partialCorrelationRegularized, ...
                                                       Settings.GroupLevel.designMatrix,   ...
                                                       Settings.GroupLevel.contrasts,[],Settings);
    %for iConG = size(p,3):-1:1,
    %    h(:,:,iConG) = false_discovery_rate(p_to_z_two_tailed(p(:,:,iConG)), ...
    %                                                Settings.FDRalpha);
    %end%for

    groupLevel(iContrast).partialCorrelationRegularized.T    = T;
    groupLevel(iContrast).partialCorrelationRegularized.p    = p;
    groupLevel(iContrast).partialCorrelationRegularized.FWEp = fweptmp;
    groupLevel(iContrast).partialCorrelationRegularized.FDRp = fdrptmp;
    groupLevel(iContrast).partialCorrelationRegularized.COPE = COPEmat;
end%if
% add in first levels for reference
if iContrast<=length(Settings.SubjectLevel.contrasts)
groupLevel(iContrast).firstLevelContrast        = Settings.SubjectLevel.contrasts{iContrast};
groupLevel(iContrast).firstLevelConditionLabels = Settings.SubjectLevel.conditionLabel;
end

% add in design for reference
groupLevel(iContrast).groupDesignMatrix = Settings.GroupLevel.designMatrix;
groupLevel(iContrast).groupContrasts    = Settings.GroupLevel.contrasts;
end