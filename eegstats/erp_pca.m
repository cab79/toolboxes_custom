function [FactorResults] = erp_pca(data, NUM_FAC, type, varargin)

% CAB NEW: varargin inputs: erp_pca(data, NUM_FAC, degree, S, Sb)
% S = covariance matrix
% Sb = for generalised eigenvector problem

% Full description: https://www.tandfonline.com/doi/full/10.1080/87565641.2012.697503
% This version implement promax with Kaiser rotation on covariance matrix

% ep_doPCA - [FactorResults] = ep_doPCA(PCAmode, ROTATION, ROTOPT, MAT_TYPE, NUM_FAC, theData, LOADING, GAVE, badData, crossFacCof) -
%          do PCA and calculate factor scores
%Inputs
%  PCAmode	 : Primary mode for the PCA ('temp': time points as the variables; 'spat': channels as the variables; 'freq': frequencies as variables; 'asis': do not rearrange input data).
%  ROTATION	:
%       PCAs and has various options for the type of rotation, including
%       UNRT - Unrotated
%       VMAX - Varimax
%       PMAX - Promax (kappa=3 seems to give the best results for ERPs)
%       IMAX - ICA Infomax (if EEGLab is also installed.  ICA is run with
%           the PCA subspace option on by default.  This program can be edited
%           to change its options.  Output is rescaled to be consistent with PCA conventions)
%       QMAX - Quartimax
%       QMIN - Quartimin
%       OMIN - Oblimin (gamma=0 tends to be recommended)
%       CRFE - Crawford-Ferguson family
%       MINE - minimum entropy
%       IPSC - Bentler's invariant pattern simplicity criterion
%       TIIC - Comrey's tandem II criterion
%       GMIN - geomin
%       MMER - McCammon minimum entropy ratio
%       VOMN - Variable Oblimin (tries to choose optimal gamma value)
%  ROTOPT   : Rotation parameter for those having an optional parameter (Promax and Oblimin)
%  MAT_TYPE	: Matrix type (SCP: sums-of-squares-cross-product matrix, COV: variance-covariance matrix, COR: correlation matrix)
%  NUM_FAC	: Number of factors retained
%  theData  : Structured array with the data and accompanying information.  See readData.
%  LOADING	: Loading normalization ('K' Kaiser normalization, 'N' no
%  normalization, 'C' covariance loadings, 'W' Cureton-Mulaik weighting).
%
%  OPTIONAL
%  GAVE     : Convert the data to grand average form, perform the analysis,
%               then reconvert the factor scores back to average form, 'Y' or 'N'.
%  badData  : sparse matrix with same dimensions as theData where 1 denotes bad data.
%  crossFacCof : factor scoring coefficients from previous PCA for cross-verification
%
%Outputs
%  FactorResults: Structured variable with results and information about analysis.
%  .FacPat	: Factor pattern matrix - produces standardized variables from scores, scaled by communality  (rows=variables, cols=factors)
%  .FacStr	: Factor structure matrix - correlation between factors and variables (rows=variables, cols=factors)
%  .FacScr	: Factor scores, variance standardized, non-mean corrected  (rows=variables, cols=factors)
%             The ordering of the data dimensions in the rows (which ones vary the most quickly) is in the order of the
%             seven data dimensions (earlier varying faster), with the one serving as the column dimension left out.
%  .FacCof   : Factor scoring coefficients (rows=variables, cols=factors)
%  .FacCor	: Factor correlations
%  .scree    : The full set of unrotated eigenvalues, prior to truncation.
%  .facVar   : %age of variance accounted for by each rotated factor, ignoring other factors.
%			  Calculated based on variance of relationship matrix (if COV or SCP, each variable weighted).
%  .facVarQ  : %age of variance uniquely accounted for by each rotated factor.
%			  Calculated based on variance of relationship matrix (if COV or SCP, each variable weighted)
%  .facVarTot   : Total variance accounted for.
%  .varSD    : Standard deviations of the variables.
%  .PCAmode     : PCAmode used
%  .ROTATION    : Rotation used
%  .MAT_TYPE    : Relationship matrix used
%  .KAISER      : Loading option used
%  .timepoints  : Number of time points in original data
%  .numchan     : Number of channels in original data
%  .numFreqs    : Number of frequency bins in the original data
%  .numSubs     : Number of subjects.
%  .numCells    : Number of cells
%  .numFacs     : Number of factors
%  .montage     : The electrode montage.
%  .chanNames   : The channel names.
%  .timeNames   : The timepoint names.
%  .freqNames   : The frequency names.
%  .subNames    : The subject names.
%  .cellNames   : The cell names.
%  .facNames  : The factor names (only for factor files)
%  .chanTypes : The type of the channel: EEG, MGA (axial MEG), MGP (planar MEG), ECG, ANS (autonomic), REG (regional average)
%  .subTypes  : The type of the subject: RAW (single trial), AVG (subject average), GAV (grand average)
%  .cellTypes : The type of the cell: SGL (one cell), CMB (combination of cells)
%  .facTypes  : The type of the factor: SGL (one factor), CMB (combination of factors)
%  .subjectSpecs     : Cell array (subject,spec) of specific information for subjects
%  .subjectSpecNames : Cell array of the name of each subject spec type.
%  .taskSpecs     : Array (subject,task,measure) of task performance for subjects.  Can be empty if no specs.
%  .taskNames : Cell array of the name of each task.  Can be empty if no specs.
%  .taskMeasNames : Cell array of the name of each task measure.  Can be empty if no specs.
%  .Fs          : The sampling frequency in Hz.
%  .baseline    : The number of samples in the baseline.
%  .ename       : The name of the experiment.
%  .dataType    : The type of the data: 'analysis'
%  .fileName  : Name of original file.
%  .history   : Command used to create current file.
%  .fileFormat: The file format.
%  .ced       : The name of the .ced file for electrode coordinates.
%  .eloc      : The electrode location information, one for each channel (see readlocs header)
%  .implicit  : The electrode information for implicit references and fiducial locations (see readlocs header)
%    .reference
%        .original    : Original recording reference(s): 1-2 numbers in array
%        .current     : Current reference channel(s): 1-2 numbers in array
%        .type        : Current reference type: REG (regular), AVG (average reference, even if no longer adds to zero), CSD (current source density)
%  .badObs   : sparse column vector of same number of rows as the factor score matrix where 0=good data and 1=bad data.
%    .noise     : 4D matrix [channels, time points, cells/trials, subjects] mirroring .data.
%                 This is the +/- reference average (Schimmel, 1967) which provides an estimate of the noise level
%                 in averaged data by flipping every other trial.  Not applicable to spectral data.  Not from combining data other than across subjects.
%                 For grand average data, every other average is flipped.
%
%  doPCA runs the PCA on the data.  It can do either spatial or temporal,
%  using a number of rotations and rotation options.
%  See accompanying readme file for information about appropriate references for the rotations.

FactorResults=[];
varSD = std(data);
if ~isreal(data) %if the data has an imaginary component, as in spectral data
    data=[real(data);imag(data)];
    isComplex=1;
else
    isComplex=0;
end
D = diag(varSD);  %diagonal matrix of standard deviations of variables

R=cov(data,'partialrows');
[~,R]=cov2corr(R);
%R=corrcoef(data);  % OLD VERSION - does not allow for NaNs

goodVars=find(std(data) ~= 0);

%covariance matrix
if nargin>3
    S = varargin{1};
else
    S = cov(data,'partialrows'); 
end
Sd = diag(sqrt(diag(S)));  %diagonal matrix of standard deviations of variables as used to generate relationship matrix

%disp('eigenvariate decomposition...')
if nargin>4
    % generalized eigenproblem - this is primary setup for CCA analysis
    Sb = varargin{2};
    Sd = diag(sqrt(diag(Sb)));  %diagonal matrix of standard deviations of variables as used to generate relationship matrix
    [V,L] = eigs(S,Sb,NUM_FAC);
%elseif size(data,2)>size(data,1)
%    [V,L] = eigs(S,NUM_FAC);
else 
    [V,L] = eig(S);
end

if ~all(isreal(V))
    %Rounding errors resulted in relationship matrix being non-symmetric, which in turn resulted in complex numbers.
    %To fix, lower half of relationship matrix will be flipped to upper half to form symmetric relationship matrix.
    disp('Making relationship matrix symmetric to address effects of rounding errors.');
    Sfix=tril(S)+tril(S)'-diag(diag(S));
    if nargin>3
        [V,L] = eigs(Sfix,Sb,NUM_FAC);
    %elseif size(data,2)>size(data,1)
    %    [V,L] = eigs(Sfix,NUM_FAC);
    else
        [V,L] = eig(Sfix);
    end
end

%since the output of eig is not always sorted in order of ascending size, first sort them
[B,IX] = sort(diag(L),'descend');
L=L(IX,IX);
V=V(:,IX);


V = V(:,1:NUM_FAC);  %truncated eigenvector matrix
scree = diag(L);
L = L(1:NUM_FAC,1:NUM_FAC);  %truncated eigenvalue matrix

%factor scores
FacScr = (data) * V;    %factor scores, not mean corrected.

ScrDiag = diag(nanstd(FacScr));
A = inv(Sd) * (V * ScrDiag);  %unrotated factor loading matrix
%Note, this equation is commonly cited in statistics books but is misleading in this form.  The full
%form is A = inv(Sd) * (inv(V') * sqrt (L))
%this means that A is not in fact a scaled form of V as is commonly implied.
%This makes sense since A is a factor loading matrix (converts scores to raw data)
%while V is a scoring coefficient matrix (converts raw data to scores).
%inv(V') does reduce down to V, though, since X'=inv(X)
%for an orthogonal matrix, X.

switch type.method

    case 'dien'
        failed =1;
        while failed==1
            C = sum((A.^2),2)';
            % Kaiser normalisaiton
            Ak = (diag(sqrt(C).^-1)) * A;  %factor loadings Kaiser-normalized by communalities

            % ROTATION
            [FacPat]= ep_doVarimax(Ak);
            FacCor = diag(ones(NUM_FAC,1));
            FacStr = FacPat;

            %renormalize factor loadings by original communalities    
            FacPat = diag(sqrt(C)) * FacPat;  
            FacStr = diag(sqrt(C)) * FacStr;

            if strcmp(type.rotation_option,'promax') && type.degree>1
                %Only apply loading weighting to the Varimax step
                [FacPat, FacCor] = ep_doPromax(FacPat, type.degree); 
                %to match SAS output and to avoid rounding errors.
                FacStr = FacPat * FacCor;	%factor structure matrix (Harman, eq. 12.19, p. 268)
            end

            if NUM_FAC == 1 %ensure that matrices have the right orientation when there is only one factor.
                if size(FacPat,1) < size(FacPat,2)
                    FacPat=FacPat';
                end
                if size(FacStr,1) < size(FacStr,2)
                    FacStr=FacStr';
                end
            end

            LargestLoading=max(max(abs(FacStr))); %factor pattern loadings can go over 1 for oblique rotations
            LargestCom=max(max(abs(sum(FacPat.*FacStr,2))));
%             if NUM_FAC > 1
% 
%                 %Deal with loadings that are too large
%                 if round(LargestLoading*100) > 120 %allow very small violation of factor loading limit due to rounding errors
%                     %disp('Loadings are over the permitted maximum of 1.  It appears this rotation has crashed.');
%                     NUM_FAC = NUM_FAC-1;
% 
%                 elseif round(LargestCom*100) > 120 %allow very small violation of communality limit due to rounding errors
%                     %disp('Communalities are over the permitted maximum of 1.  It appears this rotation has crashed.');
%                     NUM_FAC = NUM_FAC-1;
%                 else failed = 0;
%                 end
%             else
                failed = 0;
%             end

        end
    case 'matlab'
        if strcmp(type.rotation_option,'promax')
            [FacPat, FacCor] = rotatefactors(A,...
                'Method',type.rotation_option,... 
                'Power',type.degree,... % Promax only
                'Normalize','on',...
                'Reltol',.00001,...
                'Maxit',1000);
            FacStr = FacPat * FacCor;	%factor structure matrix (Harman, eq. 12.19, p. 268)
        else
            [FacPat, FacCor] = rotatefactors(A,...
                'Method',type.rotation_option,... % 'Method' is 'orthomax', 'varimax', 'quartimax', 'equamax', or 'parsimax'
                'Normalize','on',...
                'Reltol',.00001,...
                'Maxit',1000);
            FacStr = FacPat;
        end
        LargestLoading=max(max(abs(FacStr))); %factor pattern loadings can go over 1 for oblique rotations
        LargestCom=max(max(abs(sum(FacPat.*FacStr,2))));
end

invR=pinv(R);
FacCof=invR*FacStr;
FacScr=(data)*FacCof;

facVar=[];
facVarQ=[];
Comm=[];
    
if LargestCom ~=0 && LargestLoading~=0

    if nargin==2 % this is just to suppress a warning that comes up for CCA analysis
        FacScr=(FacScr)*inv(diag(std(FacScr))); %Standardize factor scores, not mean corrected.
    end
    
    Var=Sd.^2;

    %calculate communalities (equivalent to calculating R-Squared in multiple regression).

    %Cohen, J. & Cohen, P. (1983). Applied multiple regression/correlation
    %analysis for the behavioral sciences. Hillsdale, NJ: Lawrence Erlbaum Associates.

    Comm=sum((Var*FacPat.*FacStr),2)/sum(diag(Var)); % p. 100

%     if strcmp(theRotationName,'none')
%         theRotationName2='unrotated';
%     else
%         theRotationName2=theRotationName;
%     end;
%     
%     disp(['Amount of variance accounted for by the ' theRotationName2 ' solution is ' sprintf('%3.2f%%.',100*sum(Comm,1))]);

    facVar = sum((Var*FacPat.*FacStr))/sum(diag(Var));

    %calculate unique variance of factors by calculating semi-partial correlation
    %between the factor and the variable.  p. 470 of C&C (1983)
    
    facVarQ=sum(Var*(FacPat*inv(diag(sqrt(diag(inv(FacCor)))))).^2)/sum(diag(Var));
    
    %sort factors in order of size
    [dummy index]=sort(facVar); 
    index = fliplr(index);

    FacPat=FacPat(:,index);
    FacStr=FacStr(:,index);
    FacCof=FacCof(:,index);
    FacScr=FacScr(:,index);
    FacCor=FacCor(index,index);
    facVar=facVar(index);
    facVarQ=facVarQ(index);

    for f1 = 1:NUM_FAC
        if sum(FacPat(:,f1)) < 0	%flip factor loadings so mostly positive
            FacPat(:,f1) = FacPat(:,f1) .* (-1);
            FacStr(:,f1) = FacStr(:,f1) .* (-1);
            FacCof(:,f1) = FacCof(:,f1) .* (-1);
            FacScr(:,f1) = FacScr(:,f1) .* (-1); %if the loading is flipped, the scores must be flipped too.
            FacCor(:,f1) = FacCor(:,f1) .* (-1);
            FacCor(f1,:) = FacCor(f1,:) .* (-1);
        end
    end
end

facNames=[];
facTypes=[];

factorDigits=length(num2str(NUM_FAC));
for i=1:NUM_FAC
    facNames{i,1}=['TF' sprintf(['%0' num2str(factorDigits) 'd'],i)];
    facTypes{i,1}='SGL';
end

% FactorResults.FacPat=zeros(size(badData,2),size(FacPat,2));
% FactorResults.FacScr=zeros(size(badData,1),size(FacScr,2));
FactorResults.FacPat(goodVars,:)=FacPat;
% FactorResults.FacStr=zeros(size(badData,2),size(FacStr,2));
FactorResults.FacStr(goodVars,:)=FacStr;
FactorResults.FacScr=zeros(size(data,1),NUM_FAC);
if isComplex
    FactorResults.FacScr=complex(FacScr(1:size(FacScr,1)/2,:),FacScr(size(FacScr,1)/2+1:end,:));
else 
    FactorResults.FacScr=FacScr;
end
FactorResults.FacCor=FacCor;
FactorResults.FacCof(goodVars,:)=FacCof;
FactorResults.scree(1:length(scree))=scree;
FactorResults.facVar=facVar;
FactorResults.facVarQ=facVarQ;
FactorResults.facVarTot=sum(Comm,1);
FactorResults.varSD(goodVars)=varSD;
FactorResults.numFacs=NUM_FAC;
FactorResults.facNames=facNames;
FactorResults.facTypes=facTypes;

function [Ak]= ep_doVarimax(A)

%[Ak]= ep_doVarimax(A);  - Compute varimax rotation for PCA
%  Gives essentially identical results to SAS command: PROC FACTOR rot = varimax
%
%Inputs
%  A		: Unrotated factor loadings matrix (variables, factors)
%
%Outputs
%  Ak		: Rotated factor pattern matrix
%
%History
%
%  Cureton, E. E. & Mulaik, S. A. (1975).  The weighted varimax rotation and the promax rotation.
%  Psychometrika, 40(2) 183-195.
%
%  Kaiser, H. F. (1959).  Computer program for varimax rotation in factor analysis.
%  Educational and Psychological Measurement, 19(3) 413-420.
%
%  Harman, H. H.  (1976).  Modern factor analysis, 3rd edition.  Chicago:University of Chicago Press.
%
%  by Joseph Dien (4/99)
%  jdien07@mac.com
%
%  modified (9/30/00) JD
%  Modified to allow Kaiser normalization to be turned off.
%
%  bugfix (2/27/01) JD
%  Fixed bug in Kaiser correction code (not undoing it when turned on)
%
%  modified (3/1/01) JD
%  Added direct rotation of factor scores.
%
%  bugfix (3/24/01) JD
%  Factor scores reordered and flipped when factor loadings are.
%
%  modified (5/27/01) JD
%  Factor scoring coefficients added.
%
%  modified (8/14/01) JD
%  Deleted rotation of V to minimize effect of cumulative rounding errors
%
%  modified (7/26/02) JD
%  Added non-convergence warning.
%
%  modified (10/22/02) JD
%  Lowered minimum criteria for rotation to .00001 since Kaiser suggested
%  .00116 was too loose and was sometimes not rotating when it should.
%  Added warning if no rotation occurred.
%
%  modified (2/15/04) JD
%  Added accomodation for covariance loading option
%
%  modified (01/09/05) JD and Dan Beal
%  Added Weighted Varimax option (Cureton and Mulaik, 1975)
%
%  modified (2/4/08) JD
%  Removed manual rotation of factor scores as it turned out to be less
%  accurate.
%
%  modified (2/12/08) JD
%  Moved factor loading normalization and factor loading ordering
%  to doPCA so that they can be applied to all rotations.
%
%  modified (2/27/08) JD
%  Random starting rotation procedure added to avoid local minima.
%
%  bugfix (3/24/08) JD
%  Fixed orientation problem with Ak when only one factor.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Copyright (C) 1999-2020  Joseph Dien
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NumReps=10;
NumFacs=size(A,2);
NumVars=size(A,1);
FacPatReps=zeros(NumReps,NumVars,NumFacs);
ThReps=zeros(NumReps,NumFacs,NumFacs);
PhiReps=zeros(NumReps,NumFacs,NumFacs);
fReps=zeros(NumReps,1);
rotTot=zeros(NumReps,1);
didRot=zeros(NumReps,1);
for repetition =1:NumReps
    %disp(['erp_pca: repetition ' num2str(repetition) '...'])

    MaxRotations = 1000;	%Maximum number of rotations to do before calling it a day
    C = sum((A.^2)')';      %vector of communalities
    NUM_FAC = size(A,2);	%Number of factors retained
    NUM_VAR = size(A,1);	%Number of variables
    Ak=A;
    
    [T,R]=qr(randn(NUM_FAC,NUM_FAC)); %random initial rotation
    Ak=A*inv(T)';

    %  Start rotation

    counter = NUM_FAC * (NUM_FAC -1)/2;  %initialize counter
    nrot = 0;	%initialize rotation counter
    rotstats = [];
    rotated = 0;

    while (counter > 0) && (nrot < MaxRotations)  %do rotations until has gone through all w/o rotating or exceeded max rotations.
        for f1 = 1:(NUM_FAC-1)
            for f2 = (f1 + 1):NUM_FAC
                A1 = Ak(:,f1);
                A2 = Ak(:,f2);
                Ar = [A1 A2];
                ru = A1.^2 - A2.^2;
                rv = 2 * A1 .* A2;
                rA = sum(ru);
                rB = sum(rv);
                rC = sum(ru.^2 - rv.^2);
                rD = 2*sum(ru .* rv);

                NUM = rD - (2 * rA * rB)/NUM_VAR;
                DEN = rC - (rA.^2 - rB.^2)/NUM_VAR;
                if abs(NUM/DEN) > .00001 %(would rotation be enough to bother with?)
                    rotated = 1;

                    G = sqrt(NUM^2+DEN^2);	%Following computational variant due Wingersky (Harman, p. 294)
                    cos4phi = DEN/G;
                    cos2phi = sqrt((1+cos4phi)/2);
                    cosphi = sqrt((1+cos2phi)/2);
                    sinphi = sqrt((1-cos2phi)/2);

                    if NUM < 0
                        sinphi = -abs(sinphi);
                    else
                        sinphi = abs(sinphi);
                    end;

                    rot = [cosphi -sinphi; sinphi cosphi];
                    Ar = Ar * rot;
                    Ak(:,f1) = Ar(:,1);
                    Ak(:,f2) = Ar(:,2);
                    counter = NUM_FAC * (NUM_FAC -1)/2;

                else
                    counter = counter -1; %no rotation occurred so reduce counter.
                end;
            end;
        end;
        nrot = nrot +1;	%Increment number of rotations.  It is expected that after a certain point just not worth it.
    end;

    FacPatReps(repetition,:,:)=Ak;
    fReps(repetition,1)=sum(var(Ak.^2),2); %The Varimax criterion
    rotTot(repetition)=nrot;
    didRot(repetition)=rotated;
end;

[B I]=sort(fReps); %Find the global maximum of the computed rotations.

Ak=squeeze(FacPatReps(I(NumReps),:,:));

if size(Ak,1)==1 %FacPat oriented wrong way when just one factor
    Ak=Ak';
end;

if (sum(didRot) == 0)
    if NumFacs == 1
        %disp('Since there was only one factor, no rotation occurred.');
    else
        %disp(['Warning - no rotation occurred.']);
    end;
end;

if (min(rotTot) >= MaxRotations)
    disp(['Warning - solution did not converge, meaning that it was not able to reach a stable solution within a reasonable number of rotations.']);
end;

function [PmxPat, phi] = ep_doPromax(VmxPat, k);

%[PmxPat, phi] = ep_doPromax(VmxPat, k)  - Compute promax rotation for PCA
%  Gives similar results to SAS command: PROC FACTOR rot = promax
%
%Inputs
%  VmxPat	: Varimax rotated factor loading matrix (variables, factors)
%  k		: Power to raise loadings to produce target matrix.  Higher power results in more oblique solutions.
%				4 better for simple factor structures, 2 for more complex, 3 is a good compromise (SAS default value)
%
%Outputs
%  PmxPat	: Factor pattern matrix
%  phi  : Correlations between factors
%
%History
%
%  With assistance from Lew Goldberg and Jack Digman.
%
%  First proposed by:
%  Hendrickson, A. E. & White, P. O. (1964).  Promax: A quick method for
%  rotation to oblique simple structure.  The British Journal of Statistical
%  Psychology, 17:65-70.
%
%  This algorithm uses the original procedure of rotating the reference
%  vectors and is therefore an "indirect Promax" rotation.
%
%  Normalization process attributed to:
%  Cureton, E. E., & D'Agostino, R. B. (1983). Factor Analysis: An applied approach. Hillsdale, NJ: Lawrence Erlbaum and Associates.
%
%  Harman, H. H.  (1976).  Modern factor analysis, 3rd edition.  Chicago:University of Chicago Press.
%
%  Gorsuch, R. L.  (1983).  Factor analysis, 2nd edition.  Hillsdale, NJ:Lawrence Erlbaum Associates.
%
%  Dillon, W. R. & Goldstein, M. (1984).  Multivariate analysis: Methods and applications.  New York:Wiley & Sons.
%
%  by Joseph Dien (4/99)
%  jdien07@mac.com
%
%  12/7/00 JD
%  Fixed error in promax algorithm.  Modified to output factor correlations and reference structure.
%  Given the same varimax solution, now produces identical results to SAS 6 promax output.
%
%  modified (4/3/01) JD
%  Added manual rotation of factor scores
%
%  bugfix 1/12/03 JD
%  Fixed column-normalization of H to be on absolute value.
%
%  modified (2/4/08) JD
%  Removed manual rotation of factor scores as it turned out to be less
%  accurate.

%     Copyright (C) 1999-2008  Joseph Dien
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

NUM_FAC = size(VmxPat,2);	%Number of factors retained
NUM_VAR = size(VmxPat,1);	%Number of variables

H = diag(1./sqrt(sum(VmxPat'.^2)))*VmxPat;	%Normalize rows - equalize sum of squares of each variable's loadings
H = H * diag(1./max(abs(H)));			%Column-normalize by highest absolute value in each column - tends to equalize factor sizes

H = H.^k;		%compute target matrix by taking higher power of factor loadings

H = abs(H) .* sign(VmxPat);		%add the signs back to the target matrix since even powers eliminate them

lambda = inv(VmxPat'*VmxPat) * VmxPat' * H;	%nonnormalized transformation matrix between starting factor matrix and reference vector (H & W, p. 66)

lambda = lambda*diag(1./sqrt(sum(lambda .^2)));		%Normalize the columns of lambda (H & W, p. 66)

psi = lambda' * lambda;		%Correlations of reference vectors (Harman, eq. 12.29, p. 273)                               

r = inv(psi);					%Inverse of psi

D = diag(1./sqrt(diag(r)));	%relationship between reference and primary axes (Gorsuch, p. 226)

T = lambda*inv(D);  %Procrustean Transformation Matrix (Gorsuch, eq. 10.2.9, p. 226)

phi = inv(T)*inv(T)';	%Correlations between oblique primary factors (Gorsuch, eq. 10.1.15, p. 215)

PmxPat = VmxPat * T; 		%factor pattern matrix (Gorsuch, eq. 10.2.2, p. 220)


