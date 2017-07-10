%% @SearchProcessRefine
%       Search the possible best-fit after the pole search and period
%       search
%
%   Edited by LUXP
%   Date: 2016-10-09
function SearchProcessRefine
global JDT  E0 E  BRT LCm LebData LebInd
[LC, LCm] = fReadData;
JDT = LC(:,1); BRT= LC(:,2); E0 = LC(:, 3:5)'; E = LC(:, 6:8)';  

%%  Get Lebedev Normal Vectors
N = 590%1202;%590;
[LebData, LebInd] = GetLebData(N);

if exist('mexEigFunction.m','file') == 0
    addpath('./Others');
end

%% derived bestfit
cg0 = [1.0000
    1.0055
    0.1995
    0.0280
    0.5459
    0.5428
  197.6247
   12.7097
    5.2710
   44.9749];
cg0 = [1.0000
    1.0055
    0.5459
    0.5428
     0.1995
    0.0280
  17.6247
   12.7097
    5.2710
   44.9749];
cg0 = cgModify(cg0, 'M');
UpdateIndex = ones(10,1); UpdateIndex(1) = 0;

[cg_fit, Chisq] = LM_Inverse(cg0, UpdateIndex);
cg = cgModify(cg_fit, 'G');
cg=1
end

