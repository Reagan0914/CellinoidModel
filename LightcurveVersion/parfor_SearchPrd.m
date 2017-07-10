%%  @parfor_SearchPrd.m
%
%   Descriptions:
%       Apply the parallel pool to search Period
%   Edited by LUXP
%   Date: 2016-10-06
function parfor_SearchPrd

if exist('mexEigFunction.m','file') == 0
    addpath('./Others');
end
if exist('OUTPUT','dir') == 0
        system('mkdir OUTPUT');
end
%%  NOTICE: setting the searching interval in SearchPrdInterval function; 
%MaxPools = 4;
parfor i=1:8
    SearchPrdInterval(i);
end

end

function SearchPrdInterval(i)
%%  Precalculate 
global JDT  E0 E  BRT LCm LebData LebInd
[LC, LCm] = fReadData;
JDT = LC(:,1); BRT= LC(:,2); E0 = LC(:, 3:5)'; E = LC(:, 6:8)';  

%%  Get Lebedev Normal Vectors
N = 590;
[LebData, LebInd] = GetLebData(N);
UpdateIndex = ones(10,1); UpdateIndex(1) = 0; 
%fixed Pole
UpdateIndex(7) = 0;UpdateIndex(8) = 0;
UpdateIndex(2:6)=[0,0,0,0,0]';

%%  Parallex Computing
DeltaPrd = 0.001;
switch i
    case 1
        PrdInterval = 1.5:DeltaPrd:2.5;
        [cg_History, Chisq_History] = SearchPrd(UpdateIndex, PrdInterval,1);
    case 2
        PrdInterval = 2.5:DeltaPrd:3.5;
        [cg_History, Chisq_History] = SearchPrd(UpdateIndex, PrdInterval,1);
    case 3
        PrdInterval = 3.5:DeltaPrd:4.5;
        [cg_History, Chisq_History] = SearchPrd(UpdateIndex, PrdInterval,1);
    case 4
        PrdInterval = 4.5:DeltaPrd:5.5;
        [cg_History, Chisq_History] = SearchPrd(UpdateIndex, PrdInterval,1);
    case 5
        PrdInterval = 5.5:DeltaPrd:6.5;
        [cg_History, Chisq_History] = SearchPrd(UpdateIndex, PrdInterval,1);
    case 6
        PrdInterval = 6.5:DeltaPrd:7.5;
        [cg_History, Chisq_History] = SearchPrd(UpdateIndex, PrdInterval,1);
    case 7
        PrdInterval = 7.5:DeltaPrd:8.5;
        [cg_History, Chisq_History] = SearchPrd(UpdateIndex, PrdInterval,1);
    case 8
        PrdInterval = 8.5:DeltaPrd:9.5;
        [cg_History, Chisq_History] = SearchPrd(UpdateIndex, PrdInterval,1);
end

end

%%  Search period with various initial test prd
function [cg_History, Chisq_History] = SearchPrd(UpdateIndex, PrdInterval, savedata)
%  cg Test Scheme;
    cg_test = zeros(10,10);
    cg_test(:,1) = cgModify([1, 0.7, 0.85, 0.8,  0.5,0.3,  5, 30, 5.5, 10]', 'M');
    cg_test(:,2) = cgModify([1, 0.8, 0.85, 0.6,  0.7,0.6,  50, -30, 5.5, 40]', 'M');
    cg_test(:,3) = cgModify([1, 0.9, 0.75, 0.6,  0.65,0.5,  85, 60, 5.5, 70]', 'M');
    cg_test(:,4) = cgModify([1, 0.6, 0.85, 0.6,  0.7,0.5,  135, 80, 5.5, 120]', 'M');
    cg_test(:,5) = cgModify([1, 0.7, 0.8, 0.7,  0.7,0.6,  170, -30, 5.5, 160]', 'M');
    cg_test(:,6) = cgModify([1, 0.8, 0.9, 0.7,  0.8,0.6,  225, 30, 5.5, 200]', 'M');
    cg_test(:,7) = cgModify([1, 0.9, 0.7, 0.5,  0.5,0.4,  260, 50, 5.5, 250]', 'M');
    cg_test(:,8) = cgModify([1, 0.8, 0.85, 0.5,  0.5,0.3,  315, -50, 5.5, 300]', 'M');
    cg_test(:,9) = cgModify([1, 0.7, 0.95, 0.6,  0.7,0.45,  345, 75, 5.5, 330]', 'M');
    cg_test(:,10) = cgModify([1, 0.5, 0.95, 0.5,  0.5,0.43,  25, -75, 5.5, 110]', 'M');
    
TotalTest = size(PrdInterval,2);
cg_History = zeros(10, TotalTest);
Chisq_History = zeros(TotalTest,1);
Ind = 1;
for i = 1:TotalTest
    cg0 = cg_test(:, Ind); Ind = mod(Ind, 10) + 1; 
    cg0(9) = 2*pi * 24 / PrdInterval(i);
    %%Fix to Ellipsoid and  fix pole to z-axis%%%%
    cg0(1:6) = [1, 1, 0.8,0.8, 0.6,0.6]';
    cg0(7) = 0; cg0(8) = 0; %in Math Model
    %%%%%%%%%%%%%%%%%%%%
    %cg0 = cgModify(cg0, 'M');
    [cg_fit, Chisq_History(i)] = LM_Inverse(cg0, UpdateIndex);
    cg_History(:,i) = cgModify(cg_fit, 'G');
    fprintf(1,'Finishing the search for %d in %d times, Chisq:=%f. Interval:[%2.1f, %2.2f]\n',i,TotalTest,Chisq_History(i), PrdInterval(1), PrdInterval(end));
end
%%  save data
if nargin == 3 && savedata == 1
    DateMarker = datestr(now);
    FileName = ['OUTPUT/Prd_', DateMarker([1:12,13,14,16,17]), 'From_',num2str(PrdInterval(1)),'_To_',num2str(PrdInterval(end)),'.mat'];
    if exist(FileName, 'file') ~= 0
        AddTail =num2str(rand);
        FileName = [FileName(1:end-4), '_', AddTail(3:end),'_temp.mat'];
    end
    save(FileName, 'cg_History', 'Chisq_History');
end
end
