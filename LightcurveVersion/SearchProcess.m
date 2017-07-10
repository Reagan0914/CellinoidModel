%% @SearchProcess
%
%
%   Edited by LUXP
%   Date: 2016-10-06
function SearchProcess
global JDT  E0 E  BRT LCm LebData LebInd
[LC, LCm] = fReadData;
JDT = LC(:,1); BRT= LC(:,2); E0 = LC(:, 3:5)'; E = LC(:, 6:8)';  

%%  Get Lebedev Normal Vectors
N = 590;
[LebData, LebInd] = GetLebData(N);

if exist('mexEigFunction.m','file') == 0
    addpath('./Others');
end
 
cg0 = [1, 0.9, 0.75, 0.75, 0.5, 0.5, 100, 50, 5, 120]';cg0 = cgModify(cg0, 'M');
UpdateIndex = ones(10,1); UpdateIndex(1) = 0;

PrdInterval = 5.1:1:5.2;
[cg_History, Chisq_History] = SearchPrd(UpdateIndex, PrdInterval,1);

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
    cg0(9) = PrdInterval(i);
    cg0 = cgModify(cg0, 'M');
    [cg_fit, Chisq_History(i)] = LM_Inverse(cg0, UpdateIndex);
    cg_History(:,i) = cgModify(cg_fit, 'G');
    fprintf(1,'Finishing the search for %d in %d times, Chisq:=%f.\n',i,TotalTest,Chisq_History(i));
end
%%  save data
if nargin == 3 && savedata == 1
    if exist('OUTPUT','dir') == 0
        system('mkdir OUTPUT');
    end
    DateMarker = datestr(now);
    FileName = ['OUTPUT/Prd_', DateMarker([1:12,13,14,16,17]), '.mat'];
    save(FileName, 'cg_History', 'Chisq_History');
end
end