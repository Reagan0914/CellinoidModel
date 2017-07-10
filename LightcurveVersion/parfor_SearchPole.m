%%  @parfor_SearchPole.m
%
%   Descriptions:
%       Apply the parallel pool to Search Pole
%   Edited by LUXP
%   Date: 2016-10-07
function parfor_SearchPole

if exist('mexEigFunction.m','file') == 0
    addpath('./Others');
end
if exist('OUTPUT','dir') == 0
        system('mkdir OUTPUT');
end
%%  NOTICE: setting the searching interval in SearchPrdInterval function; 
%MaxPools = 4;
parfor i=1:4
    SearchPoleInterval(i);
end

end

function SearchPoleInterval(i)
%%  Precalculate 
global JDT  E0 E  BRT LCm LebData LebInd
[LC, LCm] = fReadData;
JDT = LC(:,1); BRT= LC(:,2); E0 = LC(:, 3:5)'; E = LC(:, 6:8)';  

%%  Get Lebedev Normal Vectors
N = 590;
[LebData, LebInd] = GetLebData(N);
UpdateIndex = ones(10,1); UpdateIndex(1) = 0; 

%%  Parallex Computing
DeltaPole = 2;
%   Lam and beta interval in Math models: i.e. radial and beta from 0 to 180
switch i
    case 1
        LamInterval = 0:DeltaPole:90;
        BetaInterval = 0:DeltaPole:180;
        LamInterval = LamInterval*pi/180;  BetaInterval = BetaInterval*pi/180;
        [cg_History, Chisq_History] = SearchPole(UpdateIndex, LamInterval, BetaInterval,1);
    case 2
        LamInterval = 90:DeltaPole:180;
        BetaInterval = 0:DeltaPole:180;
        LamInterval = LamInterval*pi/180;  BetaInterval = BetaInterval*pi/180;
        [cg_History, Chisq_History] = SearchPole(UpdateIndex, LamInterval, BetaInterval,1);
    case 3
        LamInterval = 180:DeltaPole:270;
        BetaInterval = 0:DeltaPole:180;
        LamInterval = LamInterval*pi/180;  BetaInterval = BetaInterval*pi/180;
        [cg_History, Chisq_History] = SearchPole(UpdateIndex, LamInterval, BetaInterval,1);
    case 4
        LamInterval = 270:DeltaPole:360;  
        BetaInterval = 0:DeltaPole:180;
        LamInterval = LamInterval*pi/180;  BetaInterval = BetaInterval*pi/180;
        [cg_History, Chisq_History] = SearchPole(UpdateIndex, LamInterval, BetaInterval,1);
end

end

%%  Search period with various initial test prd
function [cg_History, Chisq_History] = SearchPole(UpdateIndex, LamInterval, BetaInterval, savedata)
%% Derived Best-fit period and phi0 from Prd_search, in Math models
BestFitPrd = 2*pi * 24 / 5.2710; 
BestFitPhi0=57.6453 * pi/180;

%%  cg Test Scheme for semi-axes;
    cg_test = zeros(10,10);
    cg_test(:,1) = [1, 0.7, 0.85, 0.8,  0.5,0.3,  5, 30, 5.5, 10]';
    cg_test(:,2) = [1, 0.8, 0.85, 0.6,  0.7,0.6,  50, -30, 5.5, 40]';
    cg_test(:,3) = [1, 0.9, 0.75, 0.6,  0.65,0.5,  85, 60, 5.5, 70]';
    cg_test(:,4) = [1, 0.6, 0.85, 0.6,  0.7,0.5,  135, 80, 5.5, 120]';
    cg_test(:,5) = [1, 0.7, 0.8, 0.7,  0.7,0.6,  170, -30, 5.5, 160]';
    cg_test(:,6) = [1, 0.8, 0.9, 0.7,  0.8,0.6,  225, 30, 5.5, 200]';
    cg_test(:,7) = [1, 0.9, 0.7, 0.5,  0.5,0.4,  260, 50, 5.5, 250]';
    cg_test(:,8) = [1, 0.8, 0.85, 0.5,  0.5,0.3,  315, -50, 5.5, 300]';
    cg_test(:,9) = [1, 0.7, 0.95, 0.6,  0.7,0.45,  345, 75, 5.5, 330]';
    cg_test(:,10) =[1, 0.5, 0.95, 0.5,  0.5,0.43,  25, -75, 5.5, 110]';
    
TotalLam = size(LamInterval,2);
TotalBeta = size(BetaInterval,2);
TotalTest = TotalLam*TotalBeta;
cg_History = zeros(10, TotalTest);
Chisq_History = zeros(TotalTest,1);
Ind = 1;
for lam = 1:TotalLam
    for beta = 1:TotalBeta
        cg0 = cg_test(:, Ind); Ind = mod(Ind, 10) + 1; 
        cg0(9) = BestFitPrd; cg0(10) = BestFitPhi0;
        cg0(7) = LamInterval(lam);      %in math model: NOTICE!
        cg0(8) = BetaInterval(beta);
        Index_record = (lam-1)*TotalBeta + beta;
        [cg_fit, Chisq_History(Index_record)] = LM_Inverse(cg0, UpdateIndex);
        cg_History(:,Index_record) = cgModify(cg_fit, 'G');
        fprintf(1,'Finishing the search for %d in %d times, Chisq:=%f. Interval:[%2.1f, %2.2f]\n',Index_record,TotalTest,Chisq_History(Index_record), LamInterval(1), LamInterval(end));
    end
end
%%  save data
if nargin == 4 && savedata == 1
    DateMarker = datestr(now);
    FileName = ['OUTPUT/Pole_', DateMarker([1:12,13,14,16,17]), 'From_',num2str(LamInterval(1)*180/pi),'_To_',num2str(LamInterval(end)*180/pi),'.mat'];
    if exist(FileName, 'file') ~= 0
        AddTail =num2str(rand);
        FileName = [FileName(1:end-4), '_', AddTail(3:end),'_temp.mat'];
    end
    save(FileName, 'cg_History', 'Chisq_History');
end
end
