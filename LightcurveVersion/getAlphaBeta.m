%%  @getAlphaBeta
%           invoked by LM_Inverse.m, to provide efficiently the Hessian matrix \Alpha and derivative vector \Beta in L-M
%           algorithm, as well as computing the chisquare
%       NOTICE: UpdateIndex could indicate which parameter will be optimized! 1-optimize, 0-not        
%   
%   Edited by LUXP
%   Date: 2016-10-06
function [A,B,Chisq] = getAlphaBeta(cg, UpdateIndex)

global BRT LCm

%% pre-allocate the matrix and initial settings
TotalNumData = size(BRT,1);
TotalNumUpdateParas = sum(UpdateIndex);
TotalNumParas = size(cg,1);

A = zeros(TotalNumUpdateParas,TotalNumUpdateParas);
B = zeros(TotalNumUpdateParas, 1);

Delta = 0.000001;
Chisq = 0;
TotalLC = size(LCm,1); % lightcurves number
%Compute the synthetic brightnesses and the derivative W.R.T N paras
SynBrt = LM_CalSynLC(cg);  
%   Calculate dyda by diffrentiation tool
dydaLC = zeros(TotalNumData, TotalNumUpdateParas);
Ind = 1;
for j=1:TotalNumParas
    if UpdateIndex(j) == 1
        cg_da = cg; cg_da(j) = cg(j) + Delta;
        Temp_SynLC = LM_CalSynLC(cg_da);
        dydaLC(:, Ind)=(Temp_SynLC - SynBrt)/Delta;
        Ind = Ind +1;
    end
end 

%%  Calculate \Alpha and \Beta in L-M algorithm
for LC_i=1:TotalLC    
    iEnd=sum(LCm(1:LC_i));  %start to end number in Total LC
    iStart=iEnd-LCm(LC_i) + 1;
    LC_Index = iStart:iEnd;
    
    LC_Brt = SynBrt(LC_Index);
    LC_dyda = dydaLC(LC_Index,:);
    LC_M = LCm(LC_i);
    
    sumBRT=sum(LC_Brt);
    dy = BRT(LC_Index) - LC_Brt*LC_M/sumBRT;
    Chisq = Chisq + dy'*dy;
    %% calculate dydavec for one lightcurves
    dydaVec = zeros(TotalNumUpdateParas,1);
    LC_dydaVec = zeros(TotalNumUpdateParas, LC_M);
    for k = 1:TotalNumUpdateParas
        LC_dydaVec(k, :) = (LC_M/sumBRT) * (LC_dyda(:,k) - LC_Brt*sum(LC_dyda(:,k))/sumBRT);
        dydaVec(k) = sum(dot(dy, LC_dydaVec(k,:)));
    end
    B = B + dydaVec;
    
    for j = 1:LC_M
        Temp = LC_dydaVec(:, j);
        A  = A + Temp*Temp';
    end
end

end