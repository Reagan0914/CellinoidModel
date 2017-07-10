%%  @LC_CMP_Syn_Obs.m
%
%   Descriptions:
%       compare the fitted synthetic lightcurves and Observed lightcurves
%   
%   Edited by LUXP
%   Date: 2016-10-08
function LC_CMP_Syn_Obs
%% bestfit cg from period and pole search
BestFitCg = [1.0000
    1.0055
    0.1995
    0.0280
    0.5459
    0.5428
  197.6247
   12.7097
    5.2710
   57.6453];
BestFitCg = cgModify(BestFitCg,'M');

%%  observed data
global JDT  E0 E  BRT LCm LebData LebInd
[LC, LCm] = fReadData;
JDT = LC(:,1); BRT= LC(:,2); E0 = LC(:, 3:5)'; E = LC(:, 6:8)'; 
addpath('./Others');
%%  Get Lebedev Normal Vectors
N = 590;
[LebData, LebInd] = GetLebData(N);

SynLC = LM_CalSynLC(BestFitCg);

%% Show the plots
for LCi=1:size(LCm,1)
    IndEnd = sum(LCm(1:LCi));
    Ind = (IndEnd-LCm(LCi)+1):IndEnd;
    jdh = (JDT(Ind)-JDT(Ind(1)))*24;%julian date -> hours
    obsbrt = BRT(Ind);
    synbrt = SynLC(Ind);
    figure(11),
    mean(obsbrt)
    plot(jdh, obsbrt,'*');
    hold on;
    plot(jdh, synbrt/mean(synbrt),'-r');
    pause,
    close all;
    
end
end