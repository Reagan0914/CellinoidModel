%% @fReadData.m
%  Description:
%       Readin Observed Light Curve data, modify it for specific case
%
%  How To Use:
%   >> [LC, LCm] = fReadData;
%   
%   Edited by LUXP
%   Date: 2016-09-26
function [LC, LCm] = fReadData
LC = load('JPL_STD_Kaasalainen_433.txt');
LCm = load('JPL_STD_Kaasalainen_LCinfo_433.txt');

%%  for Test the first 10 lightcurves
%LC=LC(1:sum(LCm(1:10)),:);
%LCm=LCm(1:10);

%%  For Generate Test Light Curves
% Run: > SynLC = CalSynLC([9,1,4,8,1.5,2.5,0,0,5,0,-1,1,0.1])
% a1=9;a2=1;b1=4;b2=8;c1=1.5;c2=2.5;
% a1=2;a2=1;b1=0.9;b2=0.9;c1=0.8;c2=0.8;      %ellipsoid
% a1=4;a2=4;b1=4;b2=4;c1=4;c2=4; %sphere
% E = [3,1,0];
% E0 = [1,4,7];
% E=E/norm(E); 
% E0=E0/norm(E0); 
% Prd = 5;
% JDT = 0:0.01:2*Prd;
% JDT = JDT'/24;
% LCm = size(JDT,1);
% LC=[JDT, JDT, repmat(E0,LCm,1),repmat(E, LCm,1)];

%%  Statistcial The Solar Phase angle
Stat_SolarPhase = 1;
if Stat_SolarPhase == 1
    SolarPhase = zeros(size(LC,1),1);
    SolarPhaseStat = [];    %max, min, mean for each lightcurves
    for i = 1:size(LCm,1)
            EndInd = sum(LCm(1:i));
            Ind_LC = (EndInd - LCm(i) + 1): EndInd;
        for j = 1:LCm(i)
            Ind = EndInd - LCm(i) + j;
            SolarPhase(Ind) = acos(dot(LC(Ind, 3:5), LC(Ind,6:8))) * 180/pi;
        end
            SolarPhaseStat = [SolarPhaseStat; max(SolarPhase(Ind_LC)), min(SolarPhase(Ind_LC)), mean(SolarPhase(Ind_LC))];
    end
    fprintf(1, 'max deviation (in degree) of solar phase in Lightcurves: %f \n', max(abs(SolarPhaseStat(:,1) - SolarPhaseStat(:,2))));
    %%  plot the distribution of solar phase angle
    SA = SolarPhaseStat(:,3);   SA = sort(SA);
    figure, plot(SA,'*');
    set(gcf,'Position',[100,200,725,770]);
    xlabel('Lightcurves','fontsize',22);
    set(gca,'FontName','Times New Roman','FontSize',16,'fontweight','bold');
    ylabel('Solar Phase Angle (Degree)','fontsize',22);
end

end