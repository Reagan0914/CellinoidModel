%% @CalSynLC
%  Description:
%       Calculate the synthetic lightcurve By applying Lebedev Quadrature
%         cg:= [six axes, lam,beta of pole, period and intial Phase, three scattering factors] 
%                   Total 10+3 = 13 parameters
%         JDT = observation juliandate vectors in one lightcurves
%
%  How To Use:
%   >> SynLC = CalSynLC(cg, JDT);
%   
%   Edited by LUXP
%   Date: 2016-09-26
function SynLC = CalSynLC(cg)
addpath('./Others');
if nargin == 0 
    cg  = [1, 0.9, 0.85, 0.75, 0.7, 0.5, 100, 50, 5, 120, -0.4, 1, 0.15]'; 
end
cg = cgModify(cg, 'M');

%%  Get Lebedev Normal Vectors
if exist('LebData','var') ~= 1
    N = 590;
    [LebData, LebInd] = GetLebData(N);
end
Ndegree = size(LebData,1);

%%  Calculate Normal vector and facet area for Lebedev Method
[NormalVec, Gw] = CalNormalVector(LebData, LebInd, cg);

%%  Get unit directions of Sun(E0) and Earth(E) at asteroid-centric system
[JDT_E0, JDT_E, JDT, LCm] = GetSunEarth(cg);

%%  Generate the Synthetic lightcurves W.R.T. JDT and E0, E and cg
SynLC = zeros(size(JDT));
for tt = 1:size(JDT,1) 
    %E0 = JDT_E0(:,tt);E = JDT_E(:,tt);
    SynLC(tt) = CalBrtPoint(Ndegree, JDT_E(:,tt), JDT_E0(:,tt), NormalVec, Gw);
end

%%  Figure out Synthetic Light curve
%   JDT: Julian Date, 
%ShowLC(JDT, SynLC, LCm)
end