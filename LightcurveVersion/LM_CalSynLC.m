%%  @LM_CalSynLC(cg)
%   Description:
%       Invoked by getAlphaBeta.m, to calculate the synthetic lightcurve
%       from cg
%       
%   Edited by LUXP
%   Date: 2016-10-06
function SynLC = LM_CalSynLC(cg)
global JDT  E0 E  LebData LebInd

%%%%%%%%%%%%    Calculate the E and E0 in asteroid-centric frame
%%  Calculate Stable Pole and Center of Mass
[Q,G] = StablePole(cg(1:6));

%%  Calculate Rotation Matrix for pole in ecliptic frame
Lam = cg(7); Beta=cg(8);
cb = cos(Beta); cl = cos(Lam); sb = sin(Beta); sl = sin(Lam);
RotatePole = [cb*cl, cb*sl, -sb; -sl, cl, 0; sb*cl, sb*sl, cb];

%%  Calculate New E and E0 in asteroid-center frame and Generate Synthetic Lightcurves
%Calculate Normal vector and facet area for Lebedev Method
[NormalVec, Gw] = CalNormalVector(LebData, LebInd, cg);
Ndegree = size(LebData,1);
TotalPoints = size(JDT,1);

%%  Replace the comment codes by a C-Mex function to accelerate the speed
SynLC = LM_CalSynLC_mex(NormalVec, Gw, E0, E, JDT, Q, RotatePole, G, cg(9), cg(10));
% SynLC = zeros(TotalPoints,1);
% for i = 1:TotalPoints
%     %%   derive New E and E0
%     e0 = E0(:,i); e = E(:,i);
%     Phi = cg(10) - cg(9)*(JDT(i)-JDT(1));
%     cp = cos(Phi); sp = sin(Phi);
%     RotatePhase = [cp, sp, 0; -sp, cp, 0; 0 0 1];
%     
%     M = Q*RotatePhase*RotatePole;
%     e=M*e + G; e0 = M*e0 + G;
%     AsterE = e/norm(e); AsterE0 = e0/norm(e0);
%     
%     %% Compute the total brightness integration W.R.T. Normal vector, Facet area and E,E0
%     FBrt=zeros(Ndegree,1);   %Facet Brightness
%     for f=1:Ndegree
%             mu = AsterE(1)*NormalVec(1,f) + AsterE(2)*NormalVec(2,f) +AsterE(3) * NormalVec(3,f);
%             mu0 = AsterE0(1)*NormalVec(1,f) + AsterE0(2)*NormalVec(2,f) +AsterE0(3) * NormalVec(3,f);
%             if mu>0 && mu0>0
%                 FBrt(f)=(mu*mu0*(0.1+1/(mu+mu0)));
%             end
%     end
%     SynLC(i) = sum(FBrt.*Gw);
% end
end