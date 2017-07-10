%% @GetSunEarth
%  Description:
%       Calculate the Sun(E0) and Earth(E) direction at asteroid-centric coordinate system
%       Calibrated by: stable pole aixs, JDT,rotational period, Real Pole Orienatation(lam,beta)
%
%  How To Use:
%   >> [E0, E] = GetSunEarth(cg)
%   
%   Edited by LUXP
%   Date: 2016-09-26
function [E0, E, JDT, LCm] = GetSunEarth(cg)
if nargin==0
    cg  = [1, 0.9, 0.75, 0.75, 0.5, 0.5, 100, 50, 0, 120, -0.4, 1, 0.15]'; 
    cg = cgModify(cg, 'M');
end
%%  Fetch Kaasalainen's Format Data: 
%   LC:= Julian date, Relative Brt, Unit Sun, Unit Earth, [Error, Solar Phase in Radial]
[LC, LCm] = fReadData;
JDT = LC(:,1) - LC(1,1);
ObsBrt = LC(:,2);
E0= LC(:, 3:5)';
E = LC(:,6:8)';

%%  Calculate Stable Pole and Center of Mass
[Q,G] = StablePole(cg(1:6));

%%  Calculate Rotation Matrix for pole in ecliptic
Lam = cg(7); Beta=cg(8);
cb = cos(Beta); cl = cos(Lam); sb = sin(Beta); sl = sin(Lam);
RotatePole = [cb*cl, cb*sl, -sb; -sl, cl, 0; sb*cl, sb*sl, cb];

%%  Calculate New E and E0 in asteroid-center frame
for i = 1:size(JDT,1)
    e0 = E0(:,i); e = E(:,i);
    Phi = cg(10) - cg(9)*JDT(i);
    cp = cos(Phi); sp = sin(Phi);
    RotatePhase = [cp, sp, 0; -sp, cp, 0; 0 0 1];
    
    M = Q*RotatePhase*RotatePole;
    e=M*e + G; e0 = M*e0 + G;
    E(:,i) = e/norm(e); E0(:,i) = e0/norm(e0);
end

end