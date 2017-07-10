%% @FileName
%  Description:
%       Calculate synthetic brightness, given the E, E0 and other
%       parameters.
%
%  How To Use:
%   >> SynTotalBrt = CalBrtPoint(Ndegree, E, E0, NormalVec, w, G)
%   
%   Edited by LUXP
%   Date: 2016-09-26
function SynTotalBrt = CalBrtPoint(Ndegree, E, E0, NormalVec, Gw)
%% Compute the total brightness integration W.R.T. Normal vector, Facet area and E,E0
FBrt=zeros(Ndegree,1);   %Facet Brightness
for i=1:Ndegree
        mu = E(1)*NormalVec(1,i) + E(2)*NormalVec(2,i) +E(3) * NormalVec(3,i);
        mu0 = E0(1)*NormalVec(1,i) + E0(2)*NormalVec(2,i) +E0(3) * NormalVec(3,i);
        if mu>0 && mu0>0
            FBrt(i)=(mu*mu0*(0.1+1/(mu+mu0)));
        end
end
SynTotalBrt = sum(FBrt.*Gw);
end