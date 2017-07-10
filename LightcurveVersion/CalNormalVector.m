%% @FileName
%  Description:
%       Calculate Normal Vector and area of small facets for Lebedev methods   
%
%  How To Use:
%   >> [NormalVec, G, w] = CalNormalVector(LebData, LebInd, cg);
%       following:  [LebData, LebInd] = GetLebData(590);
%   
%   Edited by LUXP
%   Date: 2016-09-26
function [NormalVec, Gw] = CalNormalVector(LebData, LebInd, cg)
Ndegree = size(LebData,1);
%   Calculate NormalVector of each facets by Lebedev Quadrature
a1 = cg(1); a2=cg(2); b1=cg(3); b2=cg(4); c1=cg(5); c2=cg(6);
AA = [a1,a2,a2,a1,a1,a2,a2,a1]';
BB = [b1,b1,b2,b2,b1,b1,b2,b2]';
CC = [c1,c1,c1,c1,c2,c2,c2,c2]';
NormalVec = zeros(3, Ndegree);
w = LebData(:, 4);

%%  For Kaasalainen's Scattering Law, Disk Function: = mu*mu0*(0.1+1/mu+mu0)
G = zeros(Ndegree, 1);
for octInd = 1:8
    if octInd ==1 
        VecInd = 1:LebInd(1); 
    else
        VecInd = (LebInd(octInd-1)+1):LebInd(octInd);
    end
    a = AA(octInd); b = BB(octInd); c = CC(octInd);
    x = LebData(VecInd,1); y = LebData(VecInd, 2); z = LebData(VecInd, 3);
    GTemp = a*b*c*sqrt((x/a).^2+(y/b).^2+(z/c).^2); 
    NormalVec(:,VecInd) = [b*c*x./GTemp, a*c*y./GTemp, a*b*z./GTemp]';
    G(VecInd) = GTemp;
end
Gw = G.*w;
%% end

% %%  For Lommel-Seeliger Scattering law: Disk Function: = mu*mu0/(mu+mu0)
% for octInd = 1:8
%     if octInd ==1 
%         VecInd = 1:LebInd(1); 
%     else
%         VecInd = (LebInd(octInd-1)+1):LebInd(octInd);
%     end
%     a = AA(octInd); b = BB(octInd); c = CC(octInd);
%     x = LebData(VecInd,1); y = LebData(VecInd, 2); z = LebData(VecInd, 3);
%     NormalVec(:,VecInd) = [b*c*x, a*c*y, a*b*z]';
% end
% Gw = w;
% %%  End
end