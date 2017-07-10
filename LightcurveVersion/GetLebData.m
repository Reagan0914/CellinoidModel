%% @GetLebData
%  Description:
%       Pre-Calculation for Lebedev vector and quickly index the vector for
%       each octant
%               
%       if i ==1 VecInd = 1:LebInd(i); else VecInd = (LebInd(i-1)+1):LebInd(i)        for all vectors in i-th octant
%
%  How To Use:
%   >> [LebData, LebInd] = GetLebData(N);       %N = 590,5810 from Lebedev 
%   
%   Edited by LUXP
%   Date: 2016-09-26

function [LebData, LebInd] = GetLebData(Ndegree)
if nargin == 0
    Ndegree = 590;
end

if exist('getLebedevSphere.m','file') == 0
    path('./Lebedev',path);
end
leb=getLebedevSphere(Ndegree);
xx=leb.x;
yy=leb.y;
zz=leb.z;
ww=leb.w;

OctantIndex = zeros(Ndegree,1);

for i = 1:Ndegree
    x = xx(i); y = yy(i); z = zz(i);
   if x>=0 && y>=0 && z>=0              % 1st octant
       OctantIndex(i) = 1;
   elseif x<0 && y>=0 && z>=0           % 2 octant
       OctantIndex(i) = 2;
   elseif x<0 && y < 0 && z>=0         % 3 octant
       OctantIndex(i) = 3;
   elseif x>=0 && y<0 && z>=0            % 4 octant
       OctantIndex(i) = 4;
   elseif x>=0 && y>=0 && z<0           % 5 octant
       OctantIndex(i) = 5;
   elseif x<0 && y>=0 && z<0            % 6 octant
       OctantIndex(i) = 6;
   elseif x<0 && y<0 && z<0            % 7 octant
       OctantIndex(i) = 7;
   elseif x>=0 && y<0 && z<0             % 8 octant
       OctantIndex(i) = 8;
   end
end
LebData = [xx, yy, zz, ww, OctantIndex];
LebData = sortrows(LebData,5);
LebInd = zeros(8,1);
for i=1:8
    LebInd(i) = sum(LebData(:,5) == i);
end
LebInd = cumsum(LebInd);
end