%% compare Lebedev Quadrature and Triangularizaiton for Cellinoid
%
%   Edited by LUXP
%   Date: 2016-09-21

function CMP_Leb_Tri
a1=5;a2=4;b1=4;b2=3;c1=3.5;c2=2.5;
%a1=5;a2=5;b1=4;b2=5;c1=3;c2=3;      %ellipsoid
%a1=4;a2=4;b1=4;b2=4;c1=4;c2=4; %sphere
ShowCellinoid(a1,a2,b1,b2,c1,c2);

%%  Calulate the standard area as benchmark
StdArea = CalSurfaceAreaStd(a1,a2,b1,b2,c1,c2);
%%  Calculate By Triangularization
for Nrows = 6:50
    AreaTri(Nrows-5) = CalSurfaceAreaTri(a1,a2,b1,b2,c1,c2,Nrows);
    Nfacets(Nrows -5) = 8*Nrows^2;
end
figure,
set(gcf,'Position',[100,200,725,770]);
plot(Nfacets, AreaTri, '-r*'); 
hold on; plot(Nfacets, StdArea*ones(size(Nfacets)), '-k', 'LineWidth',3); % for benchmark
title('Surface Area of Cellinoid By Triangulation','fontsize',22);
xlabel('Number of Triangular Facets','fontsize',22);
set(gca,'FontName','Times New Roman','FontSize',16,'fontweight','bold');
ylabel('Surface Area of Cellinoid','fontsize',22);

%%  Calculate By Lebedev
DEGREE=[6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, ...
    350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, ... 
 3470, 3890, 4334, 4802, 5294, 5810];
for Nrows = 6:32
    AreaLeb(Nrows-5) = CalSurfaceAreaLeb(a1,a2,b1,b2,c1,c2,DEGREE(Nrows));
    NfacetsLeb(Nrows -5) = DEGREE(Nrows);
end
figure,
set(gcf,'Position',[100,200,725,770]);
plot(NfacetsLeb, AreaLeb, '-bo'); 
hold on; TargetNum = 12;plot(NfacetsLeb(TargetNum),AreaLeb(TargetNum),'*r'); %for N=17, i.e. 590 facets in Lebedev
hold on; plot(Nfacets, StdArea*ones(size(Nfacets)), '-k', 'LineWidth',3); % for benchmark
title('Surface Area of Cellinoid By Lebedev','fontsize',22);
xlabel('Number of Triangular Facets','fontsize',22);
set(gca,'FontName','Times New Roman','FontSize',16,'fontweight','bold');
ylabel('Surface Area of Cellinoid','fontsize',22);

%% Show the comparison figure
figure,
set(gcf,'Position',[100,200,725,770]);
plot(Nfacets, AreaTri, '-r*');hold on;
plot(NfacetsLeb, AreaLeb, '-bo');
hold on; TargetNum = 12;plot(NfacetsLeb(TargetNum),AreaLeb(TargetNum),'*r'); %for N=17, i.e. 590 facets in Lebedev
hold on; plot(Nfacets, StdArea*ones(size(Nfacets)), '-k', 'LineWidth',3); % for benchmark
title('Comparision between Triangular and Lebedev','fontsize',22);
xlabel('Number of Triangular Facets','fontsize',22);
set(gca,'FontName','Times New Roman','FontSize',16,'fontweight','bold');
ylabel('Surface Area of Cellinoid','fontsize',22);
legend('Triangularization','Lebedev', 'Best Leb:N=590','Benchmark:Standard Area');

end


%%  Calculate Surface Area of Cellinoid by Trianularization
function Area = CalSurfaceAreaTri(a1,a2,b1,b2,c1,c2,Nfacets)

aa = [a1,a2,a2,a1,a1,a2,a2,a1]';
bb = [b1,b1,b2,b2,b1,b1,b2,b2]';
cc = [c1,c1,c1,c1,c2,c2,c2,c2]';
[ifp,Theta,Phi] = TriFacetsOctant(Nfacets);
TotalFacets = size(ifp,2);   %total facets of one octant
Area = 0;
for i=1:8    
    x=aa(i)*sin(Theta(:,i)).*cos(Phi(:,i));
    y=bb(i)*sin(Theta(:,i)).*sin(Phi(:,i));
    z=cc(i)*cos(Theta(:,i));
    for f = 1:TotalFacets
        vec1 = zeros(3,1);
        vec2 = zeros(3,1);
        vec1(1)=x(ifp(2,f))-x(ifp(1,f));
        vec1(2)=y(ifp(2,f))-y(ifp(1,f));
        vec1(3)=z(ifp(2,f))-z(ifp(1,f));
        vec2(1)=x(ifp(3,f))-x(ifp(1,f));
        vec2(2)=y(ifp(3,f))-y(ifp(1,f));
        vec2(3)=z(ifp(3,f))-z(ifp(1,f));
        NormalVec = cross(vec1,vec2);
        Area = Area + norm(NormalVec)/2;        
    end
end

end

%%  Calculate Surface Area of Cellinoid by Lebedev
function Area = CalSurfaceAreaLeb(a1,a2,b1,b2,c1,c2,Ndegree)
leb=getLebedevSphere(Ndegree);
xx=leb.x;
yy=leb.y;
zz=leb.z;
w=leb.w;

G = zeros(Ndegree,1);
NumNodes = Ndegree;
% if type==1
%     G=(a*b*c./((a*x).^2+(b*y).^2+(c*z).^2)).^2;  %Kaasalainen's Curvature Function
% else
for i = 1:Ndegree
    x = xx(i); y = yy(i); z = zz(i);
   if x>=0 && y>=0 && z>=0              % 1st octant
       a = a1; b = b1; c= c1;
   elseif x<0 && y>=0 && z>=0           % 2 octant
       a = a2; b = b1; c= c1;
   elseif x<0 && y < 0 && z>=0         % 3 octant
       a = a2; b = b2; c= c1;
   elseif x>=0 && y<0 && z>=0            % 4 octant
       a = a1; b = b2; c= c1;
   elseif x>=0 && y>=0 && z<0           % 5 octant
       a = a1; b = b1; c= c2;
   elseif x<0 && y>=0 && z<0            % 6 octant
       a = a2; b = b1; c= c2;
   elseif x<0 && y<0 && z<0            % 7 octant
       a = a2; b = b2; c= c2;
   elseif x>=0 && y<0 && z<0             % 8 octant
       a = a1; b = b2; c= c2;
   end
   G(i) = a*b*c*sqrt((x/a)^2+(y/b)^2+(z/c)^2);     %LUXP's Curvature Function
end
% end
Area=sum(w.*G);
end

%%  Calcualte the standard area of Cellinoid
function Area = CalSurfaceAreaStd(a1,a2,b1,b2,c1,c2)
aa = [a1,a2,a2,a1,a1,a2,a2,a1]';
bb = [b1,b1,b2,b2,b1,b1,b2,b2]';
cc = [c1,c1,c1,c1,c2,c2,c2,c2]';
Area = 0;
for i = 1:8
    a = aa(i); b = bb(i); c = cc(i);
    Sarea=funEllipsoidSurfaceArea(a,b,c);
    Area = Area + Sarea/8;
end
end
%%  Calculate the standard area of Ellipsoid 
function Sarea=funEllipsoidSurfaceArea(a,b,c)
% this function can get the surface area of Ellipsoid with three semi-axies
% a,b,c in a fast way by MATLAB internal function.
% Ref: http://www.matlabsky.com/thread-11891-1-1.html
% Modified by LUXP
% DATE:2011-8-24
syms x y z
F=x^2/a^2+y^2/b^2+z^2/c^2-1;
dzdx=-diff(F,x)/diff(F,z);
dzdy=-diff(F,y)/diff(F,z);
z=solve(F,z);
z=z(1);
fxy=subs(sqrt(1+dzdx^2+dzdy^2));
yleft=-subs(sqrt(1-x^2/a^2)*b);
yright=subs(sqrt(1-x^2/a^2)*b);
Fxy1=char(fxy);
Fxy1=strrep(Fxy1,'/','./');
Fxy1=strrep(Fxy1,'^','.^');
Fxy1=strrep(Fxy1,'*','.*');
Fxy1=eval(['@(x,y)' (Fxy1)]);
yright=char(yright);
yright=strrep(yright,'^','.^');
yright=strrep(yright,'/','./');
yright=strrep(yright,'*','.*');
yright=eval(['@(x)' (yright)]);
S1=2*quad2d(Fxy1,-a,a,0,yright);
Sarea=S1*2;
end


%% Show Cellinoid Shape Model
function ShowCellinoid(a1,a2,b1,b2,c1,c2)
aa=[a1,a2,a2,a1,a1,a2,a2,a1]';
bb=[b1,b1,b2,b2,b1,b1,b2,b2]';
cc=[c1,c1,c1,c1,c2,c2,c2,c2]';
[ifp,Theta,Phi] = TriFacetsOctant(12);
figure,
for i=1:8    
    x=aa(i)*sin(Theta(:,i)).*cos(Phi(:,i));
    y=bb(i)*sin(Theta(:,i)).*sin(Phi(:,i));
    z=cc(i)*cos(Theta(:,i));      
%     plot3(x,y,z,'bo');
%     pause(1);
    trisurf(ifp',x,y,z, x);
    hold on;
end
    quiver3(0,0,0,1,0,0,5,'k','filled','LineWidth',2);
    quiver3(0,0,0,0,1,0,4,'k','filled','LineWidth',2);
    quiver3(0,0,0,0,0,1,3,'k','filled','LineWidth',2);
    
    text(5,0,0,'X','fontsize',14,'fontweight','bold');
    text(0,4,0,'Y','fontsize',14,'fontweight','bold');
    text(0,0,3,'Z','fontsize',14,'fontweight','bold');
    view(0,0);axis off;axis equal;colormap gray;
    title('Cellinoid Shape','fontsize',22,'fontweight','bold');
end

function [ifp,Theta,Phi] = TriFacetsOctant(nrows)
%%  Generate the eight octants data
%	Make the index of all nodes distributing in one octant surface of unit sphere
%   Total number of nodes of each octant with nrows spliting is "(nrows+1)*(nrows+2)/2"
%
%   eg:
%       >> [ifp,Theta,Phi] = TriFacetsOctant(12)
%   Edited by LUXP
%   Date: 2014-08-09

TotalNodes=(nrows+1)*(nrows+2)/2;
nnod=1; % index number of all vertices in one octant from 1 to TotalNodes
nod=zeros(TotalNodes,TotalNodes);

nod(1,1)=nnod;
for row=2:(nrows+1)
    for col=1:row
        nnod=nnod+1;
        nod(row,col)=nnod;       
    end
end

%   set index of three nodes for every trifacet, triangularization
ntri=0;  %index of trifacet with a total number "nrows^2" in one octant
ifp=zeros(3,nrows^2);

for row=1:nrows            
        ntri=ntri+1;
        ifp(1,ntri)=nod(row,1);
        ifp(2,ntri)=nod(row+1,1);
        ifp(3,ntri)=nod(row+1,2);
        
        for col=1:(row-1)
            ntri=ntri+1;
            ifp(1,ntri)=nod(row+1,col+1);
            ifp(2,ntri)=nod(row,col+1);
            ifp(3,ntri)=nod(row,col);
            
            ntri=ntri+1;
            ifp(1,ntri)=nod(row,col+1);
            ifp(2,ntri)=nod(row+1,col+1);
            ifp(3,ntri)=nod(row+1,col+2);            
        end
end

Theta=zeros(TotalNodes,8);
Phi=zeros(TotalNodes,8);
%generate the above part
[Theta(:,1),Phi(:,1)]=GenerateOctant(nrows,0,pi/2,0,pi/2);
[Theta(:,2),Phi(:,2)]=GenerateOctant(nrows,0,pi/2,pi/2,pi);
[Theta(:,3),Phi(:,3)]=GenerateOctant(nrows,0,pi/2,pi,1.5*pi);
[Theta(:,4),Phi(:,4)]=GenerateOctant(nrows,0,pi/2,1.5*pi,2*pi);
%generate the below part
[Theta(:,5),Phi(:,5)]=GenerateOctant(nrows,pi,pi/2,0,pi/2);
[Theta(:,6),Phi(:,6)]=GenerateOctant(nrows,pi,pi/2,pi/2,pi);
[Theta(:,7),Phi(:,7)]=GenerateOctant(nrows,pi,pi/2,pi,1.5*pi);
[Theta(:,8),Phi(:,8)]=GenerateOctant(nrows,pi,pi/2,1.5*pi,2*pi);

end

%% Generate the Decartician Coordinates for one octant with theta in [Theta1, Theta2] and phi in [Phi1, Phi2]
function [vTheta,vPhi]=GenerateOctant(nrows,Theta1,Theta2,Phi1,Phi2)

%Total number of nodes of each octant with nrows spliting is "(nrows+1)*(nrows+2)/2"
TotalNodes=(nrows+1)*(nrows+2)/2;

% spherical coordinates of nodes in the octant
vTheta=zeros(TotalNodes,1);vPhi=zeros(TotalNodes,1);

dth=(Theta2-Theta1)/nrows;
indNodes=1; % first node is the polar position 
vTheta(1)=Theta1;
for i=1:nrows
    dph=(Phi2-Phi1)/i;
    for j=0:i
        indNodes=indNodes+1;
        vTheta(indNodes)=Theta1+i*dth;
        vPhi(indNodes)=Phi1+dph*j;
    end
end
end