function LC_CurvVS(a,b,c)
% This function is used to generate the Fig2 in the paper:Ellipsoid
%
% Compare Triangulation and Lebedev Quadrature with increasing serials
% Evaluation points for
% Computing the surface area of Ellipsoid with different a,b,c
%
% Run Ref:-------------------------------- 
% a=8;b=7;c=6;
% Fig2_CurvVS(a,b,c)
% axis([1,10^5,10^-12,10^0])
% set(gca,'XTick',[10^0,10,10^2,10^3,10^4,10^5])
% set(gca,'YTick',[10^-12,10^-10,10^-8,10^-6,10^-4,10^-2,10^0])
% a=10;b=2;c=1.5;
% Fig2_CurvVS(a,b,c)
% axis([1,10^5,10^-10,10^2])
% set(gca,'XTick',[10^0,10,10^2,10^3,10^4,10^5])
% set(gca,'YTick',[10^-10,10^-8,10^-6,10^-4,10^-2,10^0,10^2])

% EDITED by LUXP:
% Date: 2012-7-26
%
 figure,
 [x, y, z] = ellipsoid(0,0,0,a,b,c,30);
surfl(x, y, z);
title('Ellipsoid Shape','fontsize',22,'fontweight','bold');
%MaxNrows=50;
DEGREE=[6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, ...
    350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, ... 
 3470, 3890, 4334, 4802, 5294, 5810];
MaxDegree=32;
MaxNrows=ceil(sqrt((DEGREE(MaxDegree)-2)/4));
Area_tri=zeros(MaxNrows,1);
NumNodes_tri=zeros(MaxNrows,1);

for i=1:MaxNrows
    [NumNodes_tri(i),Area_tri(i)]=EllipsoidAreaTri(a,b,c,i);
end
%standardNumNodes=NumNodes_tri(MaxNrows10);
%[stdNumNode,standardArea]=EllipsoidAreaTri(a,b,c,250,1);
standardArea=funEllipsoidSurfaceArea(a,b,c);

Area_leb=zeros(MaxDegree,1);
Area_leb_xp=zeros(MaxDegree,1); % Curvature of Formular(27)
NumNodes_leb=zeros(MaxDegree,1);

for i=1:MaxDegree
    deg=DEGREE(i);
    [NumNodes_leb(i),Area_leb(i)]=EllipsoidAreaLeb(a,b,c,deg,1);  
    [NumNodes_leb(i),Area_leb_xp(i)]=EllipsoidAreaLeb(a,b,c,deg,2);  
end
Res_tri=standardArea-Area_tri;
Res_tri=Res_tri/standardArea;
Res_leb=standardArea-Area_leb;
Res_leb=Res_leb/standardArea;
Res_leb_xp=standardArea-Area_leb_xp;
Res_leb_xp=Res_leb_xp/standardArea;

figure,
set(gcf,'Position',[100,500,525,470]);
loglog(NumNodes_tri,abs(Res_tri),'d');
xlabel('Number of Triangular Facets','fontsize',22);
ylabel('Error','fontsize',22);
set(gca,'FontName','Times New Roman','FontSize',22);
title('Error Curves: Triangulation(*) Lebedev(o)','fontsize',14,'fontweight','bold');
axis square;
hold on;
loglog(NumNodes_leb,abs(Res_leb),'o');
loglog(NumNodes_leb,abs(Res_leb_xp),'*');
legend('Triangulation', 'Kaasa Lebedev', 'LUXP Lebedev');
%annotation(103,'textbox',...
%    [0.6093984375 0.73644578313253 0.15915625 0.170180722891569],...
%    'String',{['a=',num2str(a),', b=',num2str(b),', c=',num2str(c)],['stdArea=',num2str(standardArea)]},...
%    'FitBoxToText','off');
end

function [NumNodes,totalArea]=EllipsoidAreaTri(a,b,c,nrows)
% trianglization
% this is copy from C code of Prof. Kaasalainen and Durech's program 
% Transformed by LUXP
nod=zeros(2*(nrows+1),4*(nrows+1));
nnod=1; % index number of all vertices
nod(1,1)=nnod;
for row=2:(nrows+1)
    for col=1:(4*row-4)
        nnod=nnod+1;
        nod(row,col)=nnod;
        if col==1 
            nod(row,4*row-3)=nnod;
        end
    end
end
for row=nrows:-1:1
     for col=1:(4*row-4)
         nnod=nnod+1;
         nod(2*(nrows+1)-row,col)=nnod;
        if col==1 
            nod(2*(nrows+1)-row,4*row-3)=nnod;
        end
     end
end
nnod=nnod+1;
nod(2*nrows+1,1)=nnod;

% set index of three nod point for every trifacet
ntri=0;
ifp=zeros(3,8*nrows^2);
for row=1:nrows
    for oct=1:4
        col0=(oct-1)*row+1;
        ntri=ntri+1;
        ifp(1,ntri)=nod(row,col0-oct+1);
        ifp(2,ntri)=nod(row+1,col0);
        ifp(3,ntri)=nod(row+1,col0+1);
        
        for col=col0+1:row+col0-1
            ntri=ntri+1;
            ifp(1,ntri)=nod(row+1,col);
            ifp(2,ntri)=nod(row,col-oct+1);
            ifp(3,ntri)=nod(row,col-oct);
            
            ntri=ntri+1;
            ifp(1,ntri)=nod(row,col-oct+1);
            ifp(2,ntri)=nod(row+1,col);
            ifp(3,ntri)=nod(row+1,col+1);            
        end
    end
end

for row=nrows+1:2*nrows
    for oct=1:4
        col0=(oct-1)*(2*nrows-row)+1;
        ntri=ntri+1;
        ifp(1,ntri)=nod(row+1,col0);
        ifp(2,ntri)=nod(row,col0+oct);
        ifp(3,ntri)=nod(row,col0+oct-1);
        for col=col0+1:col0+2*nrows-row
            ntri=ntri+1;
            ifp(1,ntri)=nod(row+1,col);
            ifp(2,ntri)=nod(row,col+oct-1);
            ifp(3,ntri)=nod(row+1,col-1);
            
            ntri=ntri+1;
            ifp(1,ntri)=nod(row+1,col);
            ifp(2,ntri)=nod(row,col+oct);
            ifp(3,ntri)=nod(row,col+oct-1);    
        end        
    end
end

% generate the x,y,z coordinate and plot trisurf
dth=pi/2/nrows;
t=zeros(8*nrows^2,1);
f=zeros(8*nrows^2,1);
k=1;
for i=1:nrows
    dph=pi/2/i;
    for j=0:(4*i-1)
        k=k+1;
        t(k)=i*dth;
        f(k)=j*dph;
    end
end
for i=(nrows-1):-1:1
    dph=pi/2/i;
    for j=0:(4*i-1)
        k=k+1;
        t(k)=pi-i*dth;
        f(k)=j*dph;
    end
end
ndir=k+1;
t(ndir)=pi;

t=t(1:ndir);
f=f(1:ndir);

x=sin(t).*cos(f);
y=sin(t).*sin(f);
z=cos(t);

 figure(111),
 trisurf(ifp',x,y,z);
 title('Unit Sphere with Lebedev','fontsize',22,'fontweight','bold');
 axis equal;

%    figure(101);
    x=a*x;
    y=b*y;
    z=c*z;
%     trisurf(ifp',x,y,z);
%     title(['tri---Ellipsoid NumNodes=' num2str(nnod)],'fontsize',14,'fontweight','bold');
%     axis equal;

%%%%%%%%%%%%compute the unit normal vector and area of every facet
Nor=zeros(3,ntri);      %unit normal vector of every facet
Darea=zeros(ntri,1);  % area of Deluxity triangle facet
at=zeros(ntri,1);   % theta and phi of normal vector
af=zeros(ntri,1);
vec1=zeros(3,1);
vec2=vec1;
for i=1:ntri
    vec1(1)=x(ifp(2,i))-x(ifp(1,i));
    vec1(2)=y(ifp(2,i))-y(ifp(1,i));
    vec1(3)=z(ifp(2,i))-z(ifp(1,i));
    vec2(1)=x(ifp(3,i))-x(ifp(1,i));
    vec2(2)=y(ifp(3,i))-y(ifp(1,i));
    vec2(3)=z(ifp(3,i))-z(ifp(1,i));
    Nor(:,i)=cross(vec1,vec2);
    clen=norm(Nor(:,i));
    Nor(:,i)=Nor(:,i)/clen;
    Darea(i)=clen/2;
    at(i)=acos(Nor(3,i));
    af(i)=acos(Nor(1,i)/sqrt(Nor(1,i)^2+Nor(2,i)^2));
    if Nor(2,i)<0
        af(i)=2*pi-af(i);
    end
end
NumNodes=ntri;
totalArea=sum(Darea);
end

function [NumNodes,totalArea]=EllipsoidAreaLeb(a,b,c,degree,type)
leb=getLebedevSphere(degree);
x=leb.x;
y=leb.y;
z=leb.z;
w=leb.w;
NumNodes=degree;
if type==1
    G=(a*b*c./((a*x).^2+(b*y).^2+(c*z).^2)).^2;  %Kaasalainen's Curvature Function
else
    G=a*b*c*sqrt((x/a).^2+(y/b).^2+(z/c).^2);     %LUXP's Curvature Function
end
totalArea=sum(w.*G);
end

function Sarea=funEllipsoidSurfaceArea(a,b,c)
% this function can get the surface area of Ellipsoid with three semi-axies
% a,b,c in a fast way by MATLAB internal function.
% Ref: http://www.matlabsky.com/thread-11891-1-1.html
% Modified by LUXP
% DATE:2011-8-24

%???
syms x y z
F=x^2/a^2+y^2/b^2+z^2/c^2-1;
dzdx=-diff(F,x)/diff(F,z);%??
dzdy=-diff(F,y)/diff(F,z);
z=solve(F,z);
z=z(1);
fxy=subs(sqrt(1+dzdx^2+dzdy^2));%????

%???????????????????7.x?
%% ????????????????????????????????????????
% tic
% Fxy=eval(['@(x,y)' char(fxy)]);
% %??????????????????????????????????????????????
% [X,Y]=meshgrid(linspace(-a,a),linspace(-b,b));
% dS=arrayfun(Fxy,X,Y)*(X(1,2)-X(1))*(Y(2)-Y(1));%???arrayfun,????for??
% Z=X.^2/a^2+Y.^2/b^2;
% dS(Z>1)=0;
% S=sum(dS(:))%????
% toc

%% ????????????????????????????????quad2d)
%tic
yleft=-subs(sqrt(1-x^2/a^2)*b);
yright=subs(sqrt(1-x^2/a^2)*b);
% S=int(int(fxy,y,yleft,yright),x,-a,a)
Fxy1=char(fxy);
Fxy1=strrep(Fxy1,'/','./');
Fxy1=strrep(Fxy1,'^','.^');
Fxy1=strrep(Fxy1,'*','.*');
Fxy1=eval(['@(x,y)' (Fxy1)]);%????????????
yright=char(yright);
yright=strrep(yright,'^','.^');
yright=strrep(yright,'/','./');
yright=strrep(yright,'*','.*');
yright=eval(['@(x)' (yright)]);%??????????????
%??????????????
S1=2*quad2d(Fxy1,-a,a,0,yright);%???????????
%toc
Sarea=S1*2;
end
