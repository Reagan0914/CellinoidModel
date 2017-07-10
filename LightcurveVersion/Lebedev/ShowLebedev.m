%% @ShowLebedev
%  Description:
%     Show the Lebedev Points and convex hull distribution 
%     Compare the Volumns derived by Lebedev and 4*pi*r^3/3
%  How To Use:
%   >> ShowLebedev(N)  % N=1 for 230 points, N=2...
%   >> ShowLebedev   % Show the animation of convexhull, corresoponding
%   to N or points
%   
%   Edited by LUXP
%   Date: 2016-09-20
function ShowLebedev(N)
Degree = [230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354];
if nargin() == 1
    if N < 0 || N>size(Degree, 2)
        error('Input the right N!');
    end
    leb = getLebedevSphere(Degree(N));
    x = leb.x; y = leb.y; z = leb.z;
    dt=delaunayTriangulation([x,y,z]);
    [K V] = convexHull(dt);
    figure(101);
    trisurf(K,x,y,z);
    axis equal;
    colormap gray;
    title(['Distribution of Lebedev points(Vertices) (N=',num2str(Degree(N)),')'],'fontsize',14,'fontweight','bold');
    VolInfo = ['Volume: ',num2str(V),' Sphere: ',num2str(4*pi/3)];
    text(1,1,0.8,VolInfo,'fontsize',14,'fontweight','bold');
    fprintf('For Lebedev N=%d, Facets Num: %d, Vertex Num: %d\n', Degree(N), size(K,1), size(x,1));
else
    figure(102);
    for i=1:size(Degree,2)
        leb = getLebedevSphere(Degree(i));
        x = leb.x; y = leb.y; z = leb.z;
        dt=delaunayTriangulation([x,y,z]);
        [K V] = convexHull(dt);
        trisurf(K,x,y,z);
        axis equal;
        colormap gray;
        title(['Distribution of Lebedev points(Vertices) (N=',num2str(Degree(i)),')'],'fontsize',14,'fontweight','bold');
        VolInfo = ['Volume: ',num2str(V),' Sphere: ',num2str(4*pi/3)];
        text(1,1,0.8,VolInfo,'fontsize',14,'fontweight','bold');
        pause;
        close all;
    end

end
end