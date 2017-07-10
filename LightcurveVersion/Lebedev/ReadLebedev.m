%%  Readin lebedev and generate the following format file for DAMIT/convexinv_luxp
%
%       leb.x, leb.y, leb.z, theta, phi in rad
%
%   Edited by LUXP
function ReadLebedev(Num)
if nargin == 0
    Num =590;
end
Degree = [230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354];
if(find(Degree==Num))
    leb = getLebedevSphere(Num);
    fid = fopen(['Lebedev',num2str(Num)],'w+');
    for i=1:Num
        x = leb.x(i);y = leb.y(i); z = leb.z(i); w = leb.w(i);
        %%  generate by LUXP code, to get 'Leb590'
%         theta = acos(z); 
%         if abs(z) == 1
%             phi = 0;
%         else
%             if y >= 0 
%                 phi = acos(x/sqrt(1-z^2));
%             else
%                 phi = 2*pi - acos(x/sqrt(1-z^2));
%             end
%         end
        %% generate by Matlab inline code, to get Lebedev590
        [phi, theta, r] = cart2sph(x,y,z);
        theta =pi/2 - theta;
        phi = mod(phi, 2*pi);
        fprintf(fid,'%10.7f %10.7f %10.7f %10.8f %10.8f %10.8f\n', x, y, z, w, theta, phi);
    end
    fclose(fid);
else
    disp('No this Degree in Lebedev!');
end
end