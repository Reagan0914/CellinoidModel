%%  Make the modification between cg in Geosystem and in Common Math system
%   eg: 
%       >> cg_Geo = cgModify(cg_Math, 'G');
%       >> cg_Math = cgModify(cg_Geo, 'M');
%   Edited by LUXP
%   Date: 2014-08-09

function cg = cgModify(cg,char)
    if char == 'M'        
        cg(7) = cg(7) * pi/180;
        cg(8) = (90-cg(8)) * pi/180;
        cg(9) = 2*pi * 24 /cg(9);           %Omega = 2pi/prd in day = 2pi/prd in hour/24=2pi*24/prd, i.e. rad/day  So Phi = PHI0-Omega*jdt
        cg(10) = cg(10) *pi/180;        
    elseif char == 'G'
        cg(1:6) = abs(cg(1:6));
        cg(7) = cg(7) * 180/pi;
        cg(7) = mod(cg(7), 360);
        cg(8) = 90 - cg(8) * 180/pi;
        cg(8) = mod(cg(8)+90, 180)-90;
        cg(9) = 2*pi * 24 /cg(9);
        cg(10) = cg(10) *180/pi;      
        cg(10) = mod(cg(10), 360);
    end
end
