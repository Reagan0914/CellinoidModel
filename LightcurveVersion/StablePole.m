function [V,G]=StablePole(cg)
%%  This function is applied to analyze the stability of free Rotational axis  for Cellinoid Model
%   W.R.T. the six parameters
%      Calculate: Mass, Mass Center, Stable Pole Orientation for specified
%      cellinoid shape model with the parameters:  cg = [a1,a2,b1,b2,c1,c2]
%
%           >> [V,G] = StablePole(cg);
%           V will be the 3 new standard axes and G will be the Center of
%           Mass
%   Edited by LUXP
%   Date: 2015-03-06

if nargin == 0
    cg=[1, 0.8, 0.7,0.6, 0.5,0.4 ]';
%cg=3*cg;
end
% if exist('mexEigFunction.m','file') == 0
%     path('./Others',path);
% end
%% Calculating the Mass(M) of Cellinoid
M = (pi/6) * (cg(1)+cg(2)) * (cg(3)+cg(4)) * (cg(5)+cg(6));

%% Calculating the mass center(G) of Cellinoid
G = zeros(3,1);
coef_Mass = 3/8;
G(1) = coef_Mass * (cg(1)-cg(2));
G(2) = coef_Mass * (cg(3)-cg(4));
G(3) = coef_Mass * (cg(5)-cg(6));

a1=cg(1);a2=cg(2);b1=cg(3);b2=cg(4);c1=cg(5);c2=cg(6);
%% Calculating the stable rotational axes of Cellinoid
% B = zeros(3,3);
% A = zeros(3,3);
% B(1,1) = G(2)^2+G(3)^2;
% B(2,2) = G(1)^2+G(3)^2;
% B(3,3) = G(1)^2+G(2)^2;
% B(1,2) = -G(1)*G(2); B(2,1)=B(1,2);
% B(1,3) = -G(1)*G(3); B(3,1)=B(1,3);
% B(2,3) = -G(2)*G(3); B(3,2)=B(2,3);

% coef_A = pi/30;
% A(1,1) = coef_A*((a1+a2)*(b1^3+b2^3)*(c1+c2) + (a1+a2)*(b1+b2)*(c1^3+c2^3));
% A(2,2) = coef_A*((a1^3+a2^3)*(b1+b2)*(c1+c2) + (a1+a2)*(b1+b2)*(c1^3+c2^3));
% A(3,3) = coef_A*((a1^3+a2^3)*(b1+b2)*(c1+c2) + (a1+a2)*(b1^3+b2^3)*(c1+c2));
% A(1,2) = -(1/15)*(a1^2-a2^2)*(b1^2-b2^2)*(c1+c2); A(2,1)=A(1,2);
% A(1,3) = -(1/15)*(a1^2-a2^2)*(b1+b2)*(c1^2-c2^2); A(3,1)=A(1,3);
% A(2,3) = -(1/15)*(a1+a2)*(b1^2-b2^2)*(c1^2-c2^2); A(3,2)=A(2,3);

% C = A-M*B;
% [V,D] = eig(C);
%%  Conclusion:
%   C matrix will increase k^5 with the cg increases k
%   Therefore, It MUST be adjusted to make the a1 = 1, and others adjusted!!!

%%  Another Theoretical Forumla for C:
F=zeros(3,3);
F(1,1)=pi*(19*b1^2/1920 + 19*b2^2/1920 + 13*b1*b2/960 + 19*c1^2/1920 + 19*c2^2/1920 + 13*c1*c2/960);
F(2,2)=pi*(19*a1^2/1920 + 19*a2^2/1920 + 13*a1*a2/960 + 19*c1^2/1920 + 19*c2^2/1920 + 13*c1*c2/960);
F(3,3)=pi*(19*a1^2/1920 + 19*a2^2/1920 + 13*a1*a2/960 + 19*b1^2/1920 + 19*b2^2/1920 + 13*b1*b2/960);

F(1,2)=(3*pi/128 - 1/15) * (a1-a2) * (b1-b2);
F(2,1)=F(1,2);
F(1,3)=(3*pi/128 - 1/15) * (a1-a2) * (c1-c2);
F(3,1)=F(1,3);
F(2,3)=(3*pi/128 - 1/15) * (b1-b2) * (c1-c2);
F(3,2)=F(2,3);
 
F=(a1+a2)*(b1+b2)*(c1+c2)*F;
%[V,D] = eig(F);
[V,~] = mexEigFunction(F);  %for faster speed and Right-Hand-System V
end