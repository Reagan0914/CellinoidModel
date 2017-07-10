%%  Test stability of eig.c function
%   Eg:
%   >> profile on; [DiffEigValueNum, NotRHS_eig, NotRHS_mex] = TestEigFunction(50000);profile viewer
%
%   Results:
%       Matlab Eig Function could not be rely on!! While mexEigFunction
%       could get right eige values and right eigen vectors with Right Hand
%       Rules and save 60% computational-time
%
%   Edited by LUXP
%   Date: 2015-06-25

%%  @File: TestEigC.c 
%       >cc -o test TestEigC.c
%       >./test
%   eejcb_luxp: generate 3 eigen vectors W.R.T. 3 eigen values from
%   smallest to largest, and abide by Right Hand System, 
%   i.e., 3 unit eigen vectors will be X,Y,Z-axis directions 

%%  @File: mexEigFunction.c
%       >mex mexEigFunction.c
%   generate a Matlab-invoked function following TestEigC.c
%       >[V,D] = mexEigFunction(A)
%   similar to 
%       >[V1,D1] = eig(A)
%   But V will abide by Right Hand Rule, while V1 not.

function [DiffEigValueNum, NotRHS_eig, NotRHS_mex] = TestEigFunction(TotalTestNum)
if nargin == 0
    TotalTestNum = 10;
end
DiffEigValueNum = 0;
NotRHS_eig = 0; %check whether abide by Righ Hand Rule
NotRHS_mex = 0;
for i=1:TotalTestNum
    %%  Generate Matrix F with random semi-axes ABC
    F = GenerateRandomF;
    
    %%  Calculate By eig
    [V,D] = MatlabEigFunction(F);
    
    %%  Calculate by mexEigFunction 
    [V1, D1] = mexEigFunction(F);
    
    %%  Check difference of two methods and check the right hand system
    [DiffEigValueNum, NotRHS_eig, NotRHS_mex] = CheckEigenPairs(DiffEigValueNum, NotRHS_eig, NotRHS_mex, D, V, D1, V1);
end
end

%%  invoke Matlab function:eig
function [V,D] = MatlabEigFunction(F)
    [V, D] = eig(F);
end

%%  Check difference of two methods and check the right hand system
function [DiffEigValueNum, NotRHS_eig, NotRHS_mex] = CheckEigenPairs(DiffEigValueNum, NotRHS_eig, NotRHS_mex, D, V, D1, V1)
%  Check the difference of eigen values by two methods
if (norm(diag(D) - diag(D1)) > 1e-13)
    DiffEigValueNum = DiffEigValueNum + 1;
end

%  Check the right hand system
if (norm(cross(V(:,1), V(:,2)) - V(:,3)) > 1e-13)
    NotRHS_eig = NotRHS_eig + 1;
end

if (norm(cross(V1(:,1), V1(:,2)) - V1(:,3)) > 1e-13)
    NotRHS_mex = NotRHS_mex + 1;    
end
end

%%  Generate Matrix F with random semi-axes ABC
function F=GenerateRandomF
    ABC = rand(6,1);    
    a1=ABC(1);a2=ABC(2); b1=ABC(3);b2=ABC(4); c1=ABC(5);c2=ABC(6);

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
end

