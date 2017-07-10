%% @LM_Inverse
%       Search the best-fit solution with the locally extrema by
%       Levenberg-Marquardt algorithm from the given intial test values.
%
%   How to Use:
%       > cg_bestfit = LM_Inverse(LC, cg)
%   Edited by LUXP
%   Date: 2016-10-06
function [a_cur, Chisq_cur]=LM_Inverse(a_cur,UpdateIndex)

%%  initial settings for LM
[A,B,Chisq_cur] =getAlphaBeta(a_cur, UpdateIndex); 
lam = 0.001;  
IterMax = 100; 
NDONE = 40;
DONE = 0;
tol = 1e-4;
Chisq_history = Chisq_cur;

TotalNumParas = size(a_cur,1);
TotalNumUpdateParas = sum(UpdateIndex);


for iter = 1:IterMax
%   Modify matrix \alpha A by muliplying its diagonal elements
    for j=1:TotalNumUpdateParas
        A(j,j) = A(j,j)*(1+lam);
    end
%   Compute the iteration step \delta a: da
    da_Ind = A\B;

%   Generate the right update da
    Ind = 1;
    da = zeros(TotalNumParas,1);
    for j=1:TotalNumParas
        if UpdateIndex(j) ==1
            da(j) = da_Ind(Ind);
            Ind = Ind + 1;
        end
    end
    [A_next, B_next, Chisq_next] = getAlphaBeta(a_cur+da, UpdateIndex);
    
    if (abs(Chisq_cur-Chisq_next) < max(tol, tol*Chisq_cur)) 
        DONE = DONE +1;
    end
    
    if Chisq_next < Chisq_cur   %success and update 
        a_cur = a_cur+da;
        Chisq_cur = Chisq_next;
        lam = lam * 0.1;
        A = A_next;
        B = B_next;
        Chisq_history=[Chisq_history;Chisq_cur];
    else                    % fail, change lam and go to loop
        lam = lam * 10;        
    end
%   Set the break conditions
    if DONE == NDONE
        break;
    end
%   Output the information
    fprintf(['\n Interation:',num2str(iter),', Chisq=',num2str(Chisq_cur),',lam=',num2str(lam)]);
end

ShowFigure = 0;
%%  Show the result
if ShowFigure == 1
    figure;
    plot(log(Chisq_history),'bo');
    title('history of Chisq during iteration');    
end

end