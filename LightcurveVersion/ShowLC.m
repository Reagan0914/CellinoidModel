%% @FileName
%  Description:
%    Plot light curve   
%
%  How To Use:
%   >> ShowLC(jdt, SynLC)
%   
%   Edited by LUXP
%   Date: 2016-09-26

function ShowLC(JDT, SynLC, LCm)
for i=1:size(LCm,1)
    figure,
    set(gcf,'Position',[100,200,725,770]);
    EndInd = sum(LCm(1:i));
    Ind = (EndInd-LCm(i)+1):EndInd;
    X = JDT(Ind)*24;
    Y = SynLC(Ind);
    plot(X, Y,'r*');
    title('Synthetic Lightcurves of Cellinoid By Lebedev','fontsize',16);
    xlabel('Time (Hours)','fontsize',22);
    set(gca,'FontName','Times New Roman','FontSize',16,'fontweight','bold');
    ylabel('Brightness','fontsize',22);
end
end
