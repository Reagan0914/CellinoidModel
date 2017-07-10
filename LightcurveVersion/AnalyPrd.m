%%  @AnalyPrd.m
%
%   Descriptions:
%       Analyze the period search result and get the best-fit period
%
%   Edited by LUXP
%   Date: 2016-10-06
function AnalyPrd

%% specify the directory for all prd_result mat files
clear;
Folder = './PrdSearch_FixPole_Ellipsoid_Refine/';
FileNames = dir(Folder);
TotalCg = [];
TotalChisq = [];
for i = 1:size(FileNames,1)
    FileName = FileNames(i).name;
    if length(FileName)>4 && strcmp(FileName(end-3:end), '.mat')
        fprintf('Loading Data File: %s\n', FileName);
        FileName = [Folder, FileName];
        load(FileName);
        TotalCg = [TotalCg cg_History];
        TotalChisq = [TotalChisq; Chisq_History];
        clear cg_History, Chisq_History;
    end
end

%% show the plot
figure,
semilogy(TotalCg(9,:),TotalChisq,'k*');
axis([min(TotalCg(9,:)), max(TotalCg(9,:)),0.9*min(TotalChisq),1.1*max(TotalChisq)]);
set(gca,'FontName','Times New Roman','FontSize',20,'fontweight','bold');
xlabel('Period (Hour)','fontsize',20);
ylabel('\chi^2','fontsize',20);
set(gcf,'Position',[800,500,700,640]);
%   add the marker for the smallest chisq
fitchi = min(TotalChisq); bestprd = TotalCg(9,TotalChisq==fitchi);
hold on; plot(bestprd, fitchi, 'o','MarkerSize',15);
bestprd
%%  Best Fit Period
%   run: >N=find(min(TotalChisq)==TotalChisq);TotalCg(9:10,N);
end