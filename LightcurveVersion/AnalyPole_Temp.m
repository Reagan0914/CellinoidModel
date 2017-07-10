%%  @AnalyPole.m
%
%   Descriptions:
%       Analyze the pole search result and get the best-fit pole
%       orientation
%
%   Edited by LUXP
%   Date: 2016-10-07
function AnalyPole_Temp

%% specify the directory for all pole_result mat files
clear;
Folder = './OUTPUT/';
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
ChisqBenchNum = 100;
FigPole(ChisqBenchNum, TotalCg, Chisq_History);
%%  Best Fit Period
%   run: >N=find(min(TotalChisq)==TotalChisq);TotalCg(9:10,N);
end

function FigPole(ChisqBenchNum, a_History, Chisq_History)



ClassNum = 4;
%% kick off the odd point:
%Chisq_History(Chisq_History == min(Chisq_History)) = 2;

%%  extract (lam,beta) from a_History and show them
fitlam = a_History(7,:)';
fitbeta = a_History(8,:)';

figure,
set(gcf,'Position',[30,500,850,550]);
plot(fitlam, fitbeta,'MarkerFaceColor',[1 1 1],...
    'MarkerEdgeColor',[0.800000011920929 0.800000011920929 0.800000011920929],...
    'MarkerSize',2,...
    'Marker','*',...
    'LineStyle','none',...
    'Color',[0 0 0] );
axis equal;
axis([0, 360, -90, 90]);
set(gca,'YTick',[-90 -60 -30 0 30 60 90],...
    'XTick',[0 30 60 90 120 150 180 210 240 270 300 330 360]);
%set(gca,'FontName','Times New Roman','FontSize',20);
set(gca,'FontName','Helvetica','FontSize',20);
xlabel('Ecliptic Longitude of Fitted Pole','fontsize',20);
ylabel('Ecliptic Latitude of Fitted Pole','fontsize',20);
%% only show one best fit 
if ChisqBenchNum == 1
    BestChi = min(Chisq_History);
    bestcg = a_History(:, Chisq_History==BestChi)
    bestlam = bestcg(7); bestbeta = bestcg(8);
    hold on;plot(bestlam,bestbeta, 'k', 'MarkerSize',10,'Marker','d','LineWidth',2,'LineStyle','none');
    return 
end

%%  Analyze the best-fit results with the smallest Chisq
ChisqBenchMark = sort(Chisq_History);
ChisqBenchMark = ChisqBenchMark(ChisqBenchNum);
ind_good = (Chisq_History < ChisqBenchMark);
fitlam_good = fitlam(ind_good);
fitbeta_good = fitbeta(ind_good);
fitChi_good = Chisq_History(ind_good);
FitPole = [fitlam_good, fitbeta_good, fitChi_good];

%%%%%   FOR TEMP GENERATION %%%%%%%%%
if 1
    bestchi = min(Chisq_History)
    bestlam = fitlam(Chisq_History == bestchi)
    bestbeta = fitbeta(Chisq_History == bestchi)
    reverse_lam = FitPole(FitPole(:,1)<50,1);
    reverse_beta = FitPole(FitPole(:,1)<50,2);
    hold on;
    plot(FitPole(:,1), FitPole(:,2), 'k', 'MarkerSize',10,'Marker','d','LineWidth',1.5,'LineStyle','none');
    plot(bestlam,bestbeta, 'r', 'MarkerSize',15,'Marker','^','LineWidth',2,'LineStyle','none');
    plot(reverse_lam,reverse_beta, 'r', 'MarkerSize',15,'Marker','v','LineWidth',2,'LineStyle','none');
    return;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%   Classify good fit results into two class
[ClassFit, ClassCenter] = kmeans([fitlam_good, fitbeta_good], ClassNum);
%modify the sort, bad matlab function!
[Y,II]= sort(ClassCenter);
ClassFit = ClassFit + 5;
ClassFit(ClassFit==6) = II(1,1);
ClassFit(ClassFit==7) = II(2,1);
ClassFit(ClassFit==8) = II(3,1);
ClassFit(ClassFit==9) = II(4,1);
% FitLam = zeros(size(fitlam_good,1), ClassNum/2);
% FitBeta = zeros(size(fitlam_good,1), ClassNum/2);
% FitChi = zeros(size(fitlam_good,1), ClassNum/2);
for classInd = 1:ClassNum
    FitLamT = fitlam_good(ClassFit == classInd);
    FitBetaT = fitbeta_good(ClassFit == classInd);
    FitChiT = fitChi_good(ClassFit == classInd);
    if classInd == 1
        FitPole1 = [FitLamT, FitBetaT, FitChiT]; 
    elseif classInd == 2
        FitPole2 = [FitLamT, FitBetaT, FitChiT];
    elseif classInd == 3
        FitPole3 = [FitLamT, FitBetaT, FitChiT];
    elseif classInd == 4
        FitPole4 = [FitLamT, FitBetaT, FitChiT];
    end
end
% fitlam1 = fitlam_good(ClassFit == 1);
% fitlam2 = fitlam_good(ClassFit == 2);
% 
% fitbeta1 = fitbeta_good(ClassFit == 1);
% fitbeta2 = fitbeta_good(ClassFit == 2);
% 
% fitChi1 = fitChi_good(ClassFit == 1);
% fitChi2 = fitChi_good(ClassFit == 2);

bestChi = min(Chisq_History);
bestlam = fitlam(Chisq_History == bestChi);
bestbeta = fitbeta(Chisq_History == bestChi);
bestcg = a_History(:, Chisq_History == bestChi);

disp('Best-fit results:');
FitPoleBest = [bestcg; bestChi]


%%  Mark on the figure
hold on;
%plot(bestlam, bestbeta, 'ko', 'MarkerSize',10,'Marker','o','LineWidth',3,'LineStyle','none');
plot(FitPole1(:,1), FitPole1(:,2), 'k', 'MarkerSize',10,'Marker','d','LineWidth',2,'LineStyle','none');
plot(FitPole2(:,1), FitPole2(:,2), 'k', 'MarkerSize',10,'Marker','v','LineWidth',2,'LineStyle','none');

plot(FitPole3(:,1), FitPole3(:,2), 'k', 'MarkerSize',10,'Marker','o','LineWidth',2,'LineStyle','none');
plot(FitPole4(:,1), FitPole4(:,2), 'k', 'MarkerSize',10,'Marker','^','LineWidth',2,'LineStyle','none');

disp('Num of best-fit results in each pole group from 1: 4');
[size(FitPole1,1), size(FitPole2,1), size(FitPole3,1), size(FitPole4,1)]

%%  Show Result
for i= 1:ClassNum
    if i == 1
        FitPoleT = FitPole1;
        disp('FitPole1: min,mean,max values for (lam, beta, chi)');
    elseif i==2
        FitPoleT = FitPole2;
        disp('FitPole2: min,mean,max values for (lam, beta, chi)');
    elseif i==3
        FitPoleT = FitPole3;
        disp('FitPole3: min,mean,max values for (lam, beta, chi)');
     elseif i==4
        FitPoleT = FitPole4;
        disp('FitPole4: min,mean,max values for (lam, beta, chi)');
    end
    [min(FitPoleT(:,1)), mean(FitPoleT(:,1)), max(FitPoleT(:,1)), min(FitPoleT(:,2)), mean(FitPoleT(:,2)), max(FitPoleT(:,2)), min(FitPoleT(:,3)), mean(FitPoleT(:,3)), max(FitPoleT(:,3))]
end
% disp('FitPole1: min,mean,max values for (lam, beta, chi)');
% PoleCMP1 = [min(fitlam1),mean(fitlam1),max(fitlam1), min(fitbeta1),mean(fitbeta1), max(fitbeta1), min(fitChi1), mean(fitChi1),max(fitChi1)]
% disp('FitPole2: max and min values for (lam, beta, chi) respectively');
% PoleCMP2 = [min(fitlam2),mean(fitlam2),max(fitlam2), min(fitbeta2), mean(fitbeta2),max(fitbeta2) ,min(fitChi2),mean(fitChi2), max(fitChi2)]
a_good = a_History(:, ind_good);
a_good1 = a_good(:, ClassFit==1);
a_good2 = a_good(:, ClassFit==2);

a_good3=a_good1;
for i=1:size(a_good3,2)
   a_good3(1:6,i) = a_good3(1:6,i)/max(a_good3(1:6,i)); 
   a_good3(1:6,i) = sort((a_good3(1:6,i)),'descend');
end

a_good4=a_good2;
for i=1:size(a_good4,2)
   a_good4(1:6,i) = a_good4(1:6,i)/max(a_good4(1:6,i)); 
   a_good4(1:6,i) = sort((a_good4(1:6,i)),'descend');
end

disp('Semi-Axes Shows:');
SemiAxes3 = [min(a_good3(1:6,:),[],2), mean(a_good3(1:6,:),2), max(a_good3(1:6,:),[],2)]
SemiAxes4 = [min(a_good4(1:6,:),[],2), mean(a_good4(1:6,:),2), max(a_good4(1:6,:),[],2)]
%%  Show information
%   FitPole: show the [lam, beta, Chi] with the 100 smallest chisquare
%   FitPole1, FitPole2, two classified FitPole
%   PoleCMP1,2: Show min,mean,max of lam,beta,chi
%
%   a_good1, a_good2: two classified fitted cg
%   a_good3, a_good4: two classified fitted cg with ajusted semi-axes
%   SemiAxes3,4: two classified fitted semiaxes, min, mean, max of six semi
%   FitPoleBest: best cg

if 1
a_goodFit = a_History(:, ind_good);
a_goodFit1 = a_goodFit(:, ClassFit==1);
a_goodFit2 = a_goodFit(:, ClassFit==2);

a_goodFit3=a_goodFit1;
for i=1:size(a_goodFit3,2)
   a_goodFit3(1:6,i) = a_goodFit3(1:6,i)/max(a_goodFit3(1:6,i)); 
   %a_good3(1:6,i) = sort((a_good3(1:6,i)),'descend');
end

a_goodFit4=a_goodFit2;
for i=1:size(a_goodFit4,2)
   a_goodFit4(1:6,i) = a_goodFit4(1:6,i)/max(a_goodFit4(1:6,i)); 
   %a_good4(1:6,i) = sort((a_good4(1:6,i)),'descend');
end
format long
disp('Best-fit period: [min, mean, max] of Fit pole 1');
[min(a_goodFit1(9,:)), mean(a_goodFit1(9,:)), max(a_goodFit1(9,:))]
disp('Best-fit period: [min, mean, max] of Fit pole 2');
[min(a_goodFit2(9,:)), mean(a_goodFit2(9,:)), max(a_goodFit2(9,:))]
format short
end
end