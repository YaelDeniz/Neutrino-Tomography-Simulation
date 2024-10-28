%Interacting Events Oscillogram

%% Muon-like events
clear all;  close all;  clc;
IntData_ithAzi_mu = 'IntEvents/IntcthLLVP_pancakenu1_5010050_50.txt';
IntDataMu = readmatrix(IntData_ithAzi_mu)

%create surface plot
cth = unique(IntDataMu(:,1)); %Xaxis
e   = unique(IntDataMu(:,2)); %Yaxis
numu_std= IntDataMu(:,3);
numu_alt= IntDataMu(:,4);
numu_ij= IntDataMu(:,5);
%numu_ij= IntDataMu(:,5);
[Cth,E] = meshgrid(cth,e);
Nnumu = reshape(numu_ij,length(cth),length(e));

% electron-like events
IntData_ithAzi_e = 'IntEvents/IntcthLLVP_pancakenu0_5010050_50.txt';
IntDatae = readmatrix(IntData_ithAzi_e)

%create surface plot
cth = unique(IntDatae(:,1)); %Xaxis
e   = unique(IntDatae(:,2)); %Yaxis
nue_ij= IntDatae(:,5);
[Cth,E] = meshgrid(cth,e);
Nnue = reshape(nue_ij,length(cth),length(e));

%% Oscillogram visualization

% Muon-like
% figure 
figure('Renderer', 'painters', 'Position', [10 10 1000 800])
mudatatrue=pcolor(Cth,E,Nnumu);
set(gca,'FontSize',20, 'FontName', 'Courier')
set(mudatatrue,'edgecolor','none')
title( '\nu_\mu detections [Int]' ,'FontSize',30);

xlabel('cos(\theta_z)','FontSize',30)
set(gca, 'XDir','reverse')
ylabel('E_{\nu}[GeV]','FontSize',30)


hold on 
topl = xline(cosd(138.8),'--','Top of LLVP','LineWidth',2.5) %Upper Mantle
topl.FontSize=25;
cmbl = xline(cosd(146.8),'--','CMB','LineWidth',2.5) %CMB
cmbl.FontSize=25;
hold off

h = colorbar;
caxis([-6 12])

ylabel(h,'$$\Delta N/N[\%]$$','Interpreter','latex','FontSize',30)
saveas(gcf,'NumuA0.png')


% electron-like
% figure 
figure('Renderer', 'painters', 'Position', [10 10 1000 800])
mudatatrue=pcolor(Cth,E,Nnue);
set(gca,'FontSize',20, 'FontName', 'Courier')
set(mudatatrue,'edgecolor','none')
title( '\nu_e - like events[Int]' ,'FontSize',30);

xlabel('cos(\theta_z)','FontSize',30)
set(gca, 'XDir','reverse')
ylabel('E_{\nu}[GeV]','FontSize',30)

hold on 
topl = xline(cosd(138.8),'--','Top of LLVP','LineWidth',2.5) %Upper Mantle
topl.FontSize=25;
cmbl = xline(cosd(146.8),'--','CMB','LineWidth',2.5) %CMB
cmbl.FontSize=25;
hold off

h = colorbar;
caxis([-1.5 1])

ylabel(h,'$$\Delta N/N[\%]$$','Interpreter','latex','FontSize',30)
saveas(gcf,'NueA0.png')