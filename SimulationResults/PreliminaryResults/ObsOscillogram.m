%Interacting Events Oscillogram

%% Muon-like events
clear all;  close all;  clc;
ObsData_ithAzi_mu = 'ObsEvents/ObsLLVP_pancakenu1_1010010_50.txt';
obsDataMu = readmatrix(ObsData_ithAzi_mu)

%create surface plot
cth = unique(obsDataMu(:,1)); %Xaxis
e   = unique(obsDataMu(:,2)); %Yaxis
nmu_ij= obsDataMu(:,3);
[Cth,E] = meshgrid(cth,e);
Nmu = reshape(nmu_ij,length(cth),length(e));

% % electron-like events
% ObsData_ithAzi_e = 'IntEvents/IntLLVP_pancakenu0_5010050_50.txt';
% ObsDatae = readmatrix(ObsData_ithAzi_e)
% 
% %create surface plot
% cth = unique(ObsDatae(:,1)); %Xaxis
% e   = unique(ObsDatae(:,2)); %Yaxis
% ne_ij= ObsDatae(:,3);
% [Cth,E] = meshgrid(cth,e);
% Ne = reshape(ne_ij,length(cth),length(e));

%% Oscillogram visualization

% Muon-like
% figure 
figure('Renderer', 'painters', 'Position', [10 10 1000 800])
mudatatrue=pcolor(Cth,E,Nmu);
set(gca,'FontSize',20, 'FontName', 'Courier')
set(mudatatrue,'edgecolor','none')
title( '\nu_\mu - like events[Obs]' ,'FontSize',30);

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
ylabel(h,'$$\Delta N/N[\%]$$','Interpreter','latex','FontSize',30)

% % electron-like
% % figure 
% figure('Renderer', 'painters', 'Position', [10 10 1000 800])
% mudatatrue=pcolor(Cth,E,Ne);
% set(gca,'FontSize',20, 'FontName', 'Courier')
% set(mudatatrue,'edgecolor','none')
% title( '\nu_e - like events[Obs]' ,'FontSize',30);
% 
% xlabel('cos(\theta_z)','FontSize',30)
% set(gca, 'XDir','reverse')
% ylabel('E_{\nu}[GeV]','FontSize',30)
% 
% hold on 
% topl = xline(cosd(138.8),'--','Top of LLVP','LineWidth',2.5) %Upper Mantle
% topl.FontSize=25;
% cmbl = xline(cosd(146.8),'--','CMB','LineWidth',2.5) %CMB
% cmbl.FontSize=25;
% hold off
% 
% h = colorbar;
% ylabel(h,'$$\Delta N/N[\%]$$','Interpreter','latex','FontSize',30)