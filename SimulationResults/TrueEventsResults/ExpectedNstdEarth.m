
%% MyData
clear all;  close all;  clc;

data_nue = readmatrix("IntStdEarth/IntStdEarth_4_0prem_44layers10_1001");
data_numu = readmatrix("IntStdEarth/IntStdEarth_4_1prem_44layers10_1001");
%create surface plot
cth = unique(data_nue(:,1)); %Xaxis

e   = unique(data_nue(:,2)); %Yaxis

Nnuei= data_nue(:,3);
Nnumui= data_numu(:,3);

[Cth,E] = meshgrid(cth,e);

Nnue = reshape(Nnuei,length(cth),length(e));
Nnumu = reshape(Nnumui,length(cth),length(e));

figure
OscGramMu=pcolor(Cth,E, Nnumu);
set(gca,'FontSize',20, 'FontName', 'Courier')
set(OscGramMu,'edgecolor','none')
title( 'Perfect resolution N_{\nu_\mu}' ,'FontSize',30);
xlabel('cos(\theta_z)','FontSize',30)
set(gca, 'XDir','reverse')
ylabel('E_{\nu}[GeV]','FontSize',30)
h = colorbar;
ylabel(h,'$$\Delta N/N[\%]$$','Interpreter','latex','FontSize',30)

figure
OscGrame=pcolor(Cth,E, Nnue);
set(gca,'FontSize',20, 'FontName', 'Courier')
set(OscGrame,'edgecolor','none')
title( 'Perfect resolution N_{\nu_e}' ,'FontSize',30);
xlabel('cos(\theta_z)','FontSize',30)
set(gca, 'XDir','reverse')
ylabel('E_{\nu}[GeV]','FontSize',30)
h = colorbar;
ylabel(h,'$$\Delta N/N[\%]$$','Interpreter','latex','FontSize',30)



