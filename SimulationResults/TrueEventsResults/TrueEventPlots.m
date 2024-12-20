
%% MyData
clear all;  close all;  clc;

fileyael = 'EarthSensitivity2D.csv';

fileapc = 'EarthSensitivity2D_oscprob.csv';

event_yael = readmatrix(fileyael)

event_apc = readmatrix(fileapc)

data(:,1:2) = event_yael(:,1:2);

data(:,3)  = event_yael(:,3)-event_apc(:,3); %%Difference in events

%create surface plot
cth = unique(data(:,1)); %Xaxis

e   = unique(data(:,2)); %Yaxis

dn_ij= data(:,3);

[Cth,E] = meshgrid(cth,e);


dN = reshape(dn_ij,length(cth),length(e));


%Plot events
%%figure 
%%imagesc(n_ij)
figure('Renderer', 'painters', 'Position', [10 10 1000 800])
%datatrue=pcolor(Et,Etat, dNt');
datatrue=pcolor(Cth,E, dN);
set(gca,'FontSize',20, 'FontName', 'Courier')
set(datatrue,'edgecolor','none')
title( 'Perfect resolution' ,'FontSize',30);

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




