
%% MyData
clear all;  close all;  clc;

file = 'Testdata.csv';

event_diff = readmatrix(file)

data(:,1:2) = event_diff(:,1:2);

data(:,3)  = event_diff(:,5); %%Difference in events

%create surface plot
zen = unique(data(:,1)); %Xaxis

cz = cosd(zen)

e   = unique(data(:,2)); %Yaxis

dn_ij= data(:,3);

[Zen,E] = meshgrid(zen,e);

%[E,CZ] = meshgrid(e,cz);

dN = reshape(dn_ij,length(zen),length(e));


%Plot events
%%figure 
%%imagesc(n_ij)
figure('Renderer', 'painters', 'Position', [10 10 1000 800])
%datatrue=pcolor(Et,Etat, dNt');
datatrue=pcolor(cosd(Zen),E, dN);
set(gca,'FontSize',20, 'FontName', 'Courier')
set(datatrue,'edgecolor','none')
title( 'Perfect resolution' ,'FontSize',30);

xlabel('theta_z','FontSize',30)
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




