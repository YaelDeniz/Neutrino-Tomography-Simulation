
%% MyData
clear all;   clc;

file = 'dNTrue0.csv';

TrueData = readmatrix(file)

%create surface plot
cth = unique(TrueData(:,2)); %Xaxis

e   = unique(TrueData(:,3)); %Yaxis

dn_ij= TrueData(:,4);

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

cmb = cosd(180-asind((3480)/6371) )
PileLowerTop = cosd(180-asind((3480+1000/3)/6371) )
PileMidTop = cosd(180-asind((3480+2*1000/3)/6371) )
PileTop = cosd(180-asind((3480+1000)/6371) )

cmbl = xline(cmb,'--','CMB','LineWidth',2.5) %CMB
cmbl.FontSize=25;

topl = xline(PileLowerTop,'--',' Lower Segment','LineWidth',2.5) %Upper Mantle
topl.FontSize=10;

topl = xline(PileMidTop,'--','Middle Segment','LineWidth',2.5) %Upper Mantle
topl.FontSize=10;

topl = xline(PileTop,'--','Top of LLVP','LineWidth',2.5) %Upper Mantle
topl.FontSize=25;

hold off


h = colorbar;
ylabel(h,'$$\Delta N/N[\%]$$','Interpreter','latex','FontSize',30)




