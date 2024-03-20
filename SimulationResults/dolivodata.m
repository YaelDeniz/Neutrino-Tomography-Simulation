
dolivo_data = refdata;
 
refdata(:,1:2) = dolivo_data(:,1:2);
refdata(:,3)   = dolivo_data(:,3); %%Difference in events

%create surface plot
eta = unique(refdata(:,1));
e   = unique(refdata(:,2));
n_ij= refdata(:,3);

[E,Eta] = meshgrid(e,eta);

Ndolivo = reshape(n_ij,length(e),length(eta));

 
N_dolivo_tot = sum(dolivo_data(:,3)) %sum(sum(N));

Max = max(dolivo_data(:,3))
Min = min(dolivo_data(:,3))


%Plot events
%%figure 
%%imagesc(n_ij)
%figure('Renderer', 'painters', 'Position', [10 10 1000 800])
%strue=pcolor(E,Eta, N');
figure
sdolivo=mesh(E,Eta, Ndolivo');
%set(gca,'FontSize',13, 'FontName', 'Courier')
%set(strue,'edgecolor','none')
title('\mu-like events:\Delta \rho / \rho = +1 | MC simulation% ','FontSize',20);
xlabel('E^{\nu}[GeV]','FontSize',25)
ylabel('\eta°','FontSize',25)


% ttext = text(4.5,37.75,'$$\nu_\mu + \bar{\nu}_\mu $$','Interpreter','latex','FontSize',15,'FontWeight','bold')
% ttext(1).Color = 'white';
% ttext = text(4.4,37.2,'$$\frac{\Delta\rho}{\rho}=-1\%$$','Interpreter','latex','FontSize',15,'FontWeight','bold')
% ttext(1).Color = 'white';

h1 = colorbar;
ylabel(h1,'$$\frac{N^{llsvp}_{ij} - N^{std}_{ij} }{N^{std}_{ij}}\times100\%$$','Interpreter','latex','FontSize',20)