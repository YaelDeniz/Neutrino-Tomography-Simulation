%% True events Joao Integration

clear all; clc

truestd = 'AsimovData/prem_default_Asimov_1_100_100_1.000000_10.000000.csv';
nexp_true = readmatrix(truestd);
   
truellsvp = 'AsimovData/prem_llsvp_Asimov_1_100_100_1.000000_10.000000.csv';
nobs_true = readmatrix(truellsvp);
 
data(:,1:2) = nexp_true(:,1:2);
data(:,3)  = 100*( nobs_true(:,3)-nexp_true(:,3) )./nexp_true(:,3); %%Difference in events

Nexp_s = sum(nexp_true(:,3))

Nobs_s = sum(nobs_true(:,3))

%create surface plot
eta_t = unique(data(:,1))*(180.0/pi);
ct_t = cosd(eta_t)
e   = unique(data(:,2));
dp_ij= data(:,3);

[E,Etat] = meshgrid(e,eta_t);

[E,Ct_t] = meshgrid(e,ct_t);

dNt = reshape(dp_ij,length(e),length(ct_t));


%Plot events
%%figure 
%%imagesc(n_ij)
figure('Renderer', 'painters', 'Position', [10 10 1000 800])
%datatrue=pcolor(Et,Etat, dNt');
datatrue=pcolor(Ct_t,E, dNt');
set(gca,'FontSize',20, 'FontName', 'Courier')
set(datatrue,'edgecolor','none')
title( 'Perfect resolution' ,'FontSize',30);

xlabel('cos\theta_z','FontSize',30)
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





%% Observed events

clear all; clc

truestd = 'AsimovData/prem_default_AsiObsNextGen_1_100_20_1.000000_10.000000.csv';
nexp_true = readmatrix(truestd);
%   
truellsvp = 'AsimovData/prem_llsvp_AsiObsNextGen_1_100_20_1.000000_10.000000.csv';
nobs_true = readmatrix(truellsvp);
% 
data(:,1:2) = nexp_true(:,1:2);
data(:,3)  = 100*( nobs_true(:,3)-nexp_true(:,3) )./nexp_true(:,3); %%Difference in events
%data(:,3)  =nexp_true(:,3)  ; %%Difference in events
%Total neutrinos

Nexp_s = sum(nexp_true(:,3))

Nobs_s = sum(nobs_true(:,3))


%data(:,3) = nexp_true(:,3) ; %%Difference in events

%create surface plot
eta_t = unique(data(:,1))*(180.0/pi);
ct_t = cosd(eta_t)
e   = unique(data(:,2));
dp_ij= data(:,3);

[E,Etat] = meshgrid(e,eta_t);

[E,Ct_t] = meshgrid(e,ct_t);

dNt = reshape(dp_ij,length(e),length(ct_t));

%dNt = reshape(dp_ij,length(e),length(eta_t));


%Plot events
%%figure 
%%imagesc(n_ij)
figure('Renderer', 'painters', 'Position', [10 10 1000 800])

datatrue=pcolor(Ct_t,E, dNt');


set(gca,'FontSize',20, 'FontName', 'Courier')

set(datatrue,'edgecolor','none')

title('Realistic resolution','FontSize',30);

xlabel('cos\theta^{reco}_z','FontSize',30)
set(gca, 'XDir','reverse')

ylabel('E^{reco}_{\nu}[GeV]','FontSize',30)


hold on 

topl = xline(cosd(138.8),'--','Top of LLSVP','LineWidth',2.5) %Upper Mantle
topl.FontSize=25;

cmbl = xline(cosd(146.8),'--','CMB','LineWidth',2.5) %CMB
cmbl.FontSize=25;

hold off


h = colorbar;
ylabel(h,'$$\Delta N/N[\%]$$','Interpreter','latex','FontSize',30)

%% OSCILLATION PROBABILITIES:

Prob_std = 'OscProbEarth/prem_default_Probemu_1_200_200_1.000000_10.000000.csv';
pij_std= readmatrix(Prob_std);
   
Prob_alt = 'OscProbEarth/prem_llsvp_Probemu_1_200_200_1.000000_10.000000.csv';
pij_alt = readmatrix(Prob_alt);
 
data(:,1:2) = pij_alt(:,1:2);

data(:,3)  = ( pij_std(:,3) - pij_alt(:,3) ); %% Difference in oscillation probabilities

%data(:,3)  = pij_std(:,3); %% Difference in oscillation probabilities



%create surface plot
theta = unique(data(:,1))*(180.0/pi);
e   = unique(data(:,2));
dp_ij= data(:,3);

[E,Theta] = meshgrid(e,theta);

dPij = reshape(dp_ij,length(e),length(theta));


%Plot Probabilities

figure('Renderer', 'painters', 'Position', [10 10 1000 800])

OscProb = pcolor(Theta,E, dPij');

set(gca,'FontSize',20, 'FontName', 'Courier')

set(OscProb,'edgecolor','none')

title( 'Oscillation Probabilities' ,'FontSize',30);

xlabel('\theta_z[°]','FontSize',30)

%set(gca, 'XDir','reverse')

ylabel('E_{\nu}[GeV]','FontSize',30)


hold on 

topl = xline(138.8,'--','Top of LLVP','LineWidth',2.5) %Upper Mantle
topl.FontSize=25;

cmbl = xline(146.8,'--','CMB','LineWidth',2.5) %CMB
cmbl.FontSize=25;

hold off


h = colorbar;
ylabel(h,'$$\Delta P/P[\%]$$','Interpreter','latex','FontSize',30)

%set(gca, 'XDir','reverse')


%% My Integration


% truestd = 'PoisMeans/prem_default_PoiMeansT_1_200_200_1.000000_10.000000.csv';
% nexp_true = readmatrix(truestd);
% %   
% truellsvp = 'PoisMeans/prem_llsvp_PoiMeansT_1_200_200_1.000000_10.000000.csv';
% nobs_true = readmatrix(truellsvp);
% % 
% data(:,1:2) = nobs_true(:,1:2);
% data(:,3)  = 100*abs(( nobs_true(:,3)-nexp_true(:,3) ))./nexp_true(:,3); %%Difference in events
% %data(:,3)  =nexp_true(:,3)  ; %%Difference in events
% %Total neutrinos
% 
% Nexp_s = sum(nexp_true(:,3))
% 
% Nobs_s = sum(nobs_true(:,3))
% 
% 
% %data(:,3) = nexp_true(:,3) ; %%Difference in events
% 
% %create surface plot
% eta_t = unique(data(:,1));
% e_t   = unique(data(:,2));
% dnt_ij= data(:,3);
% 
% [Et,Etat] = meshgrid(e_t,eta_t);
% 
% dNt = reshape(dnt_ij,length(e_t),length(eta_t));
% 
% 
% %Plot events
% %%figure 
% %%imagesc(n_ij)
% figure('Renderer', 'painters', 'Position', [10 10 1000 800])
% datatrue=pcolor(Et,Etat, dNt');
% set(gca,'FontSize',13, 'FontName', 'Courier')
% set(datatrue,'edgecolor','none')
% title('Percentage difference \mu-like True events: Quadrature','  \Delta \rho / \rho = +2% ','FontSize',20);
% xlabel('E^o_{\nu}[GeV]','FontSize',25)
% ylabel('\eta^o','FontSize',25)
% 
% 
% h = colorbar;
% ylabel(h,'$$\frac{|N^{llsvp,o}_{ij} - N^{std,o}_{ij}|}{N^{std,o}_{ij}}\times100\%$$','Interpreter','latex','FontSize',20)


% %% Monte Carlo Simulation result:
% 
% mceventst = 'TrueEventsResults/Mydata.csv';
% 
% mc_data_true = readmatrix(mceventst);
%  
% %create surface plot
% eta = unique(mc_data_true(:,1));
% e   = unique(mc_data_true(:,2));
% %dn_ij=abs(mc_data_true(:,5));
% 
% 
% 
% dnij = mc_data_true(:,5);
% dnij(isinf(dnij)|isnan(dnij)) = 0;
% 
% %dn_ij= 100*abs(mc_data_true(:,3)-mc_data_true(:,4))./mc_data_true(:,4) ;
% dn_ij=dnij;
% 
% Nexp = sum(mc_data_true(:,3))
% 
% Nobs = sum(mc_data_true(:,4))
% 
% %dn_ij=mc_data_true(:,3)/1000;
% 
% [E,Eta] = meshgrid(e,eta);
% 
% dN = reshape(dn_ij,length(e),length(eta));
% 
% 
% 
% 
% 
% % htrue=pcolor(E,Eta, dN');
% figure('Renderer', 'painters', 'Position', [10 10 1000 800])
% htrue=pcolor(E,Eta, dN');
% set(gca,'FontSize',13, 'FontName', 'Courier')
% set(htrue,'edgecolor','none')
% 
% title('Percentage difference in \mu-like true events','\Delta \rho / \rho = +2','FontSize',20);
% xlabel('E_{\nu}[GeV]','FontSize',25)
% ylabel('\eta','FontSize',25)
% 
% h1 = colorbar;
% ylabel(h1,'$$\frac{N^{llsvp}_{ij} - N^{std}_{ij} }{N^{std}_{ij}}\times100\%$$','Interpreter','latex','FontSize',20)
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % %% True Eventsc 
% % 
% % %%TrueEventStd = 'Events_True_LLSVPlow.csv';
% % %%TrueData = readmatrix(TrueEventStd);
% % 
% % close all; clear all; clc;
% % Events = 'Ntrue1_101101.csv';
% % dolivo_data = readmatrix(Events);
% %  
% % 
% % TrueEventLLSVP = 'NtrueAlt1_101101.csv';
% % TrueData_llsvp = readmatrix(TrueEventLLSVP);
% % 
% % refdata(:,1:2) = TrueData_llsvp(:,1:2);
% % refdata(:,3) =  100*(TrueData_llsvp(:,3) - dolivo_data(:,3))./dolivo_data(:,3) ; %%Difference in events
% % 
% % 
% % %create surface plot
% % eta = unique(refdata(:,1));
% % e   = unique(refdata(:,2));
% % n_ij= refdata(:,3);
% % 
% % [E,Eta] = meshgrid(e,eta);
% % 
% % N = reshape(n_ij,length(e),length(eta));
% % 
% %  
% % 
% % N_true_tot = sum(dolivo_data(:,3)) %sum(sum(N));
% % 
% % Max = max(dolivo_data(:,3))
% % Min = min(dolivo_data(:,3))
% % 
% % 
% % %Plot events
% % %%figure 
% % %%imagesc(n_ij)
% % figure('Renderer', 'painters', 'Position', [10 10 1000 800])
% % strue=pcolor(E,Eta, N');
% % set(gca,'FontSize',13, 'FontName', 'Courier')
% % set(strue,'edgecolor','none')
% % title('\mu-like events:\Delta \rho / \rho = +1% ','FontSize',20);
% % xlabel('E^{\nu}[GeV]','FontSize',25)
% % ylabel('\eta°','FontSize',25)
% % 
% % 
% % % ttext = text(4.5,37.75,'$$\nu_\mu + \bar{\nu}_\mu $$','Interpreter','latex','FontSize',15,'FontWeight','bold')
% % % ttext(1).Color = 'white';
% % % ttext = text(4.4,37.2,'$$\frac{\Delta\rho}{\rho}=-1\%$$','Interpreter','latex','FontSize',15,'FontWeight','bold')
% % % ttext(1).Color = 'white';
% % 
% % h1 = colorbar;
% % ylabel(h1,'$$\frac{N^{llsvp}_{ij} - N^{std}_{ij} }{N^{std}_{ij}}\times100\%$$','Interpreter','latex','FontSize',20)
% % 
% % saveas(gcf,'dEvents_True_high.png')
% % 
% % %% Observed events
% % 
% % 
% % %ObsEvents = 'Events_Obs_LLSVPlow.csv';
% % %ObsData = readmatrix(ObsEvents);
% % 
% % 
% % 
% % ObsEventStd = 'Events_Obs_Std.csv';
% % ObsData_std = readmatrix(ObsEventStd);
% %  
% % 
% % ObsEventLLSVP = 'Events_Obs_LLSVPlow.csv';
% % ObsData_llsvp = readmatrix(ObsEventLLSVP);
% % 
% % ObsData(:,1:2) = ObsData_std(:,1:2);
% % ObsData(:,3) =  100*(ObsData_llsvp(:,3)-ObsData_std(:,3))./ObsData_std(:,3) ; %%Difference in events
% % 
% % 
% % 
% % 
% % 
% % %create surface plot
% % eta_o = unique(ObsData(:,1));
% % e_o= unique(ObsData(:,2));
% % n_mn= ObsData(:,3);
% % 
% % [E_o,Eta_o] = meshgrid(e_o,eta_o);
% % 
% % No = reshape(n_mn,length(e_o),length(eta_o));
% % 
% % N_obs_tot = 100*sum(sum(No))
% % 
% % %Plot observed events
% % 
% % figure('Renderer', 'painters', 'Position', [10 10 1000 800])
% % sobs=pcolor(E_o,Eta_o, No');
% % set(sobs,'edgecolor','none')
% % set(gca,'FontSize',13, 'FontName', 'Courier')
% % 
% % title('\mu-like events:\Delta \rho / \rho = -1% ','FontSize',20);
% % xlabel('E^{\nu}_o[GeV]','FontSize',25)
% % ylabel('\eta°_o','FontSize',25)
% % 
% % 
% % % otext1 = text(4.5,37.75,'$$\nu_\mu + \bar{\nu}_\mu $$','Interpreter','latex','FontSize',15,'FontWeight','bold')
% % % otext1(1).Color = 'white';
% % % otext2 = text(4.4,37.2,'$$\frac{\Delta\rho}{\rho}=-1\%$$','Interpreter','latex','FontSize',15,'FontWeight','bold')
% % % otext2(1).Color = 'white';
% % 
% % %ttext2 = text(4.75,37.35,'$$\sigma_\eta = \frac{0.25}{\sqrt{\frac{E^\nu}{GeV}}} $$','Interpreter','latex','FontSize',10)
% % %ttext3 = text(4.75,37.1,'$$\sigma_E^\nu = 0.20 \times E^\nu $$','Interpreter','latex','FontSize',10)
% % % ttext1(1).Color = 'white';
% % %ttext2(1).Color = 'white';
% % %ttext3(1).Color = 'white';
% % 
% % h1 = colorbar
% % ylabel(h1,'$$\frac{N^{llsvp}_{mn} - N^{std}_{mn} }{N^{std}_{mn}}\times100\%$$','Interpreter','latex','FontSize',18)
% % 
% % saveas(gcf,'dEvents_Obs_low.png')
