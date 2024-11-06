
%% Data
addpath VisToolBox
clear all
close all

% Setting data

Npdiff= [];
AzimuthTable = readmatrix("IntEvents/TrueIndexTable.txt");
phi=AzimuthTable(1:101,2)+50;



resultPath = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/SimulationResults/PreliminaryResults/IntEvents/"
SimulationPath = "Events_Simulation_pancake_100Mton_1-20GeV_128-171Zen_-50-50Az/"

for j = 0:1:100;

filename   = "IntcthLLVP_nu1_10Zen100Az10Enu_"+int2str(j)+".csv";
eventFile  = resultPath+SimulationPath+filename;

Osc4Azi = readmatrix(eventFile);

cth      = unique(Osc4Azi(:,1)); %Xaxis
piminth  = 180 - acosd(cth); 
ene   = unique(Osc4Azi(:,2)); %Yaxis
nstd  = Osc4Azi(:,3);
nalt  = Osc4Azi(:,4);


[Cth,Ene] = meshgrid(cth,ene);

Nstd  = reshape(nstd,length(cth),length(ene));
Nalt  = reshape(nalt,length(cth),length(ene));

Npdiff_azi= sum(Nstd-Nalt)./sum(Nstd); 

Npdiff(end + 1, :) = [Npdiff_azi']
%cth = unique(Set_cthphi(:,1)) 


% for i = 1:length(cth)
% 
%    Nstd_thphi = sum(Set_cthphi(Set_cthphi(:,1) == cth(i),3));
%    Nalt_thphi = sum(Set_cthphi(Set_cthphi(:,1) == cth(i),4));
%    dN = ((Nstd_thphi-Nalt_thphi)/Nstd_thphi)*100;
%    N(end+1) = [dN];       
% 
% end   



end

%%

% [Phi,Cth] = meshgrid(phi,cth);
% Nmu = reshape(N,length(phi),length(cth));
% [~,~,Z] = peaks(100);
cthabs = -1.0.*cth;

[~,c] = polarPcolor(piminth',phi',Npdiff,'autoOrigin','on'); % colormap is here an optional argument
ylabel(c,"$\frac{\Delta N^{(1GeV\leq E_\nu \leq 20GeV)}}{N^{(1GeV\leq E_\nu \leq 20GeV)}}(\phi,\theta_Z)$",'FontSize', 20,Interpreter="latex");
t = title('this is  my title', 'Units', 'normalized', 'Position', [0.5, 0.75, 0]);
t.Color = 'r'; t.FontSize = 10; % with this you can change color, font name and size


