close all;
%% DOlivo data

Data2 = DData;

% filename2 = 'events_std_solmax.csv';
% Data2 = readmatrix(filename2)

eta2 = unique(Data2(:,2));
e2   = unique(Data2(:,3));
n_ij2= Data2(:,4);

[Eta2,E2] = meshgrid(eta2,e2);
N2 = reshape(n_ij2,length(eta2),length(e2));
Ntot2=sum(sum(N2))
Abin = 0.1980*0.0198;

figure
s2=surf(Eta2,E2,N2/(Abin*Ntot2));
title('Probability distribution functions of true events: DOlivo');
%s.FaceColor = 'interp';
colorbar



%% My Data
%%Data = ExpData;
filename = 'Sim_data/events_meeting.csv';
Data = readmatrix(filename)


eta = unique(Data(:,1));
e   = unique(Data(:,2));
n_ij= Data(:,3);

[Eta,E] = meshgrid(eta,e);
N = reshape(n_ij,length(eta),length(e));
Ntot = sum(sum(N)) %sum(sum(N));
Abin = 0.1980*0.0198;
figure
s=surf(Eta,E, N);
title('Probability distribution functions of true events: MyData');
%s.FaceColor = 'interp';
colorbar

% %% SOme extra stuff
% 
% figure
% plot(eta, e , '*');
% hold  on
% plot(eta2,e2, 'o');
% hold off


% figure
% [Eta,E] = meshgrid(eta,e);
% dN = N2-N;
% s=mesh(Eta,E,dN );
% colorbar
%% Observed events

clear all 
%%Data = ExpData;
filename = 'Sim_data/ObservedEvents.csv';
Data = readmatrix(filename)


eta = unique(Data(:,1));
e   = unique(Data(:,2));

n_ij= Data(:,3);

[Eta,E] = meshgrid(eta,e);
N = reshape(n_ij,length(eta),length(e));

figure
s=surf(Eta,E, N);
title('Probability distribution functions of true events: MyData');
%s.FaceColor = 'interp';
colorbar

