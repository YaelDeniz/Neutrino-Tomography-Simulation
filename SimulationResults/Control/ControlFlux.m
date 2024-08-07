%% Load data
clear all;  close all;  clc;

controlfile = 'ControlFlux/ControlFlux_std_.csv';

controlmatrix = readmatrix(controlfile)

%% SetData

cth = unique(controlmatrix(:,1));
ene = unique(controlmatrix(:,2));

dphie= controlmatrix(:,3);
dphieb= controlmatrix(:,4);
dphimu= controlmatrix(:,5);
dphimub= controlmatrix(:,6);

[Cth,E] = meshgrid(cth,ene);


dPhie = reshape(dphie,length(cth),length(ene));
dPhieb = reshape(dphieb,length(cth),length(ene));
dPhimu = reshape(dphimu,length(cth),length(ene));
dPhimub = reshape(dphimub,length(cth),length(ene));


%Plot events
%%figure 
%%imagesc(n_ij)
figure('Renderer', 'painters', 'Position', [10 10 1000 800])

subplot(2,2,1)
pcolor(Cth,E, dPhie);
colorbar
set(gca, 'YScale', 'log');

subplot(2,2,2)
pcolor(Cth,E, dPhieb);
colorbar
set(gca, 'YScale', 'log');

subplot(2,2,3)
pcolor(Cth,E, dPhimu);
colorbar
set(gca, 'YScale', 'log');

subplot(2,2,4)
pcolor(Cth,E, dPhimub);
colorbar
set(gca, 'YScale', 'log');






