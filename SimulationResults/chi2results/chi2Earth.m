
%% MyData
clear all;    clc;

prem = '/home/dehy0499/OscProb/PremTables/prem_44layers.txt'

PREMDATA = readmatrix(prem);

depth = PREMDATA(:,1)

file = 'chi2true/chi2Earth15_4040.csv';

CHI2DATA = readmatrix(file)

index = unique(CHI2DATA(:,1));

chi2_3 = CHI2DATA(1:44,2);
chi2_5 = CHI2DATA(45:88,2);
chi2_10 = CHI2DATA(89:132,2);

semilogx(chi2_3,depth);
hold on
semilogx(chi2_5,depth);
semilogx(chi2_10,depth);


xlabel("\Delta \chi^2")
ylabel("Layer radius [km]")

yline(3480,'--','CMB')
yline(5700,'--','Mantle')
hold off
title("Sensitivity to the 10% Density Contrats")


%plot(chi2,depth);


