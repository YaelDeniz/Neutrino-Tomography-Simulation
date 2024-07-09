
%% MyData
clear all;    clc;

prem = '/home/dehy0499/OscProb/PremTables/prem_44layers.txt'

PREMDATA = readmatrix(prem);

depth = PREMDATA(:,1)

file = 'chi2true/chi2Earth.csv';

CHI2DATA = readmatrix(file)

index = CHI2DATA(:,1);

chi2 = CHI2DATA(:,2);

semilogx(chi2,depth);

xlabel("\Delta \chi^2")
ylabel("Layer radius [km]")

yline(3480,'--','CMB')
yline(5700,'--','Mantle')

title("Sensitivity to the 10% Density Contrats")


%plot(chi2,depth);


