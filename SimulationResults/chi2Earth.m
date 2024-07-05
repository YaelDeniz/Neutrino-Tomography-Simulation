
%% MyData
clear all;    clc;

prem = '/home/dehy0499/OscProb/PremTables/prem_44layers.txt'

PREMDATA = readmatrix(prem);

depth = PREMDATA(:,1)

file = 'chi2results/chi2true/chi2Earth5pct100100.csv';

CHI2DATA = readmatrix(file)

index = CHI2DATA(:,1);

chi2 = CHI2DATA(:,2);


semilogx(chi2,depth);


