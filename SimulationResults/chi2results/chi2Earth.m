
%% MyData
clear all;    clc;

% prem = '/home/dehy0499/OscProb/PremTables/prem_44layers.txt'

%PREMDATA = readmatrix(prem);

% depth = PREMDATA(:,1)

file = 'chi2true/chi2Earth_yael.csv';

CHI2DATA = readmatrix(file)

index = CHI2DATA(:,1);

chi2 = CHI2DATA(:,2);

semilogx(chi2,index);

%plot(chi2,depth);


