%% INTERPOLATION
clear; clc; 
load("flux.mat")
logEi = log(flux(:,1));

logMui = log(flux(:,2)*10^-2);

logei = log(flux(:,4)*10^-2);

E = linspace(1,100,1000);

LogE = log(E); %Logarithm of points

LogMu = interp1(logEi,logMui,LogE);
Loge = interp1(logEi,logei,LogE);

Mu = exp(LogMu);
e = exp(Loge);

Y1 = Mu.*(E.^3);
Y2 = e.*(E.^3);

plot(LogE,Y1)
hold on
plot(LogE,Y2)

%semilogx(E,Y);

