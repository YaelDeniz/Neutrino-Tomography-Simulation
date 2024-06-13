%% INTERPOLATION
clear; clc; close all;
load("flux2.mat")

flux = flux2
%% Logarithmic transformstion

npts = 2000; 

%logEi = log10(flux(:,1));
logEi = log10(logspace(-1,4,101));
%%Ei = logspace(0,4,101);

%%logEi=log10(Ei);


logMui = log10(flux(:,2)*10^-4);

logMuBari = log10(flux(:,3)*10^-4);

logei = log10(flux(:,4)*10^-4);

logeBari = log10(flux(:,5)*10^-4);

%E = linspace(1,100,npts);

%LogE = log10(E); %Logarithm of points
LogE = linspace(0,2,100);
%%Inerpoaltion

LogMu = interp1(logEi,logMui,LogE,"cubic");
LogMuBar = interp1(logEi,logMuBari,LogE,"cubic");
Loge = interp1(logEi,logei,LogE,"cubic");
LogeBar = interp1(logEi,logeBari,LogE,"cubic");

Mu = 10.^(LogMu);
MuBar = 10.^(LogMuBar);
e = 10.^(Loge);
eBar = 10.^(LogeBar);


numu= Mu.*((10.^LogE).^3);
numubar = MuBar.*((10.^LogE).^3);
nue = e.*((10.^LogE).^3);
nuebar = eBar.*((10.^LogE).^3);

plot(LogE,numu,"blue")

hold on

plot(LogE,numubar,"-.blue")
plot(LogE,nue,"red")
plot(LogE,nuebar,"-.red")

hold off
grid on

%semilogx(E,Y);

