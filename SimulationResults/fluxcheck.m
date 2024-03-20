%FLUX 1
FluxOk = DOlivoFlux
FluxOk(:,2) = 10^4*DOlivoFlux(:,2)./DOlivoFlux(:,1);
%honda 


E = FluxOk(:,1);
N = FluxOk (:,2);
EE = logspace(0,4, 1000);
NN =  ppval(E,N,EE);
semilogx(E,N,'+',EE,NN)

% hold on
% E2 = Table(:,1);
% N2 = Table(:,2);
% EE2 = logspace(0,4, 1000);
% NN2 =  spline(E2,N2,EE2);
% semilogx(E2,N2,'o',EE2,NN2)

legend('DOlivo','Interpol1','Honda','Interpol2')
