

eta_min=135*pi/180.0;

eta_max=151*pi/180.0;

etaOmin = (180-45)*pi/180.0;

etaOmax = (180-28)*pi/180.0;

eta_o = eta_max;


ene = 1;

syms f(eta,E); 

%%f(eta,E) = erf( (eta-eta_min)/(0.25/sqrt(E)) ) - erf( (eta-eta_max)/(0.25/sqrt(E)) );

f(eta,E) = (1 /( sqrt(2*pi)*(0.25/sqrt(E))   ) )*exp(-0.5*( (eta- eta_o)/( 0.25/sqrt(E) ) )^2 );



close all; 
hold on
fplot( f(eta, 1 ), [0 2*pi],'-r');
fplot( f(eta, 5 ), [0 2*pi],'-b');
fplot( f(eta, 10 ), [0 2*pi],'-g');

xline( eta_o ,'--k' );

xline( 0 ,'-.k' );
xline( pi ,'-.k' );
 
%  xline( eta_o - 4*0.25/sqrt(1),'*r' );
%  xline( eta_o - 4*0.25/sqrt(5),'*b');
%  xline( eta_o - 4*0.25/sqrt(10),'*g');
% 
%   xline( eta_o + 4*0.25/sqrt(1),'*r' );
%  xline( eta_o + 4*0.25/sqrt(5),'*b');
%  xline( eta_o + 4*0.25/sqrt(10),'*g');


% 
 xline(etaOmin - 2*0.25/sqrt(1),'--r' );
 xline(etaOmin - 2*0.25/sqrt(5),'--b');
 xline(etaOmin - 2*0.25/sqrt(10),'--g');

  xline(etaOmax + 2*0.25/sqrt(1),'--r' );
 xline(etaOmax + 2*0.25/sqrt(5),'--b');
 xline(etaOmax + 2*0.25/sqrt(10),'--g');






hold off
%%

sigma = @(E) 10 + 4*(0.05*E + 0.1*sqrt(E)) - E   ;

x0 = 100;


x = fzero(sigma,x0)

%%

fplot(@(E) E, [0.1 20])
hold on
fplot(@(E) 1.0 - 4.0*(0.05*E +0.1*sqrt(E)) , [0.1 20],'--r')
fplot(@(E) 10.0 + 4.0*(0.05*E +0.1*sqrt(E)) , [0.1 20],'--b')
hold off