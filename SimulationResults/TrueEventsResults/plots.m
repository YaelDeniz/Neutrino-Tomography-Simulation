lamda = pi*1/(1.27*0.002)

dm2 = 2.5*10^-3;

theta = 8.588*pi/180;

G_F = 1.166378*10^-5 %%GeV^-2

% ne= Z/A*rho*N_A;

ne = 0;

Ve = sqrt(2)*G_F*ne;

eta = sqrt( ( sin(2*theta) )^2 + (cos(2*theta)-2*E*Ve/dm2)^2   )

m_matter = eta*dm2

