Senv2Layer = 'Sensitivity2Layers/EarthSenvMeff.txt';
PREM = readmatrix('prem_44layers.txt');
depth = PREM(:,1)

chi2data = readmatrix(Senv2Layer)

chi2_e = chi2data(:,1)
chi2_mu = chi2data(:,2)
chi2_tot= chi2data(:,3)

layer = linspace(1,44,44);

figure
semilogx(chi2_e,depth)
xlabel("\Delta \chi^2")
ylabel("Depth[km]")
title("Sensitivity of \nu_e-events to 10%")
hold on
yline(630.9,'--r',"Innermost IC")
yline(1121.5,"--","IC")
yline(3480.0,"--","CMB")
yline(5711,"--","660/LM")
yline(5961,"--","410")
yline(6346.6,"--","UM")
hold off
saveas(gcf,'Sensitivity2Layers/chi2elayers10.png')



figure
semilogx(chi2_mu,depth)
xlabel("\Delta \chi^2")
ylabel("Depth[km]")
title("Sensitivity of \nu_\mu-events to 10%")
hold on
yline(630.9,'--r',"Innermost IC")
yline(1121.5,"--","IC")
yline(3480.0,"--","CMB")
yline(5711,"--","660/LM")
yline(5961,"--","410")
yline(6346.6,"--","UM")
hold off
saveas(gcf,'Sensitivity2Layers/chi2mulayers10.png')

figure
semilogx(chi2_tot,depth)
xlabel("\Delta \chi^2")
ylabel("Depth[km]")
title("Sensitivity of joint signal (\nu_\mu + \nu_e)  to 10%")
hold on
yline(630.9,'--r',"Innermost IC")
yline(1121.5,"--","IC")
yline(3480.0,"--","CMB")
yline(5711,"--","660/LM")
yline(5961,"--","410")
yline(6346.6,"--","UM")
hold off
saveas(gcf,'Sensitivity2Layers/chi2totlayers10.png')

%% XSEc
ene = logspace(0,2,20)
nuxsec = 0.75E-38.*ene
nubxsec = 0.35E-38.*ene



