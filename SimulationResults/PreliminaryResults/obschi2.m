%ObSenvNumu = 'ObsChi2/ObsChi2LLVP_pancakenu1_2010020.txt';
%ObsSenvNue = 'ObsChi2/ObsChi2LLVP_pancakenu0_2010020.txt';

ObSenvNumu = 'ObsChi2/ObsChi2_cth_LLVP_cakenu1_2010020_100100100.txt';
ObsSenvNue = 'ObsChi2/ObsChi2_cth_LLVP_pancakenu0_2010020_6010060.txt';

ObsChi2numu = readmatrix(ObSenvNumu)
ObsChi2nue = readmatrix(ObsSenvNue)

ObsChi2(:,1) = ObsChi2numu(:,1)
ObsChi2(:,2) = ObsChi2numu(:,2) 

figure
plot(ObsChi2(:,1),ObsChi2numu(:,2))
title("Obs Sensitivity to LLVP material[\mu-signal]")
xlabel("{\Delta n_e}/{n_e}")
ylabel("\Delta \chi^2")

yline(9,'--',"3\sigma")
yline(4,'--',"2\sigma")
yline(1,'--',"1\sigma")



figure
plot(ObsChi2(:,1),ObsChi2nue(:,2))
title("Obs Sensitivity to LLVP material[e-signal]")
xlabel("{\Delta n_e}/{n_e}")
ylabel("\Delta \chi^2")

yline(9,'--',"3\sigma")
yline(4,'--',"2\sigma")
yline(1,'--',"1\sigma")


%% Full Sensitivity
neq = linspace(-4,4,100);
chi2q = interp1(ObsChi2(:,1),ObsChi2(:,2),neq,'cubic')

figure
plot(neq,chi2q,'LineWidth',2.5,'Color','#fde725')

hold on
plot(neq,2*chi2q,'LineWidth',2.5,'Color','#21918c');
plot(neq,4*chi2q,'LineWidth',2.5,'Color','#440154');
plot(neq,100*chi2q,'LineWidth',2.5,'Color','red');
hold off
axis([-3 3 0 10]);

title("Observed Sensitivity to LLVP material")
xlabel("{\Delta n_e}/{n_e}")
ylabel("\Delta \chi^2")
legend({"100 Mton-year","200 Mton-year","400 Mton-year"," 10 Gigaton-year"},'location','northwest')


yline(9,'--',"3\sigma(99%) CL",'HandleVisibility','off')
yline(4,'--',"2\sigma(95%) CL",'HandleVisibility','off')
yline(1,'--',"1\sigma(68%) CL",'HandleVisibility','off')

