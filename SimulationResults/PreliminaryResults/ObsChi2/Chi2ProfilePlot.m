

File = ['ObsChi2_Simulation_pancake_100Mton_2-10GeV_128-171Zen_-55-55Az/' ...
    'Obs_cth_Chi2LLVP_densellvp_100myrs_pancakenu1_Obs40Zen100Az40Enu_Int100Zen100Az100Enu.csv'];


SensitivityNumu = csvread(File)

ne = SensitivityNumu(:,2)
chi2numu = SensitivityNumu(:,3) 

figure
plot(ne,chi2numu,'LineWidth',2.5,'Color','#fde725')

hold on

plot(ne,3*chi2numu,'LineWidth',2.5,'Color','#21918c');
plot(ne,4*chi2numu,'LineWidth',2.5,'Color','#440154');
%plot(ne,100*chi2q,'LineWidth',2.5,'Color','red');

hold off

title("Observed Sensitivity to reject a dense LLVP material[2%;\mu-signal]")
xlabel("{\Delta n_e}/{n_e}")
ylabel("\Delta \chi^2")
legend({"100 Mton-year","300 Mton-year","400 Mton-year"},'location','northwest')


yline(9,'--',"3\sigma(99%) CL",'HandleVisibility','off')
yline(4,'--',"2\sigma(95%) CL",'HandleVisibility','off')
yline(1,'--',"1\sigma(68%) CL",'HandleVisibility','off')



