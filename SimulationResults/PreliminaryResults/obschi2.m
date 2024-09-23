ObSenvNumu = 'ObsChi2/ObsChi2LLVP_pancakenu1_2010020.txt';
ObsSenvNue = 'ObsChi2/ObsChi2LLVP_pancakenu0_2010020.txt';

ObsChi2numu = readmatrix(ObSenvNumu)
ObsChi2nue = readmatrix(ObsSenvNue)

ObsChi2(:,1) = ObsChi2numu(:,1)
ObsChi2(:,2) = ObsChi2numu(:,2) + ObsChi2nue(:,2)

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



figure
plot(ObsChi2(:,1),ObsChi2(:,2))

title("Obs Sensitivity to LLVP material")
xlabel("{\Delta n_e}/{n_e}")
ylabel("\Delta \chi^2")
yline(9,'--',"3\sigma")
yline(4,'--',"2\sigma")
yline(1,'--',"1\sigma")

axis([-3 3 0 10])
hold on

plot(ObsChi2(:,1),2*ObsChi2(:,2))
plot(ObsChi2(:,1),4*ObsChi2(:,2))


hold off