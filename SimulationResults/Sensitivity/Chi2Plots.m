TrueSenvNumu = 'True/TrueChi25010050.txt';
TrueSenvNue = 'True/TrueChi2_0_pancake5010050.txt';

TrueChi2numu = readmatrix(TrueSenvNumu)
TrueChi2nue = readmatrix(TrueSenvNue)

TrueChi2(:,1) = TrueChi2numu(:,1)
TrueChi2(:,2) = TrueChi2numu(:,2) + TrueChi2nue(:,2)

tiledlayout(2,2)
nexttile
plot(TrueChi2(:,1),TrueChi2numu(:,2))

title("True Sensitivity to LLVP material[\mu-signal]")
xlabel("{\Delta n_e}/{n_e}")
ylabel("\Delta \chi^2")

yline(9,'--',"3\sigma")
yline(4,'--',"2\sigma")
yline(1,'--',"1\sigma")



nexttile
plot(TrueChi2(:,1),TrueChi2nue(:,2))
title("True Sensitivity to LLVP material[e-signal]")
xlabel("{\Delta n_e}/{n_e}")
ylabel("\Delta \chi^2")

yline(9,'--',"3\sigma")
yline(4,'--',"2\sigma")
yline(1,'--',"1\sigma")


nexttile
plot(TrueChi2(:,1),TrueChi2(:,2))

title("True Sensitivity to LLVP material")
xlabel("{\Delta n_e}/{n_e}")
ylabel("\Delta \chi^2")
yline(9,'--',"3\sigma")
yline(4,'--',"2\sigma")
yline(1,'--',"1\sigma")

axis([-3 3 0 20])
hold on

plot(TrueChi2(:,1),2*TrueChi2(:,2))
plot(TrueChi2(:,1),4*TrueChi2(:,2))


hold off

