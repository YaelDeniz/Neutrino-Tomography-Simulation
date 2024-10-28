TrueSenvNumu = 'IntChi2/Int_cth_Chi2LLVP100myrs_cakenu1_100100100';
TrueSenvNue = 'IntChi2/Int_cth_Chi2LLVP100myrs_pancakenu0_5010050';

TrueChi2numu = readmatrix(TrueSenvNumu)
TrueChi2nue = readmatrix(TrueSenvNue)

TrueChi2(:,1) = TrueChi2numu(:,1)
TrueChi2(:,2) = TrueChi2numu(:,2) 
%TrueChi2(:,2) = TrueChi2numu(:,2) + TrueChi2nue(:,2)



figure
plot(TrueChi2(:,1),TrueChi2numu(:,2))
title("True Sensitivity to LLVP material[\mu-signal]")
xlabel("{\Delta n_e}/{n_e}")
ylabel("\Delta \chi^2")

yline(9,'--',"3\sigma")
yline(4,'--',"2\sigma")
yline(1,'--',"1\sigma")
%%saveas(gcf,'True/chi2LLVPmu.png')



figure
plot(TrueChi2(:,1),TrueChi2nue(:,2))
title("True Sensitivity to LLVP material[e-signal]")
xlabel("{\Delta n_e}/{n_e}")
ylabel("\Delta \chi^2")

yline(9,'--',"3\sigma")
yline(4,'--',"2\sigma")
yline(1,'--',"1\sigma")
%%saveas(gcf,'True/chi2LLVPe.png')

%% Full Sensitivity 
neq = linspace(-4,4,100);
chi2q = interp1(TrueChi2(:,1),TrueChi2(:,2),neq,'cubic')

figure
%%plot(TrueChi2(:,1),TrueChi2(:,2))
plot(neq,chi2q,'LineWidth',2.5,'Color','#fde725')

hold on
plot(neq,2*chi2q,'LineWidth',2.5,'Color','#21918c');
plot(neq,4*chi2q,'LineWidth',2.5,'Color','#440154');
hold off
axis([-3 3 0 10]);

title("True Sensitivity to LLVP material")
xlabel("{\Delta n_e}/{n_e}")
ylabel("\Delta \chi^2")
legend({"100 Mton-year","200 Mton-year","400 Mton-year"},'location','northwest')


yline(9,'--',"3\sigma(99%) CL",'HandleVisibility','off')
yline(4,'--',"2\sigma(95%) CL",'HandleVisibility','off')
yline(1,'--',"1\sigma(68%) CL",'HandleVisibility','off')