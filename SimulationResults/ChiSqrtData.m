% True data
chisqrt_true = 'chi2results/AsimovFixed_1_200_200_1.000000_10.000000.csv';
ChiSqrt_true = readmatrix(chisqrt_true);

n_det = 1;


rho_t = ChiSqrt_true(:,2);
chi_t = n_det*ChiSqrt_true(:,3);

figure
plot(rho_t,chi_t,'-','LineWidth',2.5)
set(gca,'FontSize',12, 'FontName', 'Courier')
hold on

plot(rho_t,2*chi_t,'-','LineWidth',2.5)

plot(rho_t,4*chi_t,'-','LineWidth',2.5)

%%Confidence levels

ones=yline(1,'--k','1\sigma= 68.3%','LineWidth',1.5)
ones.FontSize=20;

twos=yline(4,'--k','2\sigma=95.45%','LineWidth',1.5)
twos.FontSize=20;


hold off


xlabel('\fontsize{25} \Delta n_e/ n_e[%]');
ylabel('\fontsize{25} \Delta\chi^2');
title('Perfect resolution','FontSize',25)

legend({'100 Mton-years','200 Mton-years', '400 Mton-years'},'Location','northwest')

%% STATISTICAL ANALYSIS OF OBSERVED EVENTS

n_det = 1;

chisqrt_obs = 'chi2results/AsimovObsNextGen1_20_20_50_50_1.000000_10.000000.csv';
ChiSqrt_obs = readmatrix(chisqrt_obs);


rho_obs = ChiSqrt_obs(:,4);
chi_obs = n_det*ChiSqrt_obs(:,5);

figure
plot(rho_obs,chi_obs,'-','LineWidth',2.5)
set(gca,'FontSize',12, 'FontName', 'Courier')
hold on
% 
plot(rho_obs,2*chi_obs,'-','LineWidth',2.5)

plot(rho_obs,4*chi_obs,'-','LineWidth',2.5)



%%Confidence levels

ones=yline(1,'--k','1\sigma= 68.3%','LineWidth',1.5)
ones.FontSize=20;

twos=yline(4,'--w','2\sigma=95.45%','LineWidth',1.5)
twos.FontSize=20;



hold off



legend({'100 Mton-years','200 Mton-years', '400 Mton-years'},'Location','northwest')
xlabel('\fontsize{25} \Delta n_e/ n_e[%]');
ylabel('\fontsize{25} \Delta\chi^2');
title('Realistic resolution','FontSize',25)
ylim([0 4])






%% STATISTICAL ANALYSIS OF OBSERVED EVENTS:Det_res
% 
% n_det = 1;
% 
% chisqrt_obsp = 'chi2results/ObsDiff_Para_1_50_50_50_50_1.000000_10.000000.csv';
% ChiSqrt_obsp = readmatrix(chisqrt_obsp);
% 
% 
% n = 11
% 
% rho_obsp = ChiSqrt_obsp(1:n,2);
% chi_obsp = n_det*ChiSqrt_obsp(1:n,5);
% 
% figure
% plot(rho_obsp,chi_obsp,'r-','LineWidth',2.5)
% hold on
% 
% % plot(rho_obsp,10*2*chi_obsp,'r-o','LineWidth',2.5)
% % plot(rho_obsp,10*3*chi_obsp,'r-|','LineWidth',2.5)
% % plot(rho_obsp,10*4*chi_obsp,'r-*','LineWidth',2.5)
% 
% 
% % plot(rho2,chi2,'*')
% % plot(rho3,chi3)
% % plot(rho4,chi4)
% 
% yline(1,'--k','LineWidth',1.5)
% yline(4,'--k','LineWidth',1.5)
% yline(9,'--k','LineWidth',1.5)
% % yline(21.66,'--b')
% 
% hold off
% 
% 
% xlabel(' \fontsize{20} \alpha_E=\alpha_\eta');
% ylabel('\fontsize{30} \Delta\chi^2');
% title('\fontsize{20} \Delta\chi^2 profiles for observed events:Detector Resolution','Generic detector- 10MTon & NO 100%\Delta\rho/\rho=2')

