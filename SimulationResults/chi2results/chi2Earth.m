% 
% clear all;    clc;
% 
% prem = 'prem_44layers.txt'
% 
% PREMDATA = readmatrix(prem);
% 
% depth = PREMDATA(1:44,1)
% 
% file = 'chi2true/Sensitivity3Dmodel.csv';
% %file2 = 'chi2true/chi2_prem_44layers4040.csv';
% 
% CHI2DATA = readmatrix(file)
% %CHI2DATA2 = readmatrix(file2)
% 
% index = unique(CHI2DATA(:,1));
% 
% chi2_3 = CHI2DATA(1:44,2);
% %chi2_32 = CHI2DATA2(1:44,2);
% chi2_5 = CHI2DATA(45:88,2);
% chi2_10 = CHI2DATA(89:132,2);
% 
% semilogx(chi2_3,depth);
% hold on
% semilogx(chi2_5,depth,'-');
% semilogx(chi2_10,depth,'-');
% 
% 
% xlabel("\Delta \chi^2")
% 
% ylabel("Layer radius [km]")
% 
% xline(1,'-',"\Delta \chi^2 = 1")
% xline(4,'-',"\Delta \chi^2 = 4")
% xline(9,'-',"\Delta \chi^2 = 9")
% 
% yline(3480,'--','CMB')
% yline(3480+1000,'--','LLVP Limit')
% yline(6371-660,'--','660 Layer')
% yline(6371-410,'--','410 Layer')
% yline(6371-35,'--','Moho')
% 
% hold off
% 
% title("Sensitivity to the Density Contrats (True Events)")
% legend({'3%','5%','10%'},'Location','SouthEast')


%plot(chi2,depth);
%% 10 Pct


clear all;    clc;

%%prem = 'prem_44layers.txt'

%%PREMDATA = readmatrix(prem);

%%depth = PREMDATA(1:44,1)

chifile = 'chi2true/chi2test.csv';

chi2data = readmatrix(chifile)

index= chi2data(:,1)

chi2 = chi2data(:,2)


semilogx(chi2,index);
hold on

xlabel("\Delta \chi^2")

ylabel("Layer radius [km]")

xline(1,'-',"\Delta \chi^2 = 1")
xline(4,'-',"\Delta \chi^2 = 4")
xline(9,'-',"\Delta \chi^2 = 9")

% yline(3480,'--','CMB')
% yline(3480+1000,'--','LLVP Limit')
% yline(6371-660,'--','660 Layer')
% yline(6371-410,'--','410 Layer')
% yline(6371-35,'--','Moho')

hold off

title("Sensitivity to the Density Contrats (True Events)")




