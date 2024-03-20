
%% Gaussian
s = 1;
t =-3*s:0.25:3*s;
[X,Y] = meshgrid(t);

F = (1/(2*pi*s))*exp( -(X.^2 + Y.^2)/(2*s)   );

surf(X,Y,273*F)

%% 2D Convolution

% I = [1 2 3 4; 5 6 7 8; 9 10 11 12 ; 13 14 15 16]
% 
% K = [1 0 1; 0 -1 1; -2 1 1]
% 
% O = zeros(2)

I = [ 2 1 1 0 5; 4 2 5 5 0; 3 5 7 8 0; 0 7 4 3 2; 5 2 1 1 0];

K = [1 2 0; 0 0 2; 0 2 2];

O = zeros(3);



for i = 1:3
    for j = 1:3

        
            
        % Write elemtns of matrix
            for k = 1:3
                for l = 1:3

                O(i,j) =  O(i,j) + I(i+k-1,j+l-1)*K(k,l);

                end
            end
          

    end
end

tiledlayout(1,3)
nexttile
imagesc(I)
axis square
caxis([0 10])
colormap(jet(10))
colorbar
title('Input')

nexttile
imagesc(K)
title('Kernel')
axis square
caxis([0 10])
colormap(jet(10))
colorbar

nexttile
imagesc(O)
axis square

colormap(jet(10))
colorbar
title('Output')

%% Effective mass

E = logspace(0,1,10)
Th = linspace(pi/2,pi , 10);
tiledlayout(2,3)
for i=1:5
    e = linspace(E(i),E(i+1),101);
    th = linspace(Th(i),Th(i+1),101);
    
    [X,Y] = meshgrid(e,th);

    Etm = median(e) %Median true energy
    thtm = median(th)
    se = 0.266*Etm/( Etm^(0.171) - 0.604 );
    sth = 3.65/(Etm^(1.05) + 5.00);
    
% %     rE = ( 1/(sqrt(2*pi)*se) )*exp( (-1/2)*( (e-Etm)/se ).^2 );
% %     rTH = ( 1/(sqrt(2*pi)*sth) )*exp( (-1/2)*( (th-thtm)/sth ).^2 );
    
    rE = ( 1/(sqrt(2*pi)*se) )*exp( (-1/2)*( (X-Etm)/se ).^2 );
    rTH = ( 1/(sqrt(2*pi)*sth) )*exp( (-1/2)*( (Y-thtm)/sth ).^2 );

    R = rE.*rTH;
    
    nexttile
    surf(X,Y,R)
    
    %imagesc(R)
    colorbar
%     nexttile
%     plot(e,rE);
%     hold on
%     xline(Etm,'--')
%     xline(mean(e),'-')
%     hold off

    
end
%% Resolution functions
clear;
close all;
e = linspace(1,40,100);
th = linspace(0.5,1,100);

sd = 3.65./(e.^(1.05) +5.00);

curve1_1 = 0.6 - sd/pi;
curve1_2 = 0.6 + sd/pi;

curve2_1 = 0.75 - sd/pi;
curve2_2 = 0.75 + sd/pi;

curve3_1 = 0.9 - sd/pi;
curve3_2 = 0.9 + sd/pi;

syms sde(e)

sde(e) = 0.266*e/(e^(0.171) - 0.604);

hold on
% Angular
plot(e,curve1_1,'-b');
yline (0.6,'--b');
plot(e,curve1_2,'-b');

plot(e,curve2_1,'-b');
yline (0.75,'--b');
plot(e,curve2_2,'-b');

plot(e,curve3_1,'-b');
yline (0.9,'--b');
plot(e,curve3_2,'-b');

%Energy
xline(3,'--r')
val = 3 - sde(3);
xline(10,'--r')
xline(30,'--r')

hold off


set(gca, 'XScale', 'log')
xlim([1 40])
ylim([0.5 1])
xlabel('E/GeV')
ylabel('\theta / \pi')





%plot(E,meff, "-")
%% Error functions
e = linspace(1,10,100);

ei=e(1);
eii=e(2);

em = (e(1)+e(2))/2;

syms w(E) sd(x)

sd(x) = 3.65/(x^(1.05) + 5.00); %Width of angular resolution function
%sd(x) = 0.01; %Width of angular resolution function

w(E) = (1/2)*( erf( (eii-E)/(sqrt(2)*sd(E))  ) - erf( (ei-E)/(sqrt(2)*sd(E)) )  );

%%Limits
Emin = 1;
Emax = 40;

fplot(sd(E), [Emin Emax])
set(gca, 'XScale', 'log')

