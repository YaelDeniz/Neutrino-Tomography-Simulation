clear all
%% Trapezoidal Rule

%% Simpsons 1/3

xmax=1
xmin=0

ymax=1
ymin=0

h = (xmax-xmin)/2;
k = (ymax-ymin)/2;

x = xmin:h:xmax
y = ymin:k:ymax

% 
 i = F(x(1),y(1)) + F(x(1),y(3)) + F(x(3),y(1)) + F(x(3),y(3))

 ii = F(x(1),y(2)) + F(x(2),y(1)) + F(x(2),y(3)) + F(x(3),y(2)) 

 iii = F(x(2),y(2))
format long
 I = (h*k/9)*( i+ 4*ii + 16*iii)

%% Erf Energy
 Eo = linspace(4,6,10)
 a = 1;

 Et = linspace(4,6,200);
 Emin=Et(1)
 Emax=Et(2)
 e = linspace (Emin,Emax,10);

 syms w(E)
 w(E) =(1/2)*(erf( (Eo(2)-E)/(sqrt(2)*a*E) ) - erf( (Eo(1)-E)/(sqrt(2)*a*E) ) )

 fplot(w(E),[3.5 5]);

 ylim([0.01 0.025])

 %% Erf Angle
 

 
 To = linspace(10,30,6)
 a = 3.65;
 
 T = linspace(10,30,200);
 Et = linspace(4,6,200);
 Emin=Et(1)
 Emax=Et(2)

 Tmin=T(1)
 Tmax=T(2)

 e = linspace (Emin,Emax,10);
 t = linspace (Tmin,Tmax,10);

 [E,TH] = meshgrid(e,t);

 %W =erf( TH./(1./sqrt(E)) )

 W =(1/2)*(erf( (To(2)-TH)./(sqrt(2)*a./sqrt(E)) ) - erf( (To(1)-TH)./(sqrt(2)*a./sqrt(E)) ) )
 
 surf(E,TH,W)




 

function f = F(x,y)
    %f =  sin(x*y)/(1+x*y);
    f =  1/(1 - (x*y)^2);
end


