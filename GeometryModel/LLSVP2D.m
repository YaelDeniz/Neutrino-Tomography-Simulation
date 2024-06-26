
%% Geometry of 1D Earth model
 
R=6378; % km

% Create a vectortheta.

omega=linspace(0,2*pi,200); 

rc = 6378; %Layer 1: Crust
xc=(rc/R)*cos(omega)-(rc/R); % Generate x-coordinates.
yc=(rc/R)*sin(omega); % Generate y-coordinate.
 
rup = 6343; %Layer 1: Crust
xup=(rup/R)*cos(omega)-(rc/R); % Generate x-coordinates.
yup=(rup/R)*sin(omega); % Generate y-coordinate.

rlm = 5671; %Layer 2: upper mantle
xlm=(rlm/R)*cos(omega)-(rc/R); % Generate x-coordinates.
ylm=(rlm/R)*sin(omega); % Generate y-coordinate.

roc = 3486; %Layer 1: Crust
xoc=(roc/R)*cos(omega)-(rc/R); % Generate x-coordinates.
yoc=(roc/R)*sin(omega); % Generate y-coordinate.

ric = 1216; %Layer 1: Crust
xic=(ric/R)*cos(omega) - (rc/R); % Generate x-coordinates.
yic=(ric/R)*sin(omega) ; % Generate y-coordinate.

%% LLSVP vectorized
hllsvp = 700;
omega_llsvp = 50*(pi/180);
rllsvp = roc + hllsvp ;

y1=cos(omega_llsvp); % Generate x-coordinates.
x1=sin(omega_llsvp) - (rc/R); % Generate y-coordinate.

y2=cos(omega_llsvp); % Generate x-coordinates.
x2=-1.0*sin(omega_llsvp) - (rc/R); % Generate y-coordinate.


%% LLSVP arc

R=6378; % km


% Create a vectortheta.
omega_min = -50*(pi/180);

omega_max = 50*(pi/180);

omega_region=linspace(omega_min,omega_max,200);

rllsvp = roc + hllsvp ;

xllsvp=(rllsvp/R)*sin(omega_region) - (rc/R); % Generate x-coordinates.

yllsvp=(rllsvp/R)*cos(omega_region) ; % Generate y-coordinate.

%% Neutrino trayectories
 
%%zenith_min = pi-asin((roc-1800)/R); %Skimming the CMB
zenith_max = acos(-0.7095);

zenith_min = acos(-0.837);

xz1=(2)*cos(zenith_min); % Generate x-coordinates.
yz1=(2)*sin(zenith_min); % Generate y-coordinate.

xz2=(2)*cos(zenith_max); % Generate x-coordinates.
yz2=(2)*sin(zenith_max); % Generate y-coordinate.

%% Intersection points

zenith = zenith_min

% yi = -1/( tan(omega_max)-tan(zenith) );
% xi = -1*( tan(zenith) )/( tan(omega_max)-tan(zenith) );
% 
% yf = 1/( tan(omega_max) + tan(zenith) );
% xf = tan(zenith)/( tan(omega_max) + tan(zenith) );

%Last
xx2 = cot(omega_max)/(tan(zenith)-cot(omega_max));
yy2 = tan(zenith)*xx2;

%Fisrt
xx1 = -1*cot(omega_max)/(tan(zenith)+cot(omega_max));
yy1 = tan(zenith)*xx1;

%Crust
xxo = -2/( 1+tan(zenith)*tan(zenith) );
yyo = tan(zenith)*xxo;

L1 = sqrt ( (xx1-xxo)^2 + (yy1-yyo)^2 )*6386
L2 = sqrt ( (xx2-xxo)^2 + (yy2-yyo)^2 )*6386

Rtest = sqrt( (xx1 + 1 )^2 + (yy1 - 0)^2  )*6386

%% plot the circle.
plot(xc,yc,'k'); 

hold on
plot(xup,yup,'k'); 
plot(xlm,ylm,'k--'); 
plot(xoc,yoc,'k-.'); 
plot(xic,yic,'k.'); 

plot([-rc/R x1], [0 y1], 'r'); 
plot([-rc/R x2], [0 y2], 'r'); 
plot(xllsvp,yllsvp,'r')

plot([0 xz1], [0 yz1], 'b--'); 
plot([0 xz2], [0 yz2], 'b');   

%plot([xi xf], [yi yf], 'k+')
plot([xxo xx1 xx2],[yyo yy1 yy2], 'ko')
%%plot( xx1 , yy1 , 'kx')

hold off
 
 % Set equal scale on axes.
axis('equal');

title('Earth model')
