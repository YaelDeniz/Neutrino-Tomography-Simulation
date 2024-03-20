

%% Geometry of 1D Earth model
 
R=6378; % km

% Create a vectortheta.

omega=linspace(0,2*pi,200); 

rc = 6378; %Layer 1: Crust
xc=(rc/R)*cos(omega); % Generate x-coordinates.
yc=(rc/R)*sin(omega) + (rc/R); % Generate y-coordinate.
 
rup = 6343; %Layer 1: Crust
xup=(rup/R)*cos(omega); % Generate x-coordinates.
yup=(rup/R)*sin(omega) + (rc/R); % Generate y-coordinate.

rlm = 5671; %Layer 2: upper mantle
xlm=(rlm/R)*cos(omega); % Generate x-coordinates.
ylm=(rlm/R)*sin(omega) + (rc/R); % Generate y-coordinate.

roc = 3486; %Layer 1: Crust
xoc=(roc/R)*cos(omega); % Generate x-coordinates.
yoc=(roc/R)*sin(omega) + (rc/R); % Generate y-coordinate.

ric = 1216; %Layer 1: Crust
xic=(ric/R)*cos(omega); % Generate x-coordinates.
yic=(ric/R)*sin(omega) + (rc/R); % Generate y-coordinate.

%% LLSVP vectorized
hllsvp = 700;
omega_llsvp = 50*(pi/180);
rllsvp = roc + hllsvp ;

x1=cos(omega_llsvp); % Generate x-coordinates.
y1=sin(omega_llsvp) + (rc/R); % Generate y-coordinate.

x2=cos(omega_llsvp); % Generate x-coordinates.
y2=-1.0*sin(omega_llsvp) + (rc/R); % Generate y-coordinate.


%% LLSVP arc

R=6378; % km


% Create a vectortheta.
omega_min = -50*(pi/180);

omega_max = 50*(pi/180);

omega_region=linspace(omega_min,omega_max,200);

rllsvp = roc + hllsvp ;

xllsvp=(rllsvp/R)*cos(omega_region); % Generate x-coordinates.

yllsvp=(rllsvp/R)*sin(omega_region) + (rc/R); % Generate y-coordinate.

%% Neutrino trayectories
 
zenith_min = pi/2 -asin((roc-1800)/R); %Skimming the CMB
zenith_max = pi/2 - asin((rllsvp+500)/R);

xz1=(2)*cos(zenith_min); % Generate x-coordinates.
yz1=(2)*sin(zenith_min); % Generate y-coordinate.

xz2=(2)*cos(zenith_max); % Generate x-coordinates.
yz2=(2)*sin(zenith_max); % Generate y-coordinate.

%% Intersection points
xi = -1/( tan(omega_max)-tan(zenith_min) );
yi = -1*( tan(zenith_min) )/( tan(omega_max)-tan(zenith_min) );

xf = 1/( tan(omega_max) + tan(zenith_min) );
yf = tan(zenith_min)/( tan(omega_max) + tan(zenith_min) );


%% plot the circle.
plot(xc,yc,'k'); 

hold on
plot(xup,yup,'k'); 
plot(xlm,ylm,'k--'); 
plot(xoc,yoc,'k-.'); 
plot(xic,yic,'k.'); 

plot([0 x1], [rc/R y1], 'r'); 
plot([0 x2], [rc/R y2], 'r'); 
plot(xllsvp,yllsvp,'r')

plot([0 xz1], [0 yz1], 'b'); 
plot([0 xz2], [0 yz2], 'b');   

plot([xi xf], [yi yf], 'g+')

hold off
 
 % Set equal scale on axes.
axis('equal');

title('Earth model')
