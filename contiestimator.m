clear all;
clc;
%close all;
tic

fre_set = 100;
Hook    = [linspace(0.15,0.15,fre_set) linspace(0,0,fre_set) linspace(-0.15,-0.15,fre_set)];
maxtime = 30/numel(Hook); % periode

%% Parameter Ensemble Kalman Filter
P_Noise = 0.000125;
M_Noise = 0.005 ;

num_members     = 100;
num_iterations  = numel(Hook);

x_tr    = zeros(1,10); %initial value of state, 10 is the number of states in the ODE 
x_ini   = zeros(num_members,10); %ensemble of initial estimate of the state
w       = P_Noise*ones(1,10);  %process noise, zero mean, squared  standard deviation variance
z       = M_Noise*ones(1,1);  %measurement noise, zero mean, squared  standard deviation variance
uncer   = 10/100;
LER     = uncer;

%% DRAG FORCES
Loaded_Data = readtable('Drag 1KM 4 Segments.xlsx'); % in table format
drag        = Loaded_Data{:,5};                      % drag forces (lbf)
drag        = drag(~isnan(drag));

% Panjang in meter
Total = 1000;
ns    = 4;
LengthS = Total/ns;
LengthSNom = Total/ns - (LER*LengthS);
mod_el = 206842718795.3;          % Elastic modulus, N/m^2

% String Table, inches
ODm = 5 * 0.0254; %in meter now
ODj = 6.5 * 0.0254;
ODc = 6.5 * 0.0254;
ODb = 6.5 * 0.0254;
IDm = 4.276 * 0.0254;
IDj = 3.75 * 0.0254;
IDc = 2.813 * 0.0254;
IDb = 3.89 * 0.0254;
holeD = (8+(1/2)) * 0.0254;

% in SI unit
pipe_den = 7800.001722079;        %Pipe density, in kg/m^3
rho  = 1840;                      %Mud density, in kg/m^3

alp = ODm/holeD;
Kc = (alp^2 - sqrt((alp^4 +alp)/(1+alp)))/(1-alp^2);
n = 0.55; 
k = 0.13;
h1 = (holeD-ODm)/2; %The width of the annular space
h2 = (holeD-ODj)/2;
hC = (holeD-ODc)/2;
hB = (holeD-ODb)/2;
aa = 0.138; %Rheology parameter for annulus
bb = 0.312;
Au1 = 2*pi*0.95*LengthS*(ODm/2 + holeD/2);
Au2 = 2*pi*0.05*LengthS*(ODj/2 + holeD/2);
AuC = 2*pi*0.911*LengthS*(ODc/2 + holeD/2);
AuB = 2*pi*0.088*LengthS*(ODb/2 + holeD/2);

% Drillstring Parameters
Stiff_bha = mod_el*((ODb^2 - IDb^2)*pi/4)/(0.088*LengthS);
Stiff_coll= mod_el*((ODc^2 - IDc^2)*pi/4)/(0.911*LengthS);
Stiff_main= mod_el*((ODm^2 - IDm^2)*pi/4)/(0.95*LengthS);
Stiff_sec = mod_el*((ODj^2 - IDj^2)*pi/4)/(0.05*LengthS);

stiffness1 = Stiff_main*Stiff_sec/(Stiff_main+Stiff_sec);     %stiffness of drill pipe
stiffness2 = stiffness1;
stiffness3 = stiffness1;
stiffness4 = Stiff_bha*Stiff_coll/(Stiff_bha+Stiff_coll);
stiffness = [stiffness1 stiffness2 stiffness3 stiffness4];

mass_top  = 20000;
mass_bha  = ((ODb^2 - IDb^2)*pi/4)*(0.088*LengthS)*pipe_den;
mass_coll = ((ODc^2 - IDc^2)*pi/4)*(0.911*LengthS)*pipe_den;
mass_main = ((ODm^2 - IDm^2)*pi/4)*(0.95*LengthS)*pipe_den;
mass_sec  = ((ODj^2 - IDj^2)*pi/4)*(0.05*LengthS)*pipe_den;

mass1 = mass_main+mass_sec; %mass of drill pipe
mass2 = mass1;
mass3 = mass1;
mass4 = mass_bha+mass_coll;
mass  = [mass_top mass1 mass2 mass3 mass4];


r1 = stiffness1/mass_top;
r12= stiffness1/mass1;
r2 = stiffness2/mass1;
r23= stiffness2/mass2;
r3 = stiffness3/mass2;
r34= stiffness3/mass3;
r4 = stiffness4/mass3;
r45= stiffness4/mass4;
r = [r1 r12 r2 r23 r3 r34 r4 r45];

%% Drillstring Parameters Nominal
Stiff_bhaNom = mod_el*((ODb^2 - IDb^2)*pi/4)/(0.088*LengthSNom);
Stiff_collNom= mod_el*((ODc^2 - IDc^2)*pi/4)/(0.911*LengthSNom);
Stiff_mainNom= mod_el*((ODm^2 - IDm^2)*pi/4)/(0.95*LengthSNom);
Stiff_secNom = mod_el*((ODj^2 - IDj^2)*pi/4)/(0.05*LengthSNom);

stiffness1Nom = Stiff_mainNom*Stiff_secNom/(Stiff_mainNom+Stiff_secNom);     %stiffness of drill pipe
stiffness2Nom = stiffness1Nom;
stiffness3Nom = stiffness1Nom;
stiffness4Nom = Stiff_bhaNom*Stiff_collNom/(Stiff_bhaNom+Stiff_collNom);
stiffnessNom = [stiffness1Nom stiffness2Nom stiffness3Nom stiffness4Nom];

mass_bhaNom  = ((ODb^2 - IDb^2)*pi/4)*(0.088*LengthSNom)*pipe_den;
mass_collNom = ((ODc^2 - IDc^2)*pi/4)*(0.911*LengthSNom)*pipe_den;
mass_mainNom = ((ODm^2 - IDm^2)*pi/4)*(0.95*LengthSNom)*pipe_den;
mass_secNom  = ((ODj^2 - IDj^2)*pi/4)*(0.05*LengthSNom)*pipe_den;

mass1Nom = mass_mainNom+mass_secNom; %mass of drill pipe
mass2Nom = mass1Nom;
mass3Nom = mass1Nom;
mass4Nom = mass_bhaNom+mass_collNom;
massNom  = [mass_top mass1Nom mass2Nom mass3Nom mass4Nom];

r1nom = stiffness1Nom/mass_top;
r12nom= stiffness1Nom/mass1Nom;
r2nom = stiffness2Nom/mass1Nom;
r23nom= stiffness2Nom/mass2Nom;
r3nom = stiffness3Nom/mass2Nom;
r34nom= stiffness3Nom/mass3Nom;
r4nom = stiffness4Nom/mass3Nom;
r45nom= stiffness4Nom/mass4Nom;
rNom = [r1nom r12nom r2nom r23nom r3nom r34nom r4nom r45nom];

% EnKF Routine
[xtrue,xestm]=continuous_EnKF(mass,drag,ODm,holeD,Kc,rho,k,h1,n,ODj,h2,...
    ODc,hC,ODb,hB,aa,bb,Au1,Au2,AuC,AuB,Hook,num_members,x_tr,rNom,r,w,z,num_iterations,maxtime);

error = abs(xtrue)-abs(xestm);
% 
figure()
plot(xtrue(:,10),'LineWidth',1);
hold on;
plot(xestm(:,10),'LineWidth',1);
set(gca,'FontSize',24)
xlabel('Time (sec)','FontSize', 24)
ylabel('BHA Velocity (m/s)','FontSize', 24)
legend('Real Velocity','Estimated Velocity')
text(5,-3,['Number of Ensemble = ',num2str(num_members)],'fontsize',18);
text(5,-2.5,['Measurement STD = ',num2str(M_Noise)],'fontsize',18);
text(5,-2,['Process STD = ',num2str(P_Noise)],'fontsize',18);
text(5,-1.5,['Frequency = ',num2str(fre_set)],'fontsize',18);
hold off
toc

% save('data kecepatan','xtrue','xestm');