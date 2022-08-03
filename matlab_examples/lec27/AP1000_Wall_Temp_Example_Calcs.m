%% Onset of Nucleate Boiling Demo
clear
clc
close 'all'

%% Put EasyProp on Python path
% put location of EasyProp.py module on the python search path
if count(py.sys.path,' ') == 0  % <-- see if desired directory is on path
    insert(py.sys.path,int32(0),' '); %<-- if not; add it.
end

water = py.EasyProp.simpleFluid('Water','USCS');


%% AP1000 Given Data
Q_dot_total = 11.601e9*(1/3600); % BTU/s, total core heat output
P_nom = 2200; % psia, nominal system pressure
T_avg = 578.1; % F, average coolant temperature in the core (use for nominal coolant density) 
mu_nom = 0.23*(1/3600); % lbm/ft-s
sigma_nom = 0.009743; % lbf/ft, surface tension
rho_nom = 1./water.v_pT(P_nom,T_avg); % lbm/ft^3, nominal density
Pr_nom = water.Prandtl_pT(P_nom,T_avg); % dimensionless, nominal Prandtl number
k_nom = 0.32*(1/3600); % BTU/s-ft-R
T_sat = water.Tsat_p(P_nom); % F, saturation temp in core
h_sat = water.hL_p(P_nom); % BTU/lbm, enthalpy of saturated liquid at nominal pressure
h_fg = water.hV_p(P_nom) - water.hL_p(P_nom); % BTU/lbm
m_dot_total = 113.5e6*(1/3600); % lbm/s, flow core flow rate
N_rods = 157*17*17; 
D_co = 0.374/12; % ft, clad outside diameter
t_clad = 0.0225/12; % ft, clad thickness
D_gap = 6.5e-3/12; % ft, diametral gap
D_fuel = 0.3225/12; % ft, diameter of fuel (consistent with clad dimensions and diametral gap)
Pitch = 0.496/12; % ft, fuel rod pitch
H = 168/12; % ft, height of active fuel region
q_p_peak = (14.9/1.5)*3412*(1/3600); % kW/ft --> BTU/s-ft, axial peak linear heat rate
T_c = 535; % F, nominal core inlet temperature

%% calculated results
m_dot = m_dot_total/N_rods; % lbm/s, flow rate per fuel rod
A_flow = Pitch^2 - (pi/4)*D_fuel^2; % ft^2, flow area per fuel pin
G = m_dot/A_flow; %lbm/ft^2-s, average mass flux
v_avg = G/rho_nom; % ft/s

De = 4*A_flow/(pi*D_co); % ft, equivalent diameter

Re = G*De/mu_nom; 

% inlet enthalpy
h_in = water.h_pT(P_nom,T_c);

% bulk coolant enthalpy as a function of core height
h_fun = @(z) h_in + (q_p_peak/m_dot)*(H/pi)*(sin(pi*z/H)+1); % BTU/lbm

core_z = @(z) z - H/2; % so core position can go from 0 to H rather than -H/2 to H/2

q_p = @(z) q_p_peak*cos(pi*z/H); % linear heat rate as a function of position -H/2 < z < H/2
q_pp = @(z) q_p(z)./(pi*D_co); % BTU/s-ft^2, heat flux

Nz = 100;
Z = linspace(0,H,Nz);

h_b = h_fun(core_z(Z)); % bulk coolant enthalpy
T_b = nan(1,Nz);
for z = 1:Nz
    T_b(z) = water.T_ph(P_nom,h_b(z)); % bulk coolant temperature
end

% plot bulk enthalpy
figure(1)
plot(Z,h_b,'-b',...
    Z,ones(1,Nz)*h_sat,'-r','linewidth',2)
title('Enthalpy vs Core Height','fontsize',16,'fontweight','bold');
xlabel('Axial Height in Core [ft]','fontsize',14,'fontweight','bold');
ylabel('Specific Enthalpy [BTU/lbm]','fontsize',14,'fontweight','bold');
grid on
legend('bulk','sat','location','best');
set(gca,'fontsize',12,'fontweight','bold');

% plot bulk temperature
figure(2)
plot(Z,T_b,'-b',...
    Z,ones(1,Nz)*T_sat,'-r','linewidth',2)
title('Temperature vs Core Height','fontsize',16,'fontweight','bold');
xlabel('Axial Height in Core [ft]','fontsize',14,'fontweight','bold');
ylabel('Temperature [^\circ F]','fontsize',14,'fontweight','bold');
grid on
legend('bulk','sat','location','best');
set(gca,'fontsize',12,'fontweight','bold');

% compute and plot equilibrium quality
x_eq = (h_b - h_sat)/h_fg; % equilibrium quality

figure(3)
plot(Z,x_eq,'-b','linewidth',2)
title('Equilibrium Quality vs Core Height','fontsize',16,'fontweight','bold');
xlabel('Axial Height in Core [ft]','fontsize',14,'fontweight','bold');
ylabel('x_{eq}','fontsize',14,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');

%% correlations

Nu_DB = @(Re,Pr) 0.023*(Re.^0.8).*(Pr.^0.4); 

Nu = Nu_DB(Re,Pr_nom);

h_nom_DB = Nu*k_nom/De; % BTU/ft^2-sec-R, nominal convective heat transfer coefficient

% Presser correction factor
C_presser = @(P_D) 0.9217 + 0.1478*P_D - 0.1130*exp(-7*(P_D-1));
h_nom_FA = C_presser(Pitch/D_co)*h_nom_DB;

% Wall temperature using Presser correlation
T_wall = T_b + q_p(core_z(Z))/(2*pi*h_nom_FA*(D_co/2));

figure(4)
plot(Z,T_wall,'-b',...
    Z,ones(1,Nz)*T_sat,'-r','linewidth',2)
title('Wall Temperature vs Core Height','fontsize',16,'fontweight','bold');
xlabel('Axial Height in Core [ft]','fontsize',14,'fontweight','bold');
ylabel('Wall Temperature [^\circ F]','fontsize',14,'fontweight','bold');
grid on
legend('T_{wall}','T_{sat}','location','best');
set(gca,'fontsize',12,'fontweight','bold');

figure(5)
plot(Z,q_p(core_z(Z))/(pi*D_co),'-b','linewidth',2)
title('Heat Flux vs Core Height','fontsize',16,'fontweight','bold');
xlabel('Axial Height in Core [ft]','fontsize',14,'fontweight','bold');
ylabel('Heat Flux [BTU/s-ft^2]','fontsize',14,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');

% try the Bernath correlation for h_NB

De_prime = De + D_co/pi;
h_NB = 10890*(De/De_prime) + 48*v_avg/(De^0.6);%BTU/hr-ft^2-F
h_NB = h_NB*(1/3600); %BTU/s-ft^2-F

% compare h_NB to h_nom_DB
fprintf('h_NB is equal to %5.2f percent of h_nom_DB \n',...
    (h_NB/h_nom_DB)*100);

%% Use Jens-Lottes Correlation to estimate Twall where NB will start

T_wall_JL = @(z) T_sat + 1.897*exp(-P_nom/900).*(q_pp(core_z(z))*3600).^(0.25);
T_wall_Thom = @(z) T_sat + 0.072*exp(-P_nom/1260).*(q_pp(core_z(z))*3600).^(0.5);


figure(6)
plot(Z,T_wall,'-b',...
    Z,T_wall_JL(Z),'-r',...
    Z,T_wall_Thom(Z),'-k',...
    Z,ones(1,Nz)*T_sat,'-g','linewidth',2)
title('Wall Temperature vs Core Height','fontsize',16,'fontweight','bold');
xlabel('Axial Height in Core [ft]','fontsize',14,'fontweight','bold');
ylabel('Wall Temperature [^\circ F]','fontsize',14,'fontweight','bold');
grid on
legend('T_{wall}','Onset NB (JL)','Onset NB (Thom)','sat','location','best');
set(gca,'fontsize',12,'fontweight','bold');


figure(7)
plot(Z,T_b,'-b',...
    Z,T_wall,'-r',...
    Z,q_p(core_z(Z))/(pi*D_co),'-k','linewidth',3);
grid on;
title('Heat Flux and Temperature Profile',...
    'fontsize',16,'fontweight','bold')
legend('T_{bulk}','T_{wall}','q^{\prime \prime}')
xlabel('Axial Height in Core [ft]','fontsize',14,'fontweight','bold')
ylabel('Temperature [^\circ F] / Heat Flux [BTU/s-ft^2]',...
    'fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');

