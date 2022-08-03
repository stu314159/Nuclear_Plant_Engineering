%% AP1000 MDNBR Demo
clear
clc
close 'all'


H = 168/12; % ft, height of active fuel region
Nz = 100;
Z = linspace(0,H,Nz);

%function [qpp,qpp_DNB] = AP1000_model(P_fact,T_in_adj,m_dot_fact,qpp_fact)

%% nominal conditions
[qpp_nom,qpp_DNB_nom] = AP1000_model(1.0,0,1.0,1.0);
figure(1)
plot(Z,qpp_nom,'--b',...
    Z,qpp_DNB_nom,'-r','linewidth',2)
title('Nominal Conditions','fontsize',14,'fontweight','bold');
xlabel('Z [ft]','fontsize',12,'fontweight','bold');
ylabel('Heat Flux [BTU/hr-ft^2]','fontsize',12,'fontweight','bold');
set(gca,'fontsize',10,'fontweight','bold');
grid on;
legend('Heat Flux','Heat Flux (DNB)')

%% Reduced coolant mass flow rate
[qpp_low_flow,qpp_DNB_low_flow] = AP1000_model(1.0,0,0.9,1.0);
figure(2)
plot(Z,qpp_DNB_nom,'-r',...
    Z,qpp_DNB_low_flow,'--k',...
    Z,qpp_low_flow,'-k','linewidth',2)
title('Change in DNB with Low Flow','fontsize',14,'fontweight','bold');
xlabel('Z [ft]','fontsize',12,'fontweight','bold');
ylabel('Heat Flux [BTU/hr-ft^2]','fontsize',12,'fontweight','bold');
set(gca,'fontsize',10,'fontweight','bold');
grid on;

%% increased heat flux
[qpp_higher_flux,qpp_DNB_higher_flux] = AP1000_model(1.0,0,1.0,1.1);
figure(3)
plot(Z,qpp_DNB_nom,'-r',...
    Z,qpp_DNB_higher_flux,'--k',...
    Z,qpp_higher_flux,'-k',...
    Z,qpp_nom,'-b','linewidth',2)
title('Change in DNB with Linear Power','fontsize',14,'fontweight','bold');
xlabel('Z [ft]','fontsize',12,'fontweight','bold');
ylabel('Heat Flux [BTU/hr-ft^2]','fontsize',12,'fontweight','bold');
set(gca,'fontsize',10,'fontweight','bold');
grid on;

%% reduced pressure
[qpp_low_pressure,qpp_DNB_low_pressure] = AP1000_model(0.94,0,1.0,1.0);
figure(4)
plot(Z,qpp_DNB_nom,'-r',...
    Z,qpp_DNB_low_pressure,'--k',...
    Z,qpp_low_pressure,'-k','linewidth',2)
title('Change in DNB with Nominal Pressure','fontsize',14,'fontweight','bold');
xlabel('Z [ft]','fontsize',12,'fontweight','bold');
ylabel('Heat Flux [BTU/hr-ft^2]','fontsize',12,'fontweight','bold');
set(gca,'fontsize',10,'fontweight','bold');
grid on;

%% increased Tc
[qpp_high_Tin,qpp_DNB_high_Tin] = AP1000_model(1.0,10,1.0,1.0);
figure(5)
plot(Z,qpp_DNB_nom,'-r',...
    Z,qpp_DNB_high_Tin,'--k',...
    Z,qpp_high_Tin,'-k','linewidth',2)
title('Change in DNB with Inlet Temperature','fontsize',14,'fontweight','bold');
xlabel('Z [ft]','fontsize',12,'fontweight','bold');
ylabel('Heat Flux [BTU/hr-ft^2]','fontsize',12,'fontweight','bold');
set(gca,'fontsize',10,'fontweight','bold');
grid on;



function [qpp,qpp_DNB] = AP1000_model(P_fact,T_in_adj,m_dot_fact,qpp_fact)
%% Put EasyProp on Python path
% put location of EasyProp.py module on the python search path
if count(py.sys.path,' ') == 0  % <-- see if desired directory is on path
    insert(py.sys.path,int32(0),' '); %<-- if not; add it.
end

water = py.EasyProp.simpleFluid('Water','USCS');


%% AP1000 Given Data

P_nom = 2200; % psia, nominal system pressure
P_nom = P_nom*P_fact;
mu_nom = 0.23*(1/3600); % lbm/ft-s

h_sat = water.hL_p(P_nom); % BTU/lbm, enthalpy of saturated liquid at nominal pressure
h_fg = water.hV_p(P_nom) - water.hL_p(P_nom); % BTU/lbm
m_dot_total = 113.5e6*(1/3600); % lbm/s, flow core flow rate
N_rods = 157*17*17; 
D_co = 0.374/12; % ft, clad outside diameter
D_fuel = 0.3225/12; % ft, diameter of fuel (consistent with clad dimensions and diametral gap)
Pitch = 0.496/12; % ft, fuel rod pitch
H = 168/12; % ft, height of active fuel region
q_p_peak = (14.9/1.0)*3412*(1/3600); % kW/ft --> BTU/s-ft, axial peak linear heat rate
q_p_peak = q_p_peak*qpp_fact;
T_c = 535+T_in_adj; % F, nominal core inlet temperature

%% calculated results
m_dot = m_dot_total/N_rods; % lbm/s, flow rate per fuel rod
m_dot = m_dot_fact*m_dot;
A_flow = Pitch^2 - (pi/4)*D_fuel^2; % ft^2, flow area per fuel pin
G = m_dot/A_flow; %lbm/ft^2-s, average mass flux
De = 4*A_flow/(pi*D_co); % ft, equivalent diameter
% inlet enthalpy
h_in = water.h_pT(P_nom,T_c);
% bulk coolant enthalpy as a function of core height
h_fun = @(z) h_in + (q_p_peak/m_dot)*(H/pi)*(sin(pi*z/H)+1); % BTU/lbm
core_z = @(z) z - H/2; % so core position can go from 0 to H rather than -H/2 to H/2
q_p = @(z) q_p_peak*cos(pi*z/H); % linear heat rate as a function of position -H/2 < z < H/2

Nz = 100;
Z = linspace(0,H,Nz);

h_b = h_fun(core_z(Z)); % bulk coolant enthalpy
T_b = nan(1,Nz);
for z = 1:Nz
    T_b(z) = water.T_ph(P_nom,h_b(z)); % bulk coolant temperature
end

% compute and plot equilibrium quality
x_eq = (h_b - h_sat)/h_fg; % equilibrium quality


K_Tong = 1.76 - 7.433*x_eq + 12.222*(x_eq.^2);
qpp_DNB_Tong = K_Tong.*(G^0.4)*(mu_nom^0.6)*h_fg/(De^0.6);

qpp = q_p(core_z(Z))/(pi*D_co);
qpp_DNB = qpp_DNB_Tong;

end