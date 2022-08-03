%% NuScale Natural Circulation Demo

% Purpose: Illustrate the relationship between mass flow rate and reactor
% power for a PWR utilizing natural circulation.

%% Establish Environment

clear
clc
close 'all'

%% Compute Mass Flow Rate vs. Power

Q = 160e3; % MW, Design thermal power
Qmin = 1.6e3; % 1.6 MW
Qmax = 160e3; % 160 MW
nQ = 20;
Qspace = linspace(Qmin,Qmax,nQ);
m_dot = nan(nQ,1);

for i = 1:nQ
    m_dot(i) = NuScale_model(Qspace(i));
end

%% Make a nice plot
figure(1)
plot(Qspace/1e3,m_dot/1e6,'-s','linewidth',3);
grid on
title('Natural Circulation Mass Flow Rate vs. Power','fontsize',16,...
    'fontweight','bold')
xlabel('Power (MW_{th})','fontsize',14,'fontweight','bold');
ylabel('Mass Flow Rate (lb_m/hr x 10^{-6})','fontsize',14,...
    'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');

function m_dot = NuScale_model(Q_th)
%function m_dot = Nuscale_model(Q_th) - a simple natural circulation model
%of the NuScale PWR.  
% Input:
% Q_th - thermal power in kW
%
% Output:
% 
% m_dot - mass flow rate (lbm/hr)

%%
fBal = @(m) abs(NuScale_Hydraulics(m) - NuScale_Thermal(m,Q_th));

initial_guess = 5e6;
lb = 1e5;
ub = 7e6;
A = []; 
b = [];
Aeq = [];
beq = [];
options = optimoptions('fmincon','FunctionTolerance',1e-6,...
    'ConstraintTolerance',1e-9,'StepTolerance',1e-15,...
    'Display','off',...
    'MaxFunctionEvaluations',1500);

m_dot = fmincon(fBal,initial_guess,A,b,Aeq,beq,lb,ub,[],options);

end

function  dP  = NuScale_Thermal( m_dot, Q_th )
%NuScale_Thermal Simplified thermal driving head model of the NuScale PWR.
%thermal driving head (expressed as dP) given as a function of mass flow
%rate and thermal power
%
% Input:
% m_dot - mass flow rate lbm/hr
% Q_th - thermal power kW

%% check the environment and instantiate fluid object
EasyProp_path = ' '; %<- Path if EasyProp.py is in your current directory
if count(py.sys.path,EasyProp_path) == 0  % <-- see if desired directory is on path
    insert(py.sys.path,int32(0),EasyProp_path); %<-- if not; add it.
end

water = py.EasyProp.simpleFluid('Water','USCS');

%% Model

P_nom = 1850; % psia, given, nominal pressure
T_ave = 543; % F, given, core average temperature

% energy balance between core and coolant
dH = (Q_th*3.41e3)/m_dot; % BTU/lbm

% coolant properties
h_nom = water.h_pT(P_nom,T_ave); % BTU/lbm
h_h = h_nom+dH/2; % BTU/lbm, hot-leg enthalpy
h_c = h_nom-dH/2; % BTU/lbm, cold-leg enthalpy
rho_h = 1/water.v_ph(P_nom,h_h); % lbm/ft^3
rho_c = 1/water.v_ph(P_nom,h_c); % lbm/ft^3

% Reactor geometry data
L_downcomer = 46.0; % ft, given
L_lower_riser = 9.4; % ft, given
L_core = 7.9; % ft, given
% Estimated height differenec between thermal source/sink
dL_thermal = L_downcomer - L_lower_riser - L_core; % ft, estimated thermal elevation change

g = 32.2; % ft/sec^2
g_c = 32.2; % lbm/lbf * ft/sec^2

dP = (rho_c - rho_h)*(g/g_c)*dL_thermal; % lbf/ft^2 
end

function dP  = NuScale_Hydraulics( m_dot )
%NuScale_Hydraulics: provide a hydraulic model of the NuScale Reactor.
% input:
% m_dot - mass flow rate, lbm/hr
% output:
%
%  dP - differential pressure in psia
%
% For the hydraulic model, we will assume that fluid density and
% viscosity are at their nominal conditions
%% check the environment and instantiate fluid object
EasyProp_path = ' '; %<- Path if EasyProp.py is in your current directory
if count(py.sys.path,EasyProp_path) == 0  % <-- see if desired directory is on path
    insert(py.sys.path,int32(0),EasyProp_path); %<-- if not; add it.
end

water = py.EasyProp.simpleFluid('Water','USCS');

%% Model

P_nom = 1850; % psia, given, nominal pressure
T_ave = 543; % F, given, core average temperature
rho_nom = 1./water.v_pT(P_nom,T_ave);

% estimate core flow area:
N_assy = 37;
N_pins_per_assy = 17*17; 
P_assy = 8.466/12; % ft, assembly pitch
D_pin = 0.374/12; % ft, pin outside diameter 

% total core flow area taken as 
%totality of flow areas for all pins and
% assemblies.
A_flow_core = N_assy*(P_assy^2 - (pi/4)*(D_pin^2)*N_pins_per_assy); %ft^2

% velocity of flow in core
v_core = (m_dot/3600)/(A_flow_core*rho_nom); % ft/sec

hyd_term = 31; % (fL/D + sum of K), unitless, assume constant.

g_c = 32.2; % lbm/lbf*ft/sec^2

dP = hyd_term*(v_core^2)*rho_nom/(2*g_c); %<- lbf/ft^2


end