% EasyProp from matlab for Rankine Cycle Analysis

% Example problem using EasyProp for fluid properties.
% Problem is adapted from Example 8-5 of Introduction to Thermal and Fluids Engineering,
% 1st Edition, D. Kaminski and M. Jensen, John Wiley & Sons
% 
% A Rankine cycle uses water as the working fluid.  The boiler operates at 8 MPa
% and produces saturated vapor.  Saturated liquid exits the condenser at 7.5 kPa.
% The pump has an isentropic efficiency of 82 percent, and the turbine has an 
% isentropic efficiency of 88 percent.  The steam flow rate is 2.8x10**4 kg/h.
% Cooling water enters the condenser at 20 C and leaves at 38 C.
% 
% Find:
% a) the net power produced (kW)
% b) the net heat transfer rate into the boiler (kW)
% c) cycle thermal efficiency (%)
% d) cooling water flow rate (kg/h)

%% clear the environment
clear
clc
close 'all'

%% Add current directory to the Python Path
EasyProp_path = ' '; %<- Path if EasyProp.py is in your current directory
if count(py.sys.path,EasyProp_path) == 0  % <-- see if desired directory is on path
    insert(py.sys.path,int32(0),EasyProp_path); %<-- if not; add it.
end

%% Establish working fluid and unit system
fluid = 'Water';
unitSystem = 'SI';
myFluid = py.EasyProp.simpleFluid(fluid,unitSystem);

%% Problem parameters
numSp = 4;
eta_turbine = 0.88;
eta_pump = 0.82;
m_dot_steam = 2.8e4; % kg/h
Pmax = 8000; % kPa
Pmin = 7.5; % kPa
Tcool_in = 20; % C
Tcool_out = 38; % C

%% declare fluid property arrays
h = nan(numSp,1);
h_s = nan(numSp,1);
s = nan(numSp,1);
s_s = nan(numSp,1);
P = nan(numSp,1);
T = nan(numSp,1);
T_s = nan(numSp,1);

%% calculate properties at state points 
% state point 1 - condenser outlet condensate pump inlet
P(1) = Pmin;
h(1) = myFluid.hL_p(P(1));
s(1) = myFluid.sL_p(P(1));

% compression from state point 1 to state point 2
P(2) = Pmax;
s_s(2) = s(1);
h_s(2) = myFluid.h_ps(P(2),s_s(2));
h(2) = h(1) - (h(1) - h_s(2))./eta_pump;

w_pump = h(1) - h(2);

% isobaric heat addition in the boiler from state point 2 to state point 3
P(3) = P(2);
h(3) = myFluid.hV_p(P(3));
s(3) = myFluid.sV_p(P(3));

q_s = h(3) - h(2);

% expansion in turbine from state point 3 to state point 4
P(4) = P(1);
s_s(4) = s(3);
h_s(4) = myFluid.h_ps(P(4),s_s(4));
h(4) = h(3) - eta_turbine*(h(3) - h_s(4));

w_turbine = h(3) - h(4);

% isobaric heat rejection in the condenser from sp 4 to sp 1
q_r = h(1) - h(4);

%% First Law Check
w_net = w_pump + w_turbine;
q_net = q_s + q_r;

fprintf('Net work = %g kJ/kg \n',w_net);
fprintf('Net heat = %g kJ/kg \n',q_net);

%% Question a) Net power (kW)
net_power = w_net*m_dot_steam*(1/3600); % kW

fprintf('Net power = %g kW\n',net_power);

%% Question b) Net heat transfer rate into the boiler
net_Qdot_supplied = q_s*m_dot_steam*(1/3600); % kW

fprintf('Net heat transfer into the boiler = %g kW \n',net_Qdot_supplied);

%% Question c) cycle thermal efficiency
eta_th = w_net/q_s;

fprintf('Thermal efficieency = %4.1f percent \n',eta_th*100);

%% Question d) condenser cooling water flow rate

Po = 101.325; % kPa - assume Standard Pressure
h_cool_in = myFluid.h_pT(Po,Tcool_in);
h_cool_out = myFluid.h_pT(Po,Tcool_out);

m_dot_cool = m_dot_steam*(h(4) - h(1))/(h_cool_out - h_cool_in); % kg/h

fprintf('Condenser cooling water flow rate = %g kg/h \n',m_dot_cool);