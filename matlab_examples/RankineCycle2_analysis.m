%RankineCycle2_analysis.m
%Use of EasyProp from MATLAB environment
%% Problem Description
%
%  Consider the Rankine cycle depicted in the figure below:
%
% 
% <<RankineCycle2_Schematic.png>>
% 
% 
%  The following fluid conditions are given for each numbered state point:
% 
% # state point 1: P = 1.5 psia, x = 0.0
% # state point 2: P = 164 psia, MCP efficiency = 84 percent
% # state point 3: P = 164 psia, x = 0.0
% # state point 4: P = 838 psia, MFP efficiency = 89 percent
% # state point 5: P = 838 psia, x = 1.0
% # state point 6: P = 164 psia, HP Turbine efficiency = 94%
% # state point 7: P = 164 psia, x = 1.0
% # state point 8: P = 164 psia, T = 490 F
% # state point 9: P = 1.5 psia, LP Turbine efficiency = 92%
% # state point 10: P = 838 psia, x = 0.0
% # state point 11: P = 164 psia, isenthalpic expansion in trap T1
% # state point 12: P = 164 psia, x = 0.0, liquid drain from moisture
% separator
% # state point 13: P = 2190 psia, T = 610 F
% # state point 14: P = 2170 psia, T = 535 F
% # state point 15: P = 14.7 psia, T = 91 F
% # state point 16: P = 14.7 psia, T = 114.9 F
%
% 
% 
%  Assume primary coolant flow rate to be 113.5 lbm/hr.  Flow fraction 
%  f1 and f2 represent a fraction of the total steam flow 
%  _at the respective extraction point_ and is constrained to be between 0
%  and 1.
%
%  Calculate the following:
% 
% # The enthalpy at all state points and the flow fraction f1 and f2.
% # The net specific work of all turbines and pumps (BTU/lbm)
% # The net specific heat supplied at the steam generator and condenser
% (BTU/lbm)
% # Cycle thermal efficiency
% # Flow rate of steam from the Steam Generator and the rate of cooling
% water flow through the condenser.
% 
% 


%% clear the environment
clear
clc
close 'all'

%% Put EasyProp on Python path
% put location of EasyProp.py module on the python search path
% here, I assume that EasyProp.py is up one level in the file system
if count(py.sys.path,'') == 0  % <-- see if desired directory is on path
    insert(py.sys.path,int32(0),''); %<-- if not; add it.
end

%% Establish working fluid and unit system
fluid = 'Water';
unitSystem = 'USCS';
myFluid = py.EasyProp.simpleFluid(fluid,unitSystem);

%% State given parameters and initialize data arrays

numStatePoints = 16;%<-- so state point arrays can be treated as 1-based.

% some given parameters
eta_mcp = 0.84; % main condensate pump isentropic efficiency
eta_mfp = 0.89; % main feed pump isentropic efficiency
eta_hpt = 0.94; % high pressure turbine isentropic efficiency
eta_lpt = 0.92; % low pressure turbine isentropic efficiency

m_dot_p = 113.5e6; % lbm/hr, primary coolant system flow rate

P = nan(numStatePoints,1);
T = nan(numStatePoints,1);
T_S = nan(numStatePoints,1);
h = nan(numStatePoints,1);
h_s = nan(numStatePoints,1);
s = nan(numStatePoints,1);
s_s = nan(numStatePoints,1);
x = nan(numStatePoints,1);

%% Compute State Point Properties

% State point 1: condenser outlet, MCP inlet
P(1) = 1.5; % psia, given
h(1) = myFluid.hL_p(P(1)); % condenser outlet is saturated.
s(1) = myFluid.sL_p(P(1));

% compression from sp1 to 2
P(2) = 164; % psia, given
s_s(2) = s(1);
h_s(2) = myFluid.h_ps(P(2),s_s(2));
h(2) = h(1) - (h(1) - h_s(2))/eta_mcp;

% state point 3: OFWH outlet
P(3) = P(2); % isobaric mixing in the OFWH
h(3) = myFluid.hL_p(P(3));
s(3) = myFluid.sL_p(P(3));

% compression from state point 3 to 4
P(4) = 838; % psia, given
s_s(4) = s(3);
h_s(4) = myFluid.h_ps(P(4),s_s(4));
h(4) = h(3) - (h(3) - h_s(4))/eta_mfp;

% isobaric heat addition in the SG
P(5) = P(4);
h(5) = myFluid.hV_p(P(5));
s(5) = myFluid.sV_p(P(5));

% expansion from sp 5 to 6 in the HP Turbine
P(6) = 164; % psia, given
s_s(6) = s(5);
h_s(6) = myFluid.h_ps(P(6),s_s(6));
h(6) = h(5) - eta_hpt*(h(5) - h_s(6));
x(6) = myFluid.x_ph(P(6),h(6)); % will be needed for energy balance

% isobaric moisture separation from sp 6 to 7
P(7) = P(6);
h(7) = myFluid.hV_p(P(7));

% isobaric heat exchange in the re-heater
P(8) = P(7);
T(8) = 490; % F, given
h(8) = myFluid.h_pT(P(8),T(8));
s(8) = myFluid.s_pT(P(8),T(8));

% expansion in LP turbine
P(9) = P(1); % will be isobaric heat rejection in condenser
s_s(9) = s(8);
h_s(9) = myFluid.h_ps(P(9),s_s(9));
h(9) = h(8) - eta_lpt*(h(8) - h_s(9));

% state pont 10 reheater high-temp / pressure outlet
P(10) = P(5);
h(10) = myFluid.hL_p(P(10));

% state point 10 to 11: isenthalpic expansion in steam trap
h(11) = h(10);

% state point 12: moisture separator saturated liquid drain
P(12) = P(6);
h(12) = myFluid.hL_p(P(12));

% state point 13: Reactor core outlet / primary SG inlet
P(13) = 2190; % psia, given
T(13) = 610; % F, given
h(13) = myFluid.h_pT(P(13),T(13));

% state point 14: Reactor core inlet / primary SG outlet
P(14) = 2170; % psia, given
T(14) = 535; % F, given
h(14) = myFluid.h_pT(P(14),T(14));

% state point 15: condenser cooling water inlet
P(15) = 14.7; % psia, given
T(15) = 91; % F, given
h(15) = myFluid.h_pT(P(15),T(15));

% state point 16: condenser cooling water outlet
P(16) = 14.7; % psia, given
T(16) = 114.9; % F, given
h(16) = myFluid.h_pT(P(16),T(16));

fprintf('Specific enthalpy at all state points: \n\n');
StatePoint = (1:numStatePoints)';
enthalpyTable = table(StatePoint,h);
disp(enthalpyTable);



%% Heat Balance - find the flow fractions
RH_heatBalance = @(f) (f(1)*h(5) + (1-f(1))*x(6)*(1-f(2))*h(7)) - ...
    ( f(1)*h(10) + (1-f(1))*x(6)*(1-f(2))*h(8));

OFWH_heatBalance = @(f) ... 
    ((1-f(1))*(1-f(2))*x(6)*h(2) + f(1)*h(11) + (1-f(1))*x(6)*f(2)*h(7) + ...
    (1-f(1))*(1-x(6))*h(12)) - h(3);

totalFunctional = @(f) abs(RH_heatBalance(f))+...
    abs(OFWH_heatBalance(f));

% the strategy is to minimize the total functional.  The minimum value is
% when they are both equal to zero.
initialGuess = [0.1,0.1];
f = fminsearch(totalFunctional,initialGuess);

%% Calculate Specific Work and Heat Quantities
w_hp = (1-f(1))*(h(5) - h(6)); % BTU/lbm
w_lp = (1-f(1))*x(6)*(1-f(2))*(h(8) - h(9));
w_mcp = (1-f(1))*x(6)*(1-f(2))*(h(1)-h(2));
w_mfp = (h(3) - h(4));

w_net = w_hp + w_lp + w_mcp + w_mfp;

fprintf('Net Specific Work: %g BTU/lbm \n',w_net);

q_sg = h(5) - h(4);
q_cond = (1-f(1))*x(6)*(1-f(2))*(h(1)-h(9));

q_net = q_sg + q_cond;

fprintf('Net Specific Heat Transferred = %g BTU/lbm \n',q_net);

%% Compute Cycle Thermal Efficiency

eta_th = w_net/q_sg;

fprintf('Cycle Thermal Efficiency = %5.2f percent. \n',eta_th*100);

%% Mass flow rate of Steam Cycle and Condenser Cooling Water

m_dot_s = m_dot_p*(h(13) - h(14))/(h(5) - h(4));

m_dot_cool = m_dot_s*(1-f(1))*x(6)*(1-f(2))*(h(9)-h(1))/(h(16) - h(15));

fprintf('Steam flow rate = %g lbm/hr \n',m_dot_s);
fprintf('Condenser cooling flow rate = %g lbm/hr \n',m_dot_cool);
