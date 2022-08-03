%% Simple Brayton Cycle
clear
clc
close 'all'

%% Put EasyProp on Python path
% put location of EasyProp.py module on the python search path
if count(py.sys.path,' ') == 0  % <-- see if desired directory is on path
    insert(py.sys.path,int32(0),' '); %<-- if not; add it.
end

%% Initialize Fluid Property object
fluid = 'Air';
units = 'USCS';
air = py.EasyProp.simpleFluid(fluid,units);

%% Initialize state point arrays for Brayton Cycle
numSP = 4;
h = nan(numSP,1);
hs = nan(numSP,1);
s = nan(numSP,1);
ss = nan(numSP,1);
T = nan(numSP,1);
P = nan(numSP,1);

%% Problem Parameters
Pmin = 14.7; % psia, atmospheric pressure
Tin = 60; % F, inlet air temperature
r_p = 14;
r_e = 1/r_p;
T_turb_inlet = 1600; % F

eta_c = 0.87;
eta_t = 0.9;

%% Compute State Point Property Data
% state point 1
P(1) = Pmin;
T(1) = Tin;
h(1) = air.h_pT(P(1),T(1));
s(1) = air.s_pT(P(1),T(1));

% state point 2
P(2) = P(1)*r_p;
ss(2) = s(1);
hs(2) = air.h_ps(P(2),ss(2));
h(2) = h(1) - (h(1)-hs(2))./eta_c;
T(2) = air.T_ph(P(2),h(2));
s(2) = air.s_ph(P(2),h(2));

% state point 3
P(3) = P(2); %isobaric heat addition
T(3) = T_turb_inlet; % F, given
h(3) = air.h_pT(P(3),T(3));
s(3) = air.s_pT(P(3),T(3));

% state point 4
P(4) = P(3)*r_e;
ss(4) = s(3); 
hs(4) = air.h_ps(P(4),ss(4));
h(4) = h(3)-(h(3)-hs(4))*eta_t;
T(4) = air.T_ph(P(4),h(4));
s(4) = air.s_ph(P(4),h(4));

% display state point data neatly
SP = {'1','2','3','4'};
SP_table = table(P,T,h,s,'RowName',SP);
disp(SP_table);

%% First Law Analysis
w_c = h(1) - h(2);
w_t = h(3) - h(4);
w_net = w_c + w_t;

q_s = h(3) - h(2);
q_r = h(1) - h(4);
q_net = q_s + q_r;

assert(abs(w_net - q_net)<1,'Conservation of energy condition not met!');

fprintf('Net specific work: %g BTU/lbm \n',w_net);
fprintf('Net specific heat added: %g BTU/lbm \n',q_net);


eta_th = w_net/q_s;

fprintf('Thermal Efficiency: %5.2f percent \n',eta_th*100);