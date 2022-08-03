%% SCO-2 energy conversion example
clear
clc
close 'all'

%% Add current directory to the Python Path
EasyProp_path = ' '; %<- Path if EasyProp.py is in your current directory
if count(py.sys.path,EasyProp_path) == 0  % <-- see if desired directory is on path
    insert(py.sys.path,int32(0),EasyProp_path); %<-- if not; add it.
end

%% Initialize Fluid Property object
fluid = 'CO2';
units = 'SI';
gas = py.EasyProp.EasyProp(fluid,units);

To_C = 25; % C
To = To_C + 273; % K
Po = 101.3; % kPa
ho = gas.h_pT(Po,To_C);
so = gas.s_pT(Po,To_C);


% function for specific flow exergy neglecting kinetic and potential energy
ef_fun = @(h_val,s_val) h_val - ho - To*(s_val - so);

%% Initialize State Point Data Arrays

numSP = 4;
h = nan(numSP,1);
h_s = nan(numSP,1);
s = nan(numSP,1);
s_s = nan(numSP,1);
T = nan(numSP,1);
P = nan(numSP,1);
ef = nan(numSP,1);

eta_t = 1.0;
eta_c = 1.0;

Pmax = 20000; % kPa
Pmin = 7700; % kPa

%% Statepoint Properties
P(1) = Pmin;
T(1) = 32; % C
h(1) = gas.h_pT(P(1),T(1));
s(1) = gas.s_pT(P(1),T(1));
ef(1) = ef_fun(h(1),s(1));

P(2) = Pmax;
s_s(2) = s(1);
h_s(2) = gas.h_ps(P(2),s_s(2));
h(2) = h(1) - (h(1) - h_s(2))/eta_c;
s(2) = gas.s_ph(P(2),h(2));
T(2) = gas.T_ph(P(2),h(2));
ef(2) = ef_fun(h(2),s(2));

w_comp = h(1) - h(2);

P(3) = P(2);
T(3) = 550; % C
h(3) = gas.h_pT(P(3),T(3));
s(3) = gas.s_pT(P(3),T(3));
ef(3) = ef_fun(h(3),s(3));

ef_in = ef(3) - ef(2);
fprintf('Flow exergy in: %g kJ/kg \n',ef_in);

q_s = h(3) - h(2);

P(4) = P(1);
s_s(4) = s(3);
h_s(4) = gas.h_ps(P(4),s_s(4));
h(4) = h(3) - eta_t*(h(3) - h_s(4));
s(4) = gas.s_ph(P(4),h(4));
T(4) = gas.T_ph(P(4),h(4));
ef(4) = ef_fun(h(4),s(4));

w_turb = h(3) - h(4);

q_r = h(1) - h(4);
ef_out = ef(4) - ef(1);
fprintf('Flow exergy out: %g kJ/kg \n',ef_out);

w_net = w_comp + w_turb;
q_net = q_s + q_r;

Flow_ex_balance = ef_in - ef_out - w_net;
fprintf('Specific exergy balance: %g \n',Flow_ex_balance);


fprintf('Net specific heat = %g kJ/kg \n',q_net);

eta_th = w_net/q_s;



bwr = abs(w_comp)/w_turb;

fprintf('Net specific work = %g kJ/kg \n',w_net);
fprintf('Thermal efficiency = %g percent \n',eta_th*100);
fprintf('Back work ratio = %g percent \n',bwr*100);




%%
fprintf('\n\nState point data: \n\n');

% display state point data neatly
SP = {'1','2','3','4'};
SP_table = table(P,T,h,s,ef,'RowName',SP);
disp(SP_table);