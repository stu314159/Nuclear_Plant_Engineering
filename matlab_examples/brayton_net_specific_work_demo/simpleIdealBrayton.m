function [w_net,eta_th] = simpleIdealBrayton(gas,Pmin,Tmin,Tmax,rp)

%% Initialize state point arrays for Brayton Cycle
numSP = 4;
h = nan(numSP,1);
hs = nan(numSP,1);
s = nan(numSP,1);
T = nan(numSP,1);
P = nan(numSP,1);

%% Problem Parameters
P_in = Pmin;
r_p = rp; % given



eta_c = 1;
eta_t = 1;

%% Compute State Point Property Data
P(1) = P_in;
T(1) = Tmin;
h(1) = gas.h_pT(P(1),T(1));
s(1) = gas.s_pT(P(1),T(1));

P(2) = P(1)*r_p;
hs(2) = gas.h_ps(P(2),s(1));
h(2) = h(1) - (h(1) - hs(2))./eta_c;
T(2) = gas.T_ph(P(2),h(2));
s(2) = gas.s_ph(P(2),h(2));

if T(2) > Tmax %<-- if this happens, stop analyzing the cycle.
    w_net = 0;
    eta_th = 0;
    return
end

P(3) = P(2); % assume isobaric heat addition
T(3) = Tmax;
h(3) = gas.h_pT(P(3),T(3));
s(3) = gas.s_pT(P(3),T(3));

P(4) = (1./r_p)*P(3);
hs(4) = gas.h_ps(P(4),s(3));
h(4) = h(3) - (h(3) - hs(4))*eta_t;
%s(4) = air.s_ph(P(4),h(4));
%T(4) = air.T_ph(P(4),h(4));



% % display state point data neatly
% SP = {'1','2','3','4'};
% SP_table = table(P,T,h,s,ef,'RowName',SP);
% disp(SP_table);

%% First Law Analysis
w_c = h(1) - h(2);
w_t = h(3) - h(4);
w_net = w_c + w_t;

q_s = h(3) - h(2);
q_r = h(1) - h(4);
q_net = q_s + q_r;

if abs(w_net - q_net) > 1 % since I am not otherwise checking this.
    error('Conservation of energy conditions not met!');
end
% fprintf('Net specific work: %g BTU/lbm \n',w_net);
% fprintf('Net specific heat added: %g BTU/lbm \n',q_net);


eta_th = w_net/q_s;