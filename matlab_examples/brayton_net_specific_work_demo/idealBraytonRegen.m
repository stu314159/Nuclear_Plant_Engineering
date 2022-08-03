function [w_net,eta_th] = idealBraytonRegen(gas,Pmin,Tmin,Tmax,rp)

%% Initialize state point arrays for Brayton Cycle
numSP = 6;
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
eta_regen = 0.65;

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

% skip to SP 4
T(4) = Tmax;
P(4) = P(2);
h(4) = gas.h_pT(P(4),T(4));
s(4) = gas.s_pT(P(4),T(4));

% SP 5
P(5) = P(4)./r_p;
hs(5) = gas.h_ps(P(5),s(4));
h(5) = h(4) - eta_t*(h(4)-hs(5));

% back to SP3
%P(3) = P(2);
h(3) = h(2) + eta_regen*(h(5) - h(2));

% SP 6
h(6) = h(5) - (h(3) - h(2)); % conservation of energy on the regenerator




%% First Law Analysis
w_c = h(1) - h(2);
w_t = h(4) - h(5);
w_net = w_c + w_t;

q_s = h(4) - h(3);
q_r = h(1) - h(6);
q_net = q_s + q_r;

if abs(w_net - q_net) > 1 % since I am not otherwise checking this.
    error('Conservation of energy conditions not met!');
end
% fprintf('Net specific work: %g BTU/lbm \n',w_net);
% fprintf('Net specific heat added: %g BTU/lbm \n',q_net);


eta_th = w_net/q_s;