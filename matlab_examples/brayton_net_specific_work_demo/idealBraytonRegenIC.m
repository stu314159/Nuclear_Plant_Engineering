function [w_net,eta_th] = idealBraytonRegenIC(gas,Pmin,Tmin,Tmax,rp)

%% Initialize state point arrays for Brayton Cycle
numSP = 8;
h = nan(numSP,1);
hs = nan(numSP,1);
s = nan(numSP,1);
T = nan(numSP,1);
P = nan(numSP,1);

%% Problem Parameters
P_in = Pmin;
r_p = rp; % given
r_pc = sqrt(r_p);



eta_c = 1;
eta_t = 1;
eta_regen = 0.65;


%% Compute State Point Property Data
P(1) = P_in;
T(1) = Tmin;
h(1) = gas.h_pT(P(1),T(1));
s(1) = gas.s_pT(P(1),T(1));

P(2) = P(1)*r_pc;
hs(2) = gas.h_ps(P(2),s(1));
h(2) = h(1) - (h(1) - hs(2))./eta_c;
T(2) = gas.T_ph(P(2),h(2));
s(2) = gas.s_ph(P(2),h(2));

P(3) = P(2);
T(3) = Tmin;
h(3) = gas.h_pT(P(3),T(3));
s(3) = gas.s_pT(P(3),T(3));

P(4) = P(3)*r_pc;
hs(4) = gas.h_ps(P(4),s(3));
h(4) = h(3) - (h(3) - hs(4))./eta_c;
T(4) = gas.T_ph(P(4),h(4));


% skip to SP 6
T(6) = Tmax;
P(6) = P(4);
h(6) = gas.h_pT(P(6),T(6));
s(6) = gas.s_pT(P(6),T(6));

% SP 7
P(7) = P(6)./r_p;
hs(7) = gas.h_ps(P(7),s(6));
h(7) = h(6) - eta_t*(h(6)-hs(7));
T(7) = gas.T_ph(P(7),h(7));

% back to SP5
P(5) = P(4);
h(5) = h(4) + eta_regen*(h(7) - h(4));
T(5) = gas.T_ph(P(5),h(5));

if T(5) > T(7)
    w_net = 0;
    eta_th = 0;
    return;
end

% SP 8
h(8) = h(7) - (h(5) - h(4)); % conservation of energy on the regenerator




%% First Law Analysis
w_c = (h(1) - h(2)) + (h(3) - h(4));
w_t = h(6) - h(7);
w_net = w_c + w_t;


q_s = h(6) - h(5);
q_r = (h(1) - h(8))+ (h(3) - h(2));
q_net = q_s + q_r;

if abs(w_net - q_net) > 1 % since I am not otherwise checking this.
    error('Conservation of energy conditions not met!');
end
% fprintf('Net specific work: %g BTU/lbm \n',w_net);
% fprintf('Net specific heat added: %g BTU/lbm \n',q_net);


eta_th = w_net/q_s;