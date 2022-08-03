function [w_net,eta_th] = idealBraytonRegenIC_Reheat(gas,Pmin,Tmin,Tmax,rp)

%% Initialize state point arrays for Brayton Cycle
numSP = 11;
h = nan(numSP,1);
hs = nan(numSP,1);
s = nan(numSP,1);
T = nan(numSP,1);
P = nan(numSP,1);

%% Problem Parameters
P_in = Pmin;
%r_p = rp; % given
r_pc = sqrt(rp);
r_pt = rp.^(1/3);

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

P(5) = P(4); % assume isobaric

T(6) = Tmax;
P(6) = P(4);
h(6) = gas.h_pT(P(6),T(6));
s(6) = gas.s_pT(P(6),T(6));

% SP 7
P(7) = P(6)./r_pt;
hs(7) = gas.h_ps(P(7),s(6));
h(7) = h(6) - eta_t*(h(6) - hs(7));
s(7) = gas.s_ph(P(7),h(7));
T(7) = gas.T_ph(P(7),h(7));

% SP 8
P(8) = P(7)./r_pt;
hs(8) = gas.h_ps(P(8),s(7));
h(8) = h(7) - eta_t*(h(7) - hs(8));
T(8) = gas.T_ph(P(8),h(8));

T(9) = (5/10)*T(6);
P(9) = P(8);
h(9) = gas.h_pT(P(9),T(9));
s(9) = gas.s_pT(P(9),T(9));

P(10) = P(9)./r_pt;
hs(10) = gas.h_ps(P(10),s(9));
h(10) = h(9) - (h(9) - hs(10))*eta_t;
T(10) = gas.T_ph(P(10),h(10));


%eta_regen = (h(5) - h(4))/((h(10)-h(2)))
% f*(h(6) - h(5)) = (1-f)*(h(9) - h(8))

%--> unknowns: x = [ f, h(5)]
regen_1 = @(x) eta_regen*(h(10)-h(2)) - (x(2) - h(4));
regen_2 = @(x) x(1)*(h(6) - x(2)) - (1-x(1))*(h(9) - h(8));

balance = @(x) abs(regen_1(x)) + abs(regen_2(x));

% % provide inputs to fmincon to solve non-linear equations
x0 = [0.5 (h(4)+h(6))./2];
A = []; % no linear inequality constraints
b = [];
Aeq = []; % no linear equality constraints
beq = [];
lb = [0 h(4)]; % all values must be non-negative
ub = [1 h(6)]; % flow fraction limited to 1
nlcon = [];
%options = [];
options = optimoptions('fmincon','Display','none'); %<-- no output display
% %options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
% %options = optimoptions('fmincon','Display','iter');
% % default is 'interior-point'
% % other options: 'sqp', 'sqp-legacy', 'active-set',
% % 'trust-region-reflective'
[x,~,~] = fmincon(balance,x0,A,b,Aeq,beq,lb,ub,nlcon,options);

f = x(1);
h(5) = x(2);
T(5) = gas.T_ph(P(5),h(5));



% (1-f)*(h10 - h11) = (1-f)*(h(5) - h(4))
h(11) = h(10) - (h(5) - h(4));
P(11) = P(10);
T(11) = gas.T_ph(P(11),h(11));
if T(11) > T(10)
    w_net = 0;
    eta_th = 0;
    return;
end

%% Check the first law balance
w_c1 = (1-f)*(h(1) - h(2));
w_c2 = (1-f)*(h(3) - h(4));
w_t1 = (1-f)*(h(6) - h(7));
w_t2 = (1-f)*(h(7) - h(8));
w_t3 = (1-f)*(h(9) - h(10));

w_net = w_c1 + w_c2 + w_t1 + w_t2 + w_t3;

q_IC = (1-f)*(h(3) - h(2));
q_PC = (1-f)*(h(1) - h(11));
q_PB = h(6) - h(5);

q_net = q_IC + q_PC + q_PB;

q_s = q_PB;
eta_th = w_net/q_s;


