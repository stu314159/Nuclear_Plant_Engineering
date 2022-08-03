%% Lecture 22 thermal conductivity
clear
clc
close 'all'
%%
% how does k vary with temperature?
Bu = 0; gad = 0; eta_td = 0.95; fuel = 'UO2';
k_vs_T = @(T) k_frapcon(T,Bu,gad,eta_td,fuel);
figure(1)
fplot(k_vs_T,[200+273,2500+273],'linewidth',3);
axis([500 2500 0 0.065]);
title('Thermal conductivity vs Temperature','fontsize',16,...
    'fontweight','bold');
grid on
xlabel('Temperature [^{\circ}K]','fontsize',14,'fontweight','bold');
ylabel('Thermal Conductivity [W/cm-K]','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');
%%
% get average for unirradiated fuel between 200 and 1000 C
k_avg = (1/(1000-200))*integral(@(T) k_vs_T(T),200+273,1000+273);
fprintf('Average thermal conductivity: %5.4f W/cm-K \n',k_avg);
% verify against number given in Table 8.1A on page 361
%%
% compare with MOX
k_vs_T_MOX = @(T) k_frapcon(T,Bu,gad,eta_td,'MOX');
figure(2)
fplot(k_vs_T,[200+273,2500+273],'-.b','linewidth',3);
hold on
fplot(k_vs_T_MOX,[200+273,2500+273],':r','linewidth',3);
hold off
title('Thermal conductivity vs Temperature','fontsize',16,...
    'fontweight','bold');
grid on
xlabel('Temperature [^{\circ}K]','fontsize',14,'fontweight','bold');
ylabel('Thermal Conductivity [W/cm-K]','fontsize',14,'fontweight','bold');
legend('UO_2','MOX','location','best');
set(gca,'fontsize',12,'fontweight','bold');
k_avg_mox = (1/(1000-200))*integral(@(T) k_vs_T_MOX(T),200+273,1000+273);
fprintf('Average thermal conductivity for MOX: %5.4f W/cm-K \n',k_avg_mox);
%%
% thermal conductivity vs BU
k_vs_TandBu = @(T,BU) k_frapcon(T,BU,gad,eta_td,fuel);
k1 = @(T) k_vs_TandBu(T,0);
k2 = @(T) k_vs_TandBu(T,25);
k3 = @(T) k_vs_TandBu(T,55);
figure(3)
fplot(k1,[473,2773],'-g','linewidth',3);
hold on
fplot(k2,[473,2773],'--b','linewidth',3);
fplot(k3,[473,2773],':r','linewidth',3);
hold off
xlabel('Temperature [^{\circ}K]','fontsize',14,'fontweight','bold');
ylabel('Thermal Conductivity [W/cm-K]','fontsize',14,'fontweight','bold');
title('Thermal conductivity vs Temperature','fontsize',16,...
    'fontweight','bold');
xlabel('Temperature [^{\circ}K]','fontsize',14,'fontweight','bold');
ylabel('Thermal Conductivity [W/cm-K]','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');
legend('Fresh','25 GWd/MT','55 GWd/MT','location','best');

%%
% thermal conductivity vs density
k_vs_TandEta = @(T,eta_td) k_frapcon(T,Bu,gad,eta_td,fuel);
k1 = @(T) k_vs_TandEta(T,0.92);
k2 = @(T) k_vs_TandEta(T,0.95);
k3 = @(T) k_vs_TandEta(T,0.97);

figure(4)
fplot(k1,[473,2773],'-g','linewidth',3);
hold on
fplot(k2,[473,2773],'--b','linewidth',3);
fplot(k3,[473,2773],':r','linewidth',3);
hold off
xlabel('Temperature [^{\circ}K]','fontsize',14,'fontweight','bold');
ylabel('Thermal Conductivity [W/cm-K]','fontsize',14,'fontweight','bold');
title('Thermal conductivity vs Temperature','fontsize',16,...
    'fontweight','bold');
xlabel('Temperature [^{\circ}K]','fontsize',14,'fontweight','bold');
ylabel('Thermal Conductivity [W/cm-K]','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');
legend('0.92','0.95','0.97','location','best');

%% Local functions
function k = k_frapcon(T_in,Bu_in,gad_in,eta_td_in,fuel_in)

% FRAPCON model
% accounts for 
% temperature
% burnup
% gadolinia
% density other than 0.95% theoretical density is corrected with an
% additional factor.

% range of applicability
% 300K < T < 3000 K
% Bu < 62 GWd/Mt
% Gad < 10 w.t. %
% 92% < TD < 97%

a_frapcon = 115.99; % cm-K/W
f_frapcon = @(Bu) 0.187*Bu; %cm-K/W
g_frapcon = @(Bu) 3.8*(Bu.^0.28);% cm-K/W
Q_frapcon = 6380; %K
h_frapcon = @(T) 1./(1 + 396*exp(-Q_frapcon./T)); % unitless


% set constants from Table 8.4 based on fuel type.
fuel_type = fuel_in; % fuel type = 'UO2' or 'MOX'
switch fuel_type
    case 'UO2'
        A_frap = 4.52; %cm-K/W
        B_frap = 2.46e-2; %cm/W
        E_frap = 3.5e7; %W-K/cm
        F_frap = 16361; %K
       
    case 'MOX'
        O_M = 2.0; % set for MOX fuel - default 2.0
        x_frap = 2.0 - O_M;
        A_frap = 285*x_frap + 3.5; %cm-K/W
        B_frap = (2.86 - 7.15*x_frap)*1e-2; %cm/W
        E_frap = 1.5e7; %W-K/cm
        F_frap = 13520; %K
    otherwise
        % raise an error
        error('Unknown fuel type for FRAPCON model\n');
end

% equation 8.22a
k_0p95_frap = @(T,gad,Bu) 1./(A_frap + B_frap*T+a_frapcon*gad+f_frapcon(Bu)+...
    (1 - 0.9*exp(-0.04*Bu)).*g_frapcon(Bu).*h_frapcon(T))+...
    E_frap./(T.^2).*exp(-F_frap./T);
% (be careful with syntax to allow any argument to be a vector)

% correction factor for theoretical density
eta_frap = @(den) 1.0789*den/(1+0.5*(1-den));

% return the resulting thermal conductivity
k = k_0p95_frap(T_in,gad_in,Bu_in).*eta_frap(eta_td_in);
end

