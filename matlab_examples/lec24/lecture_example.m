%% Thermal Calcs
clear
clc
close 'all'

%% AP1000 Given Data

% geometry
D_co = 0.374/12; % ft, clad outside diameter
Pitch = 0.496/12; % ft, fuel rod pitch
v_avg = 15.9; % ft/s, average velocity

mu_nom = 0.23*(1/3600); % lbm/ft-s
rho_nom = 45.2; % lbm/ft^3, nominal density
k_nom = 0.32*(1/3600); % BTU/s-ft-R
c_p = 1.31; % BTU/lbm-R


heat_flux = 346104; % BTU/hr-ft^2, heat flux at location of interest.

%% calculated results

A_flow = Pitch^2 - (pi/4)*D_co^2; % ft^2, flow area per fuel pin
fprintf('Aflow = %g ft^2 \n',A_flow);
P_wetted = pi*D_co; % ft, wetted perimeter
fprintf('P_wetted = %g ft \n',P_wetted);
De = 4*A_flow/(pi*D_co); % ft, equivalent diameter.
fprintf("De = %g \n",De);

Re = rho_nom*v_avg*De/mu_nom;
fprintf('Re = %g \n',Re);

Pr = mu_nom * c_p / k_nom;
fprintf('Pr = %g \n',Pr);

Nu_DB = @(Re,Pr) 0.023*(Re.^0.8).*(Pr.^0.4); 

Nu = Nu_DB(Re,Pr);
fprintf('Nu via DB = %g \n',Nu);

% Presser correction factor
C_presser = @(P_D) 0.9217 + 0.1478*P_D - 0.1130*exp(-7*(P_D-1));

Nu_presser = C_presser(Pitch/D_co)*Nu;
fprintf('Nu via Presser = %g \n',Nu_presser);

% Markoczy correction factor
B = @(P_D) (4/pi)*(P_D)^2 - 1;

C_markoczy = @(P_D,Re,Pr) 1 + 0.9120*(Re.^(-0.1)).*(Pr.^(0.4)).*...
    (1 - 2.0043*exp(-B(P_D)));

Nu_markoczy = C_markoczy(Pitch/D_co,Re,Pr)*Nu;
fprintf('Nu via Markoczy = %g \n',Nu_markoczy);

h_mark = Nu_markoczy*k_nom*3600/De;
fprintf('h = %g \n',h_mark);

dT_cool = heat_flux/h_mark;
fprintf('DeltaT = %g \n',dT_cool);

