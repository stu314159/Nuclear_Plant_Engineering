%% LOCA Containment Analysis Solution.m

clear
clc
close 'all'

EasyProp_path = ' '; %<- Path if EasyProp.py is in your current directory
if count(py.sys.path,EasyProp_path) == 0  % <-- see if desired directory is on path
    insert(py.sys.path,int32(0),EasyProp_path); %<-- if not; add it.
end


%% establish air and water objects

units = 'SI';
water = py.EasyProp.simpleFluid('Water',units);
air = py.EasyProp.simpleFluid('Air',units);

%% set initial temperature, pressure, and volumes
V_primary = 354; % m^3, given
V_containment = 50970; % m^3, given
V_total = V_primary + V_containment; % m^3

Po_primary = 15500; % kPa, initial primary pressure, given
To_primary = 344.75; % C, initial primary temperature, given

Po_containment = 101; % kPa, initial containment pressure, given
To_containment = 26.85; % C, initial containment temperature, given

Q_total = 0; % kJ, total net energy exchange in/out of the containment during the transient.

%% compute the mass of air and water

rho_a_0 = 1./air.v_pT(Po_containment,To_containment); % kg/m^3
mass_air = V_containment*rho_a_0; % kg

rho_w_0 = 1./water.v_pT(Po_primary,To_primary); % kg/m^3
mass_water = V_primary*rho_w_0; % kg

%% compute specific internal energy

u_air_0 = air.u_pT(Po_containment,To_containment); %kJ/kg
u_water_0 = water.u_pT(Po_primary,To_primary); % kJ/kg

%% set up system of equations to solve for the final state.

v_air_1 = V_total/mass_air; % m^3/kg 
v_water_1 = V_total/mass_water; % m^3/kg

%% System of equations
% unknowns: T_1

Energy = @(T) mass_water*(water.u_Tv(T,v_water_1) - u_water_0) + ...
    mass_air*(air.u_Tv(T,v_air_1) - u_air_0) - Q_total;

T_1 = fzero(Energy,0.5*(To_containment + To_primary));
fprintf('Final system temperature = %g degrees C\n',T_1);

%% Determine final system pressure (no IRWST, no core make-up tanks)
P1_water = water.P_vT(v_water_1,T_1);
P1_air = air.P_vT(v_air_1,T_1);

P1_total = P1_water + P1_air;

% get final water quality
x_1 = water.x_pu(P1_water,water.u_Tv(T_1,v_water_1));
fprintf('Final water quality = %g \n',x_1);

fprintf('Partial pressure of water = %g kPa \n',P1_water);
fprintf('Partial pressure of air = %g kPa \n',P1_air);
fprintf('Total pressure = %g kPa \n',P1_total);

%% Part 2 - add IRWST and Core make-up water volume + decay heat:
fprintf('\n\n\n Add IRWST, Core make-up water, and decay heat \n\n');
V_irwst = 2070; % m^3
V_core_makeup = 141; % m^3

V_cooling = V_irwst + V_core_makeup;
To_cooling = To_containment;
Po_cooling = Po_containment;
rho_cooling = 1./water.v_pT(Po_cooling,To_cooling);
mass_cooling = V_cooling*rho_cooling;
u_cooling_0 = water.u_pT(Po_cooling,To_cooling);

V_total = V_containment + V_primary + V_cooling;

mass_water_total = mass_water + mass_cooling;
v_water_1 = V_total/mass_water_total;
v_air_1 = V_total/mass_air;

% include decay heat
Po = 3400e3; % kW_th - Reactor Power
Ts = 365*3600*24; % operated at 100% power for 1 year

P_decay = @(T) 0.066*Po*((T).^(-0.2) - (T+Ts).^(-0.2));
% --------------------- time since shutdown...time since start-up

% if no energy losses to the environment: Q_12 = integral of P_decay

Tmin = 1;
Tmax = 24*3600; % 1 day
Nt = 15;

Tspace = linspace(Tmin,Tmax,Nt);

T_vec = nan(Nt,1);
P_vec = nan(Nt,1);

for i = 1:Nt
    
    Q_total = integral(P_decay,0.0001,Tspace(i)); % get heat input
    Energy = @(T) mass_water*(water.u_Tv(T,v_water_1) - u_water_0) + ...
        mass_air*(air.u_Tv(T,v_air_1) - u_air_0) + ...
        mass_cooling*(water.u_Tv(T,v_water_1) - u_cooling_0) ...
        - Q_total;
    
    T_vec(i) = fzero(Energy,0.5*(To_containment+To_primary));
    P_vec(i) = water.P_vT(v_water_1,T_vec(i)) + ...
        air.P_vT(v_air_1,T_vec(i));
    
    
end

% plot results
hFig = figure(2);
set(hFig,'Position',[600 250 700 500]);

yyaxis left
plot(Tspace/3600,T_vec,'linewidth',2);
ylabel('Temperature ({^\circ}C)','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');

yyaxis right
plot(Tspace/3600,P_vec,'linewidth',2);
ylabel('Pressure (kPa)','fontsize',14,'fontweight','bold');

title('Pressure and Temperature vs Time','fontsize',16,'fontweight','bold');
xlabel('Time after accident (hr)','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');
grid on

% get pressure and temperature at the end of the first day

fprintf('Temperature at the end of one day = %g degrees C.\n',T_vec(end));
fprintf('Pressure at the end of one day = %g kPa.\n',P_vec(end));



