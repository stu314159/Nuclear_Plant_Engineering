%% de Stordeur's Correlation Example

clear;
clc;
close 'all'

%% Prepare the environment for EasyProp

EasyProp_path = ' '; %<- Path if EasyProp.py is in your current directory
if count(py.sys.path,EasyProp_path) == 0  % <-- see if desired directory is on path
    insert(py.sys.path,int32(0),EasyProp_path); %<-- if not; add it.
end

%% Establish working fluid and unit system
fluid = 'Air';
unitSystem = 'SI';

water = py.EasyProp.simpleFluid('Water',unitSystem);

%% Problem Parameters
Press = 15.51e3; % kPa, nominal primary pressure
T_ave = 0.5*(281 + 321); % C, average coolant temperature
m_dot = 51.4e6*(1/3600); % kg/s, core mass flow rate
n_channels = 157*289; %157 assemblies, 289 rods per assembly
m_dot_c = m_dot/n_channels; % kg/s, channel mass flow rate

mu = water.mu_pT(Press,T_ave); % N-s/m^2, nominal coolant viscosity
rho = (water.v_pT(Press,T_ave))^(-1); % kg/m^3, nominal coolant density

% geometric parameters
P = 12.60e-3; % m, pin pitch
D = 9.5e-3; % m, pin diameter

t = 4.5e-4; % m, grid strap thickness
H = 4e-2; % m, grid strap height

%% Calculations

Av = P^2 - (pi/4)*D^2; % m^2, flow area in channel
As = 2*P*t - t^2; % m^2, projected area of spacer

v_avg = m_dot_c/(rho*Av); % m/s, avg velocity in channel away from spacer
v_s = v_avg*(Av/(Av - As)); % m/s, velocity in spacer region.

Re_s = t*v_s*rho/mu;
fprintf('Re_s = %g \n',Re_s);

%%
% 
%  Here the user of the correlation needs to read the
%  graph to obtain C_s
% 

Cs = 1.7; % approximate drag coefficient based on Re_s.  
% observe that for high Re_s, the drag coefficient 
% settles to around 1.7

dP = Cs*(0.5*rho*v_s^2)*(As/Av);
fprintf('dP using de Stordeur: %g Pa\n',dP);


