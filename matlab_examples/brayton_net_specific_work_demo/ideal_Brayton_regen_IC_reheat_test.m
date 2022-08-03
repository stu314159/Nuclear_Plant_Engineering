%% ideal Brayton Regen IC Reheat test
clear
clc
close 'all'

%% Add current directory to the Python Path
EasyProp_path = ' '; %<- Path if EasyProp.py is in your current directory
if count(py.sys.path,EasyProp_path) == 0  % <-- see if desired directory is on path
    insert(py.sys.path,int32(0),EasyProp_path); %<-- if not; add it.
end
units = 'SI';


Pmin = 101;
Tmin = 100;
Tmax = 500;
rp = 4.14286;
fluid = 'He';

gas = py.EasyProp.simpleFluid(fluid,units);

[w_net,eta_th] = idealBraytonRegenIC_Reheat(gas,Pmin,Tmin,Tmax,rp);