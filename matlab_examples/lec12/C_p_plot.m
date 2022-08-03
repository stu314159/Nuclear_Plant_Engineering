%% C_p plot
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

%%
Pmax = 20e3; % kPa
Pmin = 7700; % kPa

N = 500;
Th = linspace(34,424,N);
Tc = linspace(60,390,N);


Cp_hp = nan(1,N);
Cp_lp = nan(1,N);

for i = 1:N
    Cp_hp(i) = gas.Cp_pT(Pmax,Th(i));
    Cp_lp(i) = gas.Cp_pT(Pmin,Tc(i));
end

figure (1)
plot(Th,Cp_hp,'-b',...
    Tc,Cp_lp,'-r','linewidth',3);
grid on
title('Variation of C_p with Temperature and Pressure')
xlabel('T [^{\circ}C]','fontsize',14,'fontweight','bold');
ylabel('C_p [kJ/kg-^{\circ}K]','FontSize',14,'FontWeight','bold');
set(gca,'fontsize',12,'fontweight','bold')
legend('P = 20MPa','P=7.7MPa');