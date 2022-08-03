%% Brayton Cycle Net Specific Work Driver Script

clear
clc
close 'all'

%% Add current directory to the Python Path
EasyProp_path = ' '; %<- Path if EasyProp.py is in your current directory
if count(py.sys.path,EasyProp_path) == 0  % <-- see if desired directory is on path
    insert(py.sys.path,int32(0),EasyProp_path); %<-- if not; add it.
end
units = 'SI';

%% Analyze for different cycles and fluids
cycle_choice = 1;
% 1 = Simple Ideal Brayton
% 2 = Ideal Regenerating Brayton
% 3 = Ideal Regenerating and Intercooling Brayton
% 4 = Ideal Regenerating, Intercooling, and Reheating Brayton

switch cycle_choice
    
    case 1
        cycle_fun = @(g,Pmin,Tmin,Tmax,r) ...
            simpleIdealBrayton(g,Pmin,Tmin,Tmax,r);
        
    case 2
        cycle_fun = @(g,Pmin,Tmin,Tmax,r) ...
            idealBraytonRegen(g,Pmin,Tmin,Tmax,r);
        
    case 3
        cycle_fun = @(g,Pmin,Tmin,Tmax,r) ...
            idealBraytonRegenIC(g,Pmin,Tmin,Tmax,r);
        
    case 4
        cycle_fun = @(g,Pmin,Tmin,Tmax,r) ...
            idealBraytonRegenIC_Reheat(g,Pmin,Tmin,Tmax,r);
        
end

Pmin = 101; %kPa
Tmin = 100; %C
Tmax = 500; %C
rp_min = 1;
rp_max = 6;
N_rp = 50;
rp = linspace(rp_min,rp_max,N_rp);
numFluids = 4;


w_net = nan(N_rp,numFluids);
eta_th = nan(N_rp,numFluids);
w_net_max = nan(numFluids,1);
rp_wnet_max = nan(numFluids,1);

for f = 1:4
    fluid_choice = f;
    % 1 = Air
    % 2 = CO2
    % 3 = He
    % 4 = N2
    switch fluid_choice
        
        case 1
            fluid = 'Air';
        case 2
            fluid = 'CO2';
        case 3
            fluid = 'He';
        case 4
            fluid = 'N2';
    end
    gas = py.EasyProp.simpleFluid(fluid,units);
    
    for r = 1:N_rp
        [w_net(r,f),eta_th(r,f)] = cycle_fun(gas,Pmin,Tmin,Tmax,rp(r));
    end
    
    % get r_p for max net work and the associated max net work.
    [w_net_max(f),idx] = max(w_net(:,f));
    rp_wnet_max(f) = rp(idx);
    fprintf('Max net work for %s is %g at r_p = %g \n',fluid,w_net_max(f),...
        rp_wnet_max(f));
end

%% Plot the result
figure(1)
plot(rp,w_net(:,1),'-r',...
    rp,w_net(:,2),'-g',...
    rp,w_net(:,3),'-c',...
    rp,w_net(:,4),'-k','linewidth',3);
title_text = sprintf('Net Specific Work vs. Pressure Ratio');
title(title_text,'fontsize',16,'fontweight','bold');
xlabel('Pressure Ratio','fontsize',14,'fontweight','bold');
ylabel('Net Specific Work (kJ/kg)','fontsize',14,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
legend('Air','CO_2','He','N_2');



