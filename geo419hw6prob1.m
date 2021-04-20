% Script to plot the requested P-T curves requested in GEO 419 HW 6 Problem 1b

% Originally written by Terance Schuh, 4/19/2021

% load-in Brown and Shankland mantle P-T profile
load PTdata

% ambient-pressure properties:
% melting temperature [K]
T0 = 2000;
% density [g/cm^3]
rho0 = 4.1;
% bulk modulus [GPa]
K0 = 260;

% density range, 4.1 - 5 [g/cm^3]
rho = linspace(4.1,5.485,24);

% pressure equation
P = (3/2)*K0*((rho./rho0).^(7/3) - (rho./rho0).^(5/3));

% temperature equation for gamma0 = 1
T1 = T0.*(rho0./rho).^(2/3).*exp(2.*1.*(1 - (rho0./rho)));

% temperature equation for gamma0 = 1.5
T2 = T0.*(rho0./rho).^(2/3).*exp(2.*1.5.*(1 - (rho0./rho)));

% plotting
plot(T1,P,'LineWidth',1.5)
hold on
plot(T2,P,'LineWidth',1.5)
hold on
plot(PTdata(:,3),PTdata(:,2),'LineWidth',1.5)
hold on
yline(135.7,'LineWidth',1)
title('P-T Diagram for (Mg,Fe)SiO_3 Perovskite')
xlabel('Temperature [K]')
ylabel('Pressure [GPa]')
legend({'\gamma_0 = 1','\gamma_0 = 1.5',...
        ['Brown & Shankland' char(10) 'Mantle P-T Data'],...
        ['Pressure at Base' char(10) 'of Lower Mantle']},'Location','southeast')
grid on