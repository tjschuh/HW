% This is a script to plot the graph requested in HW 5
% Written by tschuh@princeton.edu, 4/1/2021

% Problem 2

% core density and surface density [kg/m^3]
rhoc = 13000;
rho0 = 3000;

% gravitational constant [N*m^2/kg^2]
G = 6.67e-11;

% radius [m]
R = 6371e3;

r = linspace(0,R,10000);

Pc = 2*pi*G*(((rhoc^2)/3) - r.*((7*rhoc*(rhoc-rho0))/(18*R)) + ...
    (r.^2).*(((rhoc-rho0)^2))/(8*(R^2))).*(r.^2);

g = 4*pi*G*((rhoc/3) - r.*((rhoc-rho0)/(4*R))).*r;

plot(r,Pc)
title('Pressure vs. Depth Below Earth''s Surface')
xlabel('Depth [km]')
ylabel('Pressure [GPa]')
xlim([0 R])
xticklabels({'0','1000','2000','3000','4000','5000','6000'})
yticklabels({'0','50','100','150','200','250','300','350'})
grid on

figure
plot(r,g)
title('Acceleration Due To Gravity vs Depth Below Earth''s Surface')
xlabel('Depth [km]')
ylabel('Acceleration due to Gravity [m/s^2]')
xlim([0 R])
xticks([0.371e6 1.371e6 2.371e6 3.371e6 4.371e6 5.371e6 R])
xticklabels({'6000','5000','4000','3000','2000','1000','0'})
grid on
set(gca, 'XDir','reverse')
