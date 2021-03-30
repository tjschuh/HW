% This is a script to plot the graph requested in HW 4
% Written by tschuh@princeton.edu, 3/29/2021

% Problem 3

% Pressures in lower mantle with increasing depth [GPA]
P = [23.8 28.3 32.8 37.3 41.9 46.5 51.2 55.9 60.7 65.5 70.4 ...
     75.4 80.4 85.5 90.6 95.8 101.1 106.4 111.9 117.4 123 127 128.8 ...
     134.6 135.8];
% Densities in lower mantle with increasing depth [g/cm^3]
rho = [4.38 4.44 4.5 4.56 4.62 4.68 4.73 4.79 4.84 4.9 4.95 ...
       5.00 5.05 5.11 5.16 5.21 5.26 5.31 5.36 5.41 5.46 5.49 5.51 ...
       5.56 5.57];
rho0 = 4;

F = (2.*P)./(3.*((rho./rho0).^(7/3) - (rho./rho0).^(5/3)));
f = (1/2).*((rho./rho0).^(2/3) - 1);

scatter(f,F)
title('Linearized Birch-Murnaghan Equation','fontsize',14)
xlabel('${\it} f = \frac{1}{2}[(\frac{\rho}{\rho_0})^{\frac{2}{3}} - 1]$',...
    'Interpreter','Latex','fontsize',14)
ylabel('${\it} F = \frac{P}{\frac{3}{2}[(\frac{\rho}{\rho_0})^{\frac{7}{3}} - (\frac{\rho}{\rho_0})^{\frac{5}{3}}]}$',...
    'Interpreter','Latex','fontsize',14)
grid on
hold on

% Least Squares Fit
ls = polyfit(f,F,1); 
xfit = linspace(0,0.14,100); 
yfit = polyval(ls,xfit);
plot(xfit,yfit)
hold off
m = ls(1)
yinter = ls(2)