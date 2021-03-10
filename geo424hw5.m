%This is a script to produce a graph for problem 1b

% P-wave velocity [km/s]
alpha = 6.8;

% SV-wave velocity [km/s]
beta = 4.0;

% incidence angle in degrees between 0 and 90
i = linspace(0,90);

ra = cotd(i);
rb = sqrt((alpha^2/beta^2)*((sind(i)).^-2) - 1);

% reflection coefficient for P-wave
Rp = (4*ra.*rb - (rb.^2 - 1).^2)./(4*ra.*rb + (rb.^2 - 1).^2);

% reflection coefficient for SV-wave
Rsv = (4*ra.*(1 - rb.^2))./(4*ra.*rb + (rb.^2 - 1).^2);

% plotting
plot(i,Rp)
hold on
plot(i,Rsv)
title('Reflection Coefficients R_P and R_{SV} with respect to Angle of Incidence{\it i}')
xlabel('Angle of Incidence{\it i} [degrees]')
ylabel('Reflection Coefficient Value')
legend('R_P','R_{SV}')
