% This is a script to plot the graphs requested in HW 7
% Written by tschuh@princeton.edu, 3/22/2021

% Problem 1

% period [s]
T = [20 30 40 50];

% group velocity [km/s] computed
% from wave train arrivals given 
U = [3.15 3.48 3.68 3.89];

scatter(T,U,'filled','^','r')
hold on
plot(T,U,'k')
title('U-T Diagram between PAS and NEE Seismic Stations')
xlabel('Period{\it T} [s]')
ylabel('Group Velocity{\it U} [km/s]')
ylim([3.0 4.0])
grid on

% Problem 2

% period [s]
T = [20.08 30.12 40.96 51.20];

% phase velocity [km/s] computed
% from pairs of harmonic waves given
c = [3.62 3.75 3.85 3.98];

figure
scatter(T,c,'filled','^','r')
hold on
plot(T,c,'k')
title('c-T Diagram between PAS and NEE Seismic Stations')
xlabel('Period{\it T} [s]')
ylabel('Phase Velocity{\it c} [km/s]')
xlim([20 T(4)])
grid on
