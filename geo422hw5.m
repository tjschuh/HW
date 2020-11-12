function T = geo422hw5()
% T = GEO422HW5()
%
% INPUT:
%
% OUTPUT:
%
% T  forward model time prediction
%
% TESTED ON: 9.8.0.1451342 (R2020a) Update 5
%
% Written by tschuh@princeton.edu, 11/11/2020

load geiger_student.mat %contains mediumvelocity and stationlocations
v = mediumvelocity * 10^3;
S = [stationlocations noisyarrivaltimes];
M = [1 2 3 4]; %initial earthquake location guess [x y z t]
T = forward(S,M,v);

% plotting
subplot(2,2,1)
plot(S(:,1),S(:,2),'^')
title('XY positioning of stations')
xlabel('X')
ylabel('Y')
grid on

subplot(2,2,2)
plot(S(:,1),S(:,3),'^')
title('XZ positioning of stations')
xlabel('X')
ylabel('Z')
grid on

subplot(2,2,3)
plot(S(:,2),S(:,3),'^')
title('YZ positioning of stations')
xlabel('Y')
ylabel('Z')
grid on

figure
plot3(S(:,1),S(:,2),S(:,3),'^')

function T = forward(S,M,v)
T = M(4) + sqrt((S(:,1) - M(1)).^2 + (S(:,2) - M(2)).^2 + (S(:,3) - M(3)).^2)/v;