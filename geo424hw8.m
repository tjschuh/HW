function geo424hw8(x,y,z,t,numit)
% GEO424HW8(x,y,z,t,numit)
%
% INPUT:
%
% x      initial earthquake location x-coordinate guess in km (can be +/-)
% y      initial earthquake location y-coordinate guess in km (can be +/-)
% z      initial earthquake location z-coordinate guess in km (can be +/-)
% t      initial earthquake location t-coordinate guess in s (can be +/-)
% numit  desired number of iterations, usually 6 is enough
%
% TESTED ON: 9.8.0.1451342 (R2020a) Update 5
%
% Written by tschuh@princeton.edu, 4/6/2021

% initial earthquake guess [x y z t], space in km, time in s
M = [x; y; z; t];

% load in station data
load stationdata.mat
% station locations [x y z] [km]
S = stationdata(:,1:3);
% arrival times [t] [s]
T = stationdata(:,4);
% medium velocity [km/s]
v = 6;

iteration = 1

%calculate everything once to get a baseline result
t = forward(S,M,v); %predicted time
K = sensitivity(S,M,v,t);
dM = invert(K,T,t);
M = M + dM

%calculate a new location M numit-1 more times
for i = 1:numit-1
    iteration = iteration + 1
    t = forward(S,M,v);
    K = sensitivity(S,M,v,t);
    dM = invert(K,T,t);
    M = M + dM
end

% plotting
figure
subplot(2,2,1)
scatter(S(:,1),S(:,2),'^','filled')
hold on
scatter(M(1),M(2),'*','LineWidth',1.25)
title('XY 2D Positioning Of Stations')
xlabel('X [km]')
ylabel('Y [km]')
grid on

subplot(2,2,2)
scatter(S(:,1),S(:,3),'^','filled')
hold on
scatter(M(1),M(3),'*','LineWidth',1.25)
title('XZ 2D Positioning Of Stations')
xlabel('X [km]')
ylabel('Z [km]')
grid on

subplot(2,2,3)
scatter(S(:,2),S(:,3),'^','filled')
hold on
scatter(M(2),M(3),'*','LineWidth',1.25)
title('YZ 2D Positioning Of Stations')
xlabel('Y [km]')
ylabel('Z [km]')
grid on

subplot(2,2,4)
scatter3(S(:,1),S(:,2),S(:,3),'^','filled')
hold on
scatter3(M(1),M(2),M(3),'*','LineWidth',1.25)
title('XYZ 3D Positioning Of Stations')
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
grid on

%forward model to calculate predicted travel time t
function t = forward(S,M,v)
t = M(4) + sqrt((S(:,1) - M(1)).^2 + (S(:,2) - M(2)).^2 ...
    + (S(:,3) - M(3)).^2)/v;

%sensitivity matrix K is essentially the gradient of the forward model
function K = sensitivity(S,M,v,t)
for j = 1:size(S,1)
    partial1 = -v^(-2) * ((S(j,1)-M(1))/(t(j)-M(4)));
    partial2 = -v^(-2) * ((S(j,2)-M(2))/(t(j)-M(4)));
    partial3 = -v^(-2) * ((S(j,3)-M(3))/(t(j)-M(4)));
    partial4 = 1;
    K(j,:) = [partial1 partial2 partial3 partial4];
end

%compute delta m using K and delta d
function dM = invert(K,T,t)
deltad = T - t; %residuals (observed - predicted)
dM = inv(transpose(K)*K)*transpose(K)*deltad;