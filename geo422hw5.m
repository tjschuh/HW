function M = geo422hw5(x,y,z,t)
% M = GEO422HW5(x,y,z,t)
%
% INPUT:
%
% x   initial earthquake location x-coordinate guess (can be +/-)
% y   initial earthquake location y-coordinate guess (can be +/-)
% z   initial earthquake location z-coordinate guess (can be +/-)
% t   initial earthquake location time-coordinate guess (can be +/-)
%
% OUTPUT:
%
% M   x,y,z,t coordinates of final calculated location of earthquake
%
% TESTED ON: 9.8.0.1451342 (R2020a) Update 5
%
% Written by tschuh@princeton.edu, 11/11/2020

load geiger_student.mat %contains mediumvelocity and stationlocations
v = mediumvelocity; %in km/sec
S = stationlocations; %in km
T = noisyarrivaltimes; %observed time
M = [x; y; z; t]; %initial earthquake location guess

%F = misfit(T,t); %use this somehow to know how many iterations to do

for i = 1:100
t = forward(S,M,v); %predicted time (fancy g in notes)
K = sensitivity(S,M,v,t);
dM = invert(K,T,t);
M = M + dM;
end

% plotting
subplot(2,2,1)
scatter(S(:,1),S(:,2),'^','filled')
hold on
scatter(M(1),M(2),'*','LineWidth',1.25)
title('XY 2D positioning of stations')
xlabel('X [km]')
ylabel('Y [km]')
grid on

subplot(2,2,2)
scatter(S(:,1),S(:,3),'^','filled')
hold on
scatter(M(1),M(3),'*','LineWidth',1.25)
title('XZ 2D positioning of stations')
xlabel('X [km]')
ylabel('Z [km]')
grid on

subplot(2,2,3)
scatter(S(:,2),S(:,3),'^','filled')
hold on
scatter(M(2),M(3),'*','LineWidth',1.25)
title('YZ 2D positioning of stations')
xlabel('Y [km]')
ylabel('Z [km]')
grid on

figure
scatter3(S(:,1),S(:,2),S(:,3),'^','filled')
hold on
scatter3(M(1),M(2),M(3),'*','LineWidth',1.25)
title('XYZ 3D positioning of stations')
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
grid on

function t = forward(S,M,v)
t = M(4) + sqrt((S(:,1) - M(1)).^2 + (S(:,2) - M(2)).^2 + (S(:,3) - M(3)).^2)/v;

function F = misfit(T,t)
%sums of sqaures of observed - predicted

function K = sensitivity(S,M,v,t)
%partial derivatives from handout
for j = 1:size(S,1)
partial1 = -v^(-2) * ((S(j,1)-M(1))/(t(j)-M(4)));
partial2 = -v^(-2) * ((S(j,2)-M(2))/(t(j)-M(4)));
partial3 = -v^(-2) * ((S(j,3)-M(3))/(t(j)-M(4)));
partial4 = 1;
K(j,:) = [partial1 partial2 partial3 partial4];
end

function dM = invert(K,T,t)
%look at notes from class
%take generalized inverse
%equation 13 in notes, deltaD is epsilon = observed - predicted
deltad = T - t;
dM = inv(transpose(K)*K)*transpose(K)*deltad;