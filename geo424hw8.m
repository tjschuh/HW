function [M,iteration] = geo424hw8(x,y,z,t)
% [M,iteration] = GEO422HW5(x,y,z,t)

% Written by tschuh@princeton.edu, 4/6/2021

% initial earthquake guess [x y z t]
M = [x; y; z; t];

% load in station data
load stationdata.mat
% station locations [x y z] [km]
S = stationdata(:,1:3);
% arrival times [t] [s]
T = stationdata(:,4);
% medium velocity [km/s]
v = 6;

% number of iterations
numit = 10;
iteration = 1;

%calculate everything once to get a baseline result
t = forward(S,M,v); %predicted time (fancy g in notes)
K = sensitivity(S,M,v,t);
[dM,deltad] = invert(K,T,t);
phi(iteration) = misfit(deltad); %store each phi into an array
M = M + dM;

%calculate a new location M numit-1 more times
for i = 1:numit-1
    iteration = iteration + 1;
    t = forward(S,M,v);
    K = sensitivity(S,M,v,t);
    [dM,deltad] = invert(K,T,t);
    phi(iteration) = misfit(deltad);
    M = M + dM;
end

% plotting
figure
subplot(2,3,1)
scatter(S(:,1),S(:,2),'^','filled')
hold on
scatter(M(1),M(2),'*','LineWidth',1.25)
title('XY 2D Positioning Of Stations')
xlabel('X [km]')
ylabel('Y [km]')
grid on

subplot(2,3,2)
scatter(S(:,1),S(:,3),'^','filled')
hold on
scatter(M(1),M(3),'*','LineWidth',1.25)
title('XZ 2D Positioning Of Stations')
xlabel('X [km]')
ylabel('Z [km]')
grid on

subplot(2,3,3)
scatter(S(:,2),S(:,3),'^','filled')
hold on
scatter(M(2),M(3),'*','LineWidth',1.25)
title('YZ 2D Positioning Of Stations')
xlabel('Y [km]')
ylabel('Z [km]')
grid on

subplot(2,3,4)
scatter3(S(:,1),S(:,2),S(:,3),'^','filled')
hold on
scatter3(M(1),M(2),M(3),'*','LineWidth',1.25)
title('XYZ 3D Positioning Of Stations')
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
grid on

subplot(2,3,5)
semilogy(phi,'-o','MarkerSize',8)
title('Misfit Evolving Over Iterations')
xlabel('Iteration')
ylabel('Misfit')
xlim([0.5 size(phi,2)+0.5])
ylim([0.1*min(phi) 10*max(phi)])
yticks([.1 1 10 100 1000 10000 100000 1000000])
grid on

subplot(2,3,6)
stem(transpose(deltad))
title({'Current Residuals Between Observed','And Predicted Travel Times'})
xlabel('Station Number')
ylabel('Residual [s]')
xlim([0 size(S,1)+1])
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
function [dM,deltad] = invert(K,T,t)
deltad = T - t; %residuals (observed - predicted)
dM = inv(transpose(K)*K)*transpose(K)*deltad;

%phi tells us how "good" our current prediction is
function phi = misfit(deltad)
phi = (sum(deltad))^2; %sums of squares of observed - predicted