% This is a script for GEO 424 HW 9

% Originally written by Terance Schuh, 4/12/2021

% load in ISC data that is given
load ISC_data

% load in jeffrey's model that is given
load jmodel

% distance in degrees up to 24
deltad = ISC_data(:,1);

% convert delta to radians
deltar = (pi/180)*deltad; 

% travel times in seconds
T = ISC_data(:,2);

% ISC ray parameter given
iscp = ISC_data(:,3);

% Parts 1 & 2

% computed ray parameter p, p = dT/d(delta)
p = diff(T)./diff(deltar);

figure
plot(deltar,T,'LineWidth',1.5)
title('Travel Time{\it T} as a function of Distance{\it \Delta}')
xlabel('Distance{\it \Delta} [radians]')
ylabel('Travel Time{\it T} [s]')
grid on

figure
plot(deltar(2:end),p,'LineWidth',1.5)
hold on
plot(deltar,iscp,'LineWidth',1.5)
title('Ray Parameter{\it p} as a function of Distance{\it \Delta}')
xlim([-0.1 1.8])
xlabel('Distance{\it \Delta} [radians]')
ylabel('Ray Parameter{\it p}')
legend('Calculated Ray Paramter','ISC Ray Parameter')
grid on

% Parts 3, 4, & 5
  
% radius of earth [km]
r0 = 6371;

dxd = 2;
dxr = 2*(pi/180);
eta = p;
counter = 2;

% special case for L = 0
r(1,1) = r0;

% calculating radius using p and the trapezoidal integral method
for L = 2:dxd:96
  N = (L/dxd) + 1;
  for j = 2:N
    f1 = log(p(j-1,1)/eta(N,1) + ((p(j-1,1)/eta(N,1))^2 - 1)^(1/2));
    f2 = log(p(j,1)/eta(N,1) + ((p(j,1)/eta(N,1))^2 - 1)^(1/2));
    trap(j-1) = ((f1+f2)*dxr)/2;
  end
  area = sum(trap);
  r(counter,1) = r0*exp((-1/pi)*area);
  counter = counter + 1;
end

% velocity v [km/s] --> r [km], eta [s]
v = r./eta;

figure
plot(v,r,'LineWidth',1.5)
hold on
plot(jmodel(1:22,2),r0-jmodel(1:22,1),'LineWidth',1.5)
title('P-Wave Velocity{\it v_p} as a function of Radius{\it r}')
ylim([3400 6500])
xlabel('P-Wave Velocity{\it v_p} [km/s]')
ylabel('Radius{\it r} [km]')
legend('Calculated P-Wave Velocity','Jeffrey''s Model')
grid on
