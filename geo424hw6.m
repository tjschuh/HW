% This is a script to plot the graphs requested in HW 6

% height [km]
H = 35;

% shear moduli [dyne/cm^2]
u1 = 3e11;
u2 = 7e11;

% s-wave speeds [km/s]
B1 = 3.5;
B2 = 4.6;

% phase velocity [km/s]
% B1 < c < B2
n = 1100;
c = linspace(B1,B2,n);

% period [s]
T = [20 35 50 80];

S1 = sqrt((c.^2./B1^2) - 1);
S2 = sqrt(1 - (c.^2./B2^2));

RHS = (u2.*S2)./(u1.*S1);

for i = 1:size(T,2)
  LHS = tan((2*pi*H.*S1)./(c.*T(i)));

  for j = 1:size(LHS,2)
    if LHS(j) >= 150 || LHS(j) <= -150
       LHS(j) = NaN;
    end
  end

  % find specific points to highlight
  diff = LHS - RHS;
  if i == 1
    intersect(i) = find(diff<0.011 & diff>-0.011);
  else
    intersect(i) = find(diff<0.0025 & diff>-0.0025);
  end
  c1 = find(c<3.5201 & c>3.5199);
  c2 = find(c<3.6001 & c>3.5999);
  c3 = find(c<3.8004 & c>3.8002);
  c4 = find(c<4.0006 & c>4.0004);
  c5 = find(c<4.1997 & c>4.1995);
  c6 = find(c<4.5001 & c>4.4999);

  subplot(2,2,i)
  plot(c,RHS,'b')
  hold on
  plot(c,LHS,'r')
  hold on
  scatter(c(intersect(i)),RHS(intersect(i)),75,'*','k')
  hold on
  scatter([c(c1) c(c2) c(c3) c(c4) c(c5) c(c6)],[RHS(c1) RHS(c2) RHS(c3) RHS(c4) RHS(c5) RHS(c6)],'filled','^','b')
  hold on
  scatter([c(c1) c(c2) c(c3) c(c4) c(c5) c(c6)],[LHS(c1) LHS(c2) LHS(c3) LHS(c4) LHS(c5) LHS(c6)],'filled','^','r')
  hold on
  if i == 1 %draw asymptote
    xline(4.0415,'--k');
  end
  grid on
  title({['RHS and LHS values as a fucntion of phase velocity{\it c}'];['T = ' num2str(T(i)) ' sec']})
  legend('RHS','LHS')
  xlabel('phase velocity{\it c} [km/s]')
  ylabel('RHS and LHS values')
  xlim([B1-0.1 B2+0.1])
  if i == 1
    ylim([-17.5 22.5])
  else
    ylim([-1 16])  
  end
end

% number 3

for m = 1:4      
  cT(m) = c(intersect(m));
end

figure
plot(T,cT,'k')
hold on
scatter(T,cT,'filled','o','k')
title('Love Wave Phase Velocity Dispersion')
xlabel('Period{\it T} [s]')
ylabel('phase velocity{\it c} [km/s]')
grid on
