function m=geo543hw8prob1(a,b)
% m=GEO543HW8PROB1(a,b)
%
% Plot sigmaxx and sigmayy vs x on normal scale and log-log scale
%
% INPUT:
%
% a       positive integer
% b       positive integer (b != 0 --> m != 0)
%
% OUTPUT:
%
% m       m value
%
% EXAMPLES:
%
% m=geo543hw8(10,0.01) --> m = 0.9980
% m=geo543hw8(1,0.01)  --> m = 0.9802
%
% Originally written by tschuh-at-princeton.edu, 11/02/2021

defval('a',1)
defval('b',0)

m = (a-b)/(a+b);
R = (a+b)/2;

count = 1;
for i=-4:0.25:2
    y(count) = 10^i;
    count = count + 1;
end
x1 = a*y + a;

% these are the solutions for positive x 
zeta1 = (x1 - sqrt((x1.^2) - 4*(R^2)))/(2*R);

sigxx1 = (-2.*m)./(m - real(zeta1).^(-2)) ...
         - (2.*m.*(real(zeta1).^(-4) + m.*real(zeta1).^(-2)))./((m - real(zeta1).^(-2)).^3) ...
         + ((1 + m^2).*(1 + m.*real(zeta1).^2))./(((1 - m.*real(zeta1).^2).^2).*(m - real(zeta1).^(-2)));

sigyy1 = (-2.*m)./(m - real(zeta1).^(-2)) ...
         + (2.*m.*(real(zeta1).^(-4) + m.*real(zeta1).^(-2)))./((m - real(zeta1).^(-2)).^3) ...
         - ((1 + m^2).*(1 + m.*real(zeta1).^2))./(((1 - m.*real(zeta1).^2).^2).*(m - real(zeta1).^(-2)));

x2 = logspace(log(.1*a),log(.1*(a+9)),100);

% these are the solutions for negative x
%zeta2 = (x2 + sqrt((x2.^2) - 4*(R^2)))/(2*R);

%sigxx2 = (-2.*m)./(m - zeta2.^(-2)) ...
%         - (2.*m.*(zeta2.^(-4) + m.*zeta2.^(-2)))./((m - zeta2.^(-2)).^3) ...
%         + ((1 + m^2).*(1 + m.*zeta2.^2))./(((1 - m.*zeta2.^2).^2).*(m - zeta2.^(-2)));

%sigyy2 = (-2.*m)./(m - zeta2.^(-2)) ...
%         + (2.*m.*(zeta2.^(-4) + m.*zeta2.^(-2)))./((m - zeta2.^(-2)).^3) ...
%         - ((1 + m^2).*(1 + m.*zeta2.^2))./(((1 - m.*zeta2.^2).^2).*(m - zeta2.^(-2)));

%x = [x2 x1];
%sigxx = [sigxx2 sigxx1];
%sigyy = [sigyy2 sigyy1];

% plotting
% only plotting x1 and sigxx1/sigyy1 to match plot from notes
figure
scatter(log(y),log(abs(sigxx1)),'b','filled')
hold on
scatter(log(y),log(abs(sigyy1)),'r','filled')
grid on
longticks
xlabel('log(x-a/a)')
ylabel('log(\sigma/P)')
title(sprintf('m = %d',m))