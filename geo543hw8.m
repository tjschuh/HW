function m=geo543hw8(a,b)
% m=GEO543HW8(a,b)
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

%x1 = linspace(a,a+9,100);
x1 = a:0.01:a+9;

% these are the solutions for positive x 
zeta1 = (x1 - sqrt((x1.^2) - 4*(R^2)))/(2*R);

sigxx1 = (-2.*m)./(m - zeta1.^(-2)) ...
         - (2.*m.*(zeta1.^(-4) + m.*zeta1.^(-2)))./((m - zeta1.^(-2)).^3) ...
         + ((1 + m^2).*(1 + m.*zeta1.^2))./(((1 - m.*zeta1.^2).^2).*(m - zeta1.^(-2)));

sigyy1 = (-2.*m)./(m - zeta1.^(-2)) ...
         + (2.*m.*(zeta1.^(-4) + m.*zeta1.^(-2)))./((m - zeta1.^(-2)).^3) ...
         - ((1 + m^2).*(1 + m.*zeta1.^2))./(((1 - m.*zeta1.^2).^2).*(m - zeta1.^(-2)));

%x2 = linspace(-a,-a-9,100);
x2 = -a:-0.01:-(a+9);

% these are the solutions for negative x
zeta2 = (x2 + sqrt((x2.^2) - 4*(R^2)))/(2*R);

sigxx2 = (-2.*m)./(m - zeta2.^(-2)) ...
         - (2.*m.*(zeta2.^(-4) + m.*zeta2.^(-2)))./((m - zeta2.^(-2)).^3) ...
         + ((1 + m^2).*(1 + m.*zeta2.^2))./(((1 - m.*zeta2.^2).^2).*(m - zeta2.^(-2)));

sigyy2 = (-2.*m)./(m - zeta2.^(-2)) ...
         + (2.*m.*(zeta2.^(-4) + m.*zeta2.^(-2)))./((m - zeta2.^(-2)).^3) ...
         - ((1 + m^2).*(1 + m.*zeta2.^2))./(((1 - m.*zeta2.^2).^2).*(m - zeta2.^(-2)));

x = [x2 x1];
sigxx = [sigxx2 sigxx1];
sigyy = [sigyy2 sigyy1];

% plotting
subplot(2,2,1)
scatter(x,sigxx,'filled','r')
grid on
longticks
xlabel('x')
ylabel('sigma_{xx}')
subplot(2,2,2)
scatter(x,sigyy,'filled','b')
grid on
longticks
xlabel('x')
ylabel('sigma_{yy}')

% only plotting x1 and sigxx1/sigyy1 to match plot from notes
subplot(2,2,3)
scatter(log((x1-a)/a),log(sigxx1),'filled')
hold on
plot(log((x1-a)/a),log(sigxx1))
grid on
longticks
xlabel('log(x-a/a)')
ylabel('log(sigma_{xx})')
subplot(2,2,4)
scatter(log((x1-a)/a),log(sigyy1),'filled')
hold on
plot(log((x1-a)/a),log(sigyy1))
grid on
longticks
xlabel('log(x-a/a)')
ylabel('log(sigma_{yy})')