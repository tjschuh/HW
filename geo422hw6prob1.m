function Nrate = geo422hw6prob1(w,amp,xshift,yshift)
% Nrate = GEO422HW5(w,amp,xshift,yshift)
%
% INPUT:
%
% w       angular frequency
% amp     amplitude of signal
% xshift  right/left shift of signal
% yshift  up/down shift of signal
%
% OUTPUT:
%
% Nrate   Nyquist rate computed from twice
%         the highest frequency in signal
%
% TESTED ON: 9.8.0.1451342 (R2020a) Update 5
%
% Written by tschuh@princeton.edu, 12/1/2020

%defining variables
defval('w',3*pi/8);
defval('amp',7);
defval('xshift',0);
defval('yshift',1);

%only use this section if I decide to use more complex function
%w1 = 13*pi/16;
%w2 = 15*pi/16;
%if w1 > w2
%    w = w1;
%elseif w1 < w2
%    w = w2;
%else
%    w = w1;
%end

%other variables needed
xmax = 19; %upper bound of samples
hpsrate = 0.1; %high precision sampling rate

%Nyquist rate and min sampling rate computed from Nyquist rate
Nrate = w/pi;
deltax = 1/(2*Nrate);

%plot the true function using high precision
x = 0:hpsrate:xmax;
y = amp*sin(w*x - xshift) + yshift;
%y = sin(w1*x).*cos(w2*x);
subplot(4,2,1)
plot(x,y)
title('True Plot of Function')
subplot(4,2,2)
periodogram(y)

%plot using the Nyquist rate
srate = deltax; %Nrate goes as 1/srate
x = 0:srate:xmax;
y = amp*sin(w*x - xshift) + yshift;
%y = sin(w1*x).*cos(w2*x);
subplot(4,2,3)
plot(x,y)
title('Sampling Rate = Nyquist Rate')
subplot(4,2,4)
periodogram(y)

%plot using a larger Nyquist rate (should still work)
srate = deltax - (0.5*deltax); %Nrate goes as 1/srate
x = 0:srate:xmax;
y = amp*sin(w*x - xshift) + yshift;
%y = sin(w1*x).*cos(w2*x);
subplot(4,2,5)
plot(x,y)
title('Sampling Rate > Nyquist Rate')
subplot(4,2,6)
periodogram(y)

%plot using a smaller Nyquist rate (should blow-up)
srate = deltax + (0.5*deltax); %Nrate goes as 1/srate
x = 0:srate:xmax;
y = amp*sin(w*x - xshift) + yshift;
%y = sin(w1*x).*cos(w2*x);
subplot(4,2,7)
plot(x,y)
title('Sampling Rate < Nyquist Rate')
subplot(4,2,8)
periodogram(y)