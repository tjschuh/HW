w = pi/4; %angular frequency
Nrate = pi/w; %Nyquist rate
deltax = 1/Nrate; %max sampling rate, computed from Nyquist rate

%plot the true function using high precision
n = 0:0.001:99;
x = 7*sin(w*n);
plot(n,x)
hold on

%plot using the Nyquist rate
srate = deltax + .4; %Nrate goes as 1/srate
n = 0:srate:99;
x = 7*sin(w*n);
plot(n,x)
hold on

%plot using a larger Nyquist rate
srate = deltax - .24; %Nrate goes as 1/srate
n = 0:srate:99;
x = 7*sin(w*n);
plot(n,x)
hold on

%plot using a smaller Nyquist rate
%srate = deltax + 1; %Nrate goes as 1/srate
%n = 0:srate:19;
%x = 7*sin(w*n);
%plot(n,x)