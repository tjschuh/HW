% TESTED ON: 9.8.0.1451342 (R2020a) Update 5
%
% Written by tschuh@princeton.edu, 10/31/2020

% Defining some parameters

%1 frame = 0.04 seconds
frame = 0.04;
%1 pixel = 0.2 cm = 0.002 m
pix = .002;
%diameter of golfball = 4 cm, 20 pixels
%frame is 128 cm by 96 cm
%first and last images we are looking at
firstfile = 31;
lastfile = 54;
%how many times are we running the experiment
totexp = 5;

% PART 1
%if we already did the experiments, then we dont repeat them
if exist('ycom','var') == 0 && exist('xcom','var') == 0 && exist('time','var') == 0
for exp = 1:totexp
   counter = 1;
   for index = firstfile:lastfile
   image(imread(sprintf('http://geoweb.princeton.edu/people/simons/GOLFBALL/000000%2.2i.jpg',index)))
   [X,Y] = ginput;
   xcom(counter,exp) = X;
   ycom(counter,exp) = 480 - Y;
   time(counter,exp) = (counter - 1);
   counter = counter + 1;
   end
end
   % putting displacements and time in meters and seconds, respectively
   xcom = pix.*xcom;
   ycom = pix.*ycom;
   time = frame.*time;
else
   %xcom and ycom already exist, load them in
end

%this is just to visualize the data
plot(xcom(:,1)',ycom(:,1)')
xlabel('Displacement (m)')
ylabel('Displacement (m)')
title('Sample Trajectory of Golfball')

% PART 2
Y = ycom;
G = [ones(size(time,1),1) time(:,1) ((time(:,1)).^2)./-2];
for exp = 1:totexp
mhat(:,exp) = inv(transpose(G) * G) * transpose(G) * Y(:,exp);
end

% PART 3
%assuming our variance is a constant
sigma = 0.5;
%calculate covariance matrices
Cy = (sigma)^2 * eye(size(Y,1));
Gminusg = inv(transpose(G) * G) * transpose(G);
Cm = Gminusg * Cy * transpose(Gminusg);

% PART 4
%calculate a new mhat using covariance matrix Cm
mhatprime = inv(transpose(G) * inv(Cy) * G) ...
            * transpose(G) * inv(Cy) * Y;
        
% PART 5
%perform a chi squared test for the significance of our fits
for exp = 1:totexp
    X2obs(:,exp) = (Y(:,exp) - (G*mhatprime(:,exp))).^2/(sigma^2); 
end

%# of model parameters is 3
modparam = 3;
%calculate the # of degrees of freedom
DoF = size(Y,1) - modparam;
alpha = 0.95;
cutoff = chi2inv(alpha,DoF);

for exp = 1:totexp
acccounter = 0;
for j = 1:size(X2obs,1)
    if cutoff > X2obs(j,exp)
       acccounter = acccounter + 1;
    else
       %this means it was below icdf threshold
    end
end
%calculate percentage of values accepted per experiment 
accpercent(exp) = acccounter/size(X2obs,1); 
end

% PART 6
%construct confidence intervals around my estimates
for i = 1:size(mhat,1)
    mhatmean(i,1) = mean(mhat(i,:));
    mhatstd(i,1) = std(mhat(i,:));
end

conint = sprintf('m1 = %.4f +/- %.4f, \n m2 = %.4f +/- %.4f, \n m3 = %.4f +/- %.4f',...
      mhatmean(1),mhatstd(1),mhatmean(2),mhatstd(2),mhatmean(3),mhatstd(3))