% 1 frame = 0.04 seconds
frame = 0.04;
% 1 pixel = 0.2 cm = 0.002 m
pix = .002;
% diameter of golfball = 4 cm, 20 pixels
% frame is 128 cm by 96 cm
if exist('xcom','var') == 0
   counter = 1;
   for index = 31:54 %ideal numbers
   image(imread(sprintf('http://geoweb.princeton.edu/people/simons/GOLFBALL/000000%2.2i.jpg',index)))
   [X,Y] = ginput;
   xcom(counter) = X;
   ycom(counter) = 480 - Y;
   time(counter) = (counter - 1);
   counter = counter + 1;
   end
   % putting displacements and time in meters and seconds, respectively
   xcom = pix.*xcom;
   ycom = pix.*ycom;
   time = frame.*time;
else
   %xcom and ycom already exist, load them in
end

plot(xcom,ycom)

Y = ycom';
G = [ones(size(time,2),1) time' ((time').^2)./-2];
m = linsolve(G,Y);