function pphi=geo422hw01(N)
% pphi=GEO422HW01(N)
%
% INPUT:
%
% N    is the number of convolutions
% 
% OUTPUT:
%
% pphi  is some measure of Gaussiness
%
% EXAMPLE:
% 
% for index=1:4; subplot(2,2,index); geo422hw01(index); end
%
% TESTED ON: 9.8.0.1451342 (R2020a) Update 5
%
% Written by fjsimons@princeton.edu, 09/29/2020

% This is a script
% This is not a function

% Ceci n'est pas un pdf
pdf1=[0 0 0 0 1 2 3 4 5 0 0 0 0];
pdf2=[0 0 0 0 pi pi log(2) log(2) 23 0 0 0 0];

% This is going to be the distribution of 
% the sum random variables from those two pdfs
pdff=conv(pdf1,pdf2);
% Make a variable axis
x=linspace(-4,4,length(pdff));
% Normalize by the zeroth moment
pdff=pdff/trapz(x,pdff);
% Calculate the first raw moment
epdf=trapz(x,x.*pdff);
% Calculet the second central moment
vpdf=trapz(x,(x-epdf).^2.*pdff);
% Let's generate a normal with the same expectation and variance
cpdf=normpdf(x,epdf,sqrt(vpdf));
% Invent some kind of reasonable distance to judge fit
pphi=sum((pdff-cpdf).^2)/length(pdff);
for index=1:N-1
    % Here is the successive convolution with the second one
   pdff=conv(pdff,pdf2);
   x=linspace(-4,4,length(pdff));
   pdff=pdff/trapz(x,pdff);
   epdf=trapz(x,x.*pdff);
   vpdf=trapz(x,(x-epdf).^2.*pdff);
   cpdf=normpdf(x,epdf,sqrt(vpdf));
   pphi=sum((pdff-cpdf).^2)/length(pdff);
end

% End with a plot only if you didn't request output
if nargout==0
    plot(x,pdff)
    hold on
    % Plot a vertical bar at the expectation
    plot([epdf  epdf],ylim)
    % Plot a horizontal bar at half height between +/1 standard deviation
    plot(epdf+[-1 1]*sqrt(vpdf),[1 1]*max(pdff)/2)
    % Plot the comparative norma
    plot(x,cpdf,'k')
    hold off
    grid on
    % Title it with the measure of Gaussiness
    title(sprintf('%i convolutions, %s = %8.6f',N,'\Phi',pphi))
    
    % Cosmetology
    print('-dpdf',mfilename)
end




