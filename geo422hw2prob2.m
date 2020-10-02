function pphi=geo422hw2prob2(N,L,mx)
% pphi=GEO422HW2PROB2(N,L,mx)
%
% INPUT:
%
% N is the number of convolutions
% L is the number of array elements used to create the PDF
% mx is the upper bound for the random numbers chosen to create the PDF
%
% OUTPUT:
%
% pphi is some measure of Gaussiness
%
% TESTED ON: 9.8.0.1451342 (R2020a) Update 5
%
% Written by tschuh@princeton.edu, 10/01/2020

% Build a random Probability Density Function (PDF)
% always positive, includes decimal-point values
% length L, between 0 and max, 
randsz = [1,L]; 
bounds = [0,mx]; 
randpdf = rand(randsz)*range(bounds)+bounds(1); 

% Convolution of pdf1 and pdf2
convpdf = conv(randpdf,randpdf);

% Variable axis
x = linspace(-mx,mx,length(convpdf));

% Normalize by the zeroth moment
convpdf = convpdf/trapz(x,convpdf);

% Calculate the first raw moment (expectation)
epdf = trapz(x,x.*convpdf);

% Calculate the second central moment (variance)
vpdf = trapz(x,(x-epdf).^2.*convpdf);

% Normal PDF with the same expectation and variance
npdf = normpdf(x,epdf,sqrt(vpdf));

% A reasonable distance to judge fit
pphi = sum((convpdf-npdf).^2)/length(convpdf);

for index = 1:N-1
   % Successive convolutions
   convpdf = conv(convpdf,randpdf);
   x = linspace(-mx,mx,length(convpdf));
   convpdf = convpdf/trapz(x,convpdf);
   epdf = trapz(x,x.*convpdf);
   vpdf = trapz(x,(x-epdf).^2.*convpdf);
   npdf = normpdf(x,epdf,sqrt(vpdf));
   pphi = sum((convpdf-npdf).^2)/length(convpdf);
end

% End with a plot only if you didn't request output
if nargout==0
    p1 = plot(x,convpdf,'r','LineWidth',1.5);
    hold on
    % Plot a vertical bar at the expectation
    p2 = plot([epdf  epdf],ylim,'g','LineWidth',1.5);
    % Plot a horizontal bar at half height between +/-1 standard deviation
    p3 = plot(epdf+[-1 1]*sqrt(vpdf), ...
        [1 1]*max(convpdf)/2,'b','LineWidth',1.5);
    % Plot the comparative norm
    p4 = plot(x,npdf,'k','LineWidth',1.5);
    hold off
    grid on
    xlim([x(1) x(end)])
    % Title it with the measure of Gaussiness
    title(sprintf('%i convolutions, %s = %8.6f',N,'\Phi',pphi))
    xlabel('arbitrary x');
    ylabel('arbitrary y');
    legend([p1 p2 p3 p4],{'convolved pdf','expectation',...
        '1 std','normalized pdf'})
    
    % Cosmetology
    %print('-dpdf',mfilename)
end