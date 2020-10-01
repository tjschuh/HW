function pphi=geo422hw2prob1(N)
% pphi=GEO422HW2PROB1(N)
%
% INPUT:
%
% N is the number of convolutions 
%
% OUTPUT:
%
% pphi is some measure of Gaussiness
%
% EXAMPLE:
% 
% for index=1:4; subplot(2,2,index); geo422hw2prob1(index); end
%
% TESTED ON: 9.8.0.1451342 (R2020a) Update 5
%
% Written by tschuh@princeton.edu, 09/30/2020

% Probability Density Functions (PDFs)
pdf1=[0 0 0 0 5 4 3 2 1 0 0 0 0];
pdf2=[0 0 0 0 pi pi log(3) log(2) sqrt(5) 0 0 0 0];

% Convolution of pdf1 and pdf2
convpdf=conv(pdf1,pdf2);

% Variable axis
x=linspace(-4,4,length(convpdf));

% Normalize by the zeroth moment
convpdf=convpdf/trapz(x,convpdf);

% Calculate the first raw moment (expectation)
epdf=trapz(x,x.*convpdf);

% Calculate the second central moment (variance)
vpdf=trapz(x,(x-epdf).^2.*convpdf);

% Normal PDF with the same expectation and variance
npdf=normpdf(x,epdf,sqrt(vpdf));

% A reasonable distance to judge fit
pphi=sum((convpdf-npdf).^2)/length(convpdf);

for index=1:N-1
   % Here is the successive convolution with the second one
   convpdf=conv(convpdf,pdf2);
   x=linspace(-4,4,length(convpdf));
   convpdf=convpdf/trapz(x,convpdf);
   epdf=trapz(x,x.*convpdf);
   vpdf=trapz(x,(x-epdf).^2.*convpdf);
   npdf=normpdf(x,epdf,sqrt(vpdf));
   pphi=sum((convpdf-npdf).^2)/length(convpdf);
end

% End with a plot only if you didn't request output
if nargout==0
    p1 = plot(x,convpdf,'r','LineWidth',1.5)
    hold on
    % Plot a vertical bar at the expectation
    p2 = plot([epdf  epdf],ylim,'g','LineWidth',1.5)
    % Plot a horizontal bar at half height between +/-1 standard deviation
    p3 = plot(epdf+[-1 1]*sqrt(vpdf),[1 1]*max(convpdf)/2,'b','LineWidth',1.5)
    % Plot the comparative norm
    p4 = plot(x,npdf,'k','LineWidth',1.5)
    hold off
    grid on
    % Title it with the measure of Gaussiness
    title(sprintf('%i convolutions, %s = %8.6f',N,'\Phi',pphi))
    legend([p1 p2 p3 p4],{'convolved pdf','expectation value',...
        '1 standard deviation','normalized pdf'})
    
    % Cosmetology
    print('-dpdf',mfilename)
end