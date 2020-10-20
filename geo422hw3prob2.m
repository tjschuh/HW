function malpha = geo422hw3prob2(N,tests,bins,p,alpha)
% malpha = GEO422HW3PROB2(N,tests,bins,p,alpha)
%
% INPUT:
%
% N       user-defined set of numbers of size N
% tests   number of X2 values to compute/tests to run
% bins    number of bins to parse data into (default 10)
% p       number of independent linear constraints imposed on data (default 3)
% alpha   confidence interval (default 0.95) 
%
% OUTPUT:
%
% malpha  measured, resultant confidence interval computed by testing your test
%
% TESTED ON: 9.8.0.1451342 (R2020a) Update 5
%
% Written by tschuh@princeton.edu, 10/13/2020

% number of bins
defval('bins',10);

% Degrees of Freedom
defval('p',3);

% Confidence interval
defval('alpha',0.95);

for index = 1:tests
% Generate an Nx1 array with random numbers drawn from a standard normal distribution
R = rand(N,1);

% Calculate 0:10:100 -th percentiles of R
percents = prctile(R,0:bins:100);

% Calculate the frequencies
f = histc(R,[percents(1:end-1) percents(end) + 1]);
f = f(1:end-1);

% So now percents contains our integration limits
F = [[normcdf(percents(2:end))-normcdf(percents(1:end-1))]*N]';

% Compute X2 which is a single value statistic
% Put X2 in an array called X2mat
for i = 1:size(f)
    X2(i) = ((f(i) - F(i))^2)/F(i);
end
X2 = sum(X2);
X2mat(index) = X2;
end

% inverse cumulative distribution function (icdf)
icdf = chi2inv(alpha,bins-p);

% testing the test
% accounter is the number of X2's below confidence interval
% rejcounter is the number of X2's above the confidence interval
acccounter = 0;
rejcounter = 0;
for i = 1:size(X2mat,2)
    if icdf > X2mat(i)
       acccounter = acccounter + 1;
    else
       rejcounter = rejcounter + 1;
    end
end

% measured alpha, testing the test
malpha = acccounter/size(X2mat,2);

% plot a histogram of X2 values
histogram(X2mat,20,'Normalization','probability')
%hold on
% create a X2 distribution to compare with our data 
%chixval = 1:max(X2mat) + 10;
%chidist = chi2pdf(chixval,bins-p);
%legend1 = plot(chixval,chidist);
title('Chi Squared Histogram using Arbitrary Probability Distribution')
%legend([legend1],['DoF = ' num2str(bins-p)],'FontSize', 12)
ylabel('Frequency')
xlabel('X2')

disp(sprintf('Theoretical Confidence Interval is %d%%, Measured Confidence Interval is %g%%',100*alpha,100*malpha))
